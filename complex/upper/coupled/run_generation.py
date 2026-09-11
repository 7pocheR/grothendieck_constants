"""Generate and freshly replay all retained modes into an external directory."""
import argparse
import fcntl
from fractions import Fraction as Q
import os
from pathlib import Path
import signal
import sys
import time

from package_support import (ROOT, SOURCE, atomic_new, candidate, code_binding, compressed,
                             digest, environment, external, load, object_hash,
                             require, safe_path, source_binding, unoptimized)
sys.path.insert(0, str(ROOT/'validation'))
from linux_lifecycle import arm_parent_death, protected_spawn, stop_group
from full_common import conditioning_rows, inventory


def scheduled_production(row, d, settings, binding, directory):
    """Evaluate the prescribed successful radius, degree and precision.

    The complete-arc, formal-series and exact-validator functions are the
    supplied unchanged mathematical kernels. The schedule avoids making the
    mathematical input depend on the time spent in rejected radius attempts.
    """
    import gzip
    from flint import ctx
    from coupled_exact import primitive, encode_primitives
    from coupled_arb import geometry, formal_coefficients, suggested_degree
    from full_worker import finite_width
    from validate_records import validate_mode
    schedule = next(r for r in compressed('reference_schedule.json.gz') if r['label'] == row['label'])
    degree = schedule['degree']-settings['replay_extra_degree']
    bits = schedule['bits']-settings['replay_extra_bits']
    radius = schedule['radius']
    require(settings['minimum_degree'] <= degree <= settings['maximum_degree'] and
            bits in settings['precision_bits'] and radius in settings['radii'], 'Invalid prescribed parameters')
    directory.mkdir(parents=True, exist_ok=False)
    polynomials = primitive(d, row['key'])
    ctx.prec = settings['precision_bits'][0]
    name = 'geometry_'+radius.replace('/', '_')+'.jsonl.gz'
    with gzip.open(directory/name, 'wt', encoding='utf8') as stream:
        proof = geometry(d, polynomials, Q(radius), 0, dict(settings, bits=ctx.prec, degree=0), stream)
    minimum = suggested_degree(Q(row['weight']), Q(proof['boundary_modulus_upper']), Q(radius),
                               Q(settings['per_mode_tail']), cap=settings['maximum_degree'])
    require(minimum is not None and degree >= max(minimum, settings['minimum_degree']),
            'Prescribed degree does not meet the exact complete-tail budget')
    proof['scalar_tail_upper'] = str(Q(proof['boundary_modulus_upper'])*Q(radius)**(-2*degree-3)/(1-Q(radius)**-2))
    ctx.prec = bits
    coefficients = formal_coefficients(d, polynomials, degree)
    width = finite_width(coefficients, row['weight'])
    require(width <= Q(settings['per_mode_width']), 'Prescribed precision does not meet the exact width budget')
    record = dict(row, source_binding=binding, status='MODE_ARITHMETIC_PASS', radius_attempts=[],
                  arithmetic_attempts=[dict(radius=radius, degree=degree, bits=bits, weighted_width=str(width))],
                  primitives=encode_primitives(polynomials), degree=degree, coefficients=coefficients,
                  geometry=proof, geometry_file=name,
                  weighted_finite_coefficient_interval_width=str(width), settings=dict(settings, degree=degree, bits=bits))
    atomic_new(directory/'mode.json', record)
    validate_mode(record, directory, d, binding)
    return record


def checkpoint(raw, row, phase, envelope):
    marker = raw/'generation_checkpoints'/(phase+'_'+row['label']+'.json')
    if not marker.exists():
        return None
    saved = load(marker)
    require(saved['fingerprint'] == envelope['fingerprint'] and saved['row'] == row,
            'Generation checkpoint has different source, inputs or environment')
    for name, expected in saved['files'].items():
        require(digest(safe_path(raw, name)) == expected, 'Generation raw dependency changed')
    directory = raw/phase/row['label']
    mode = load(directory/'mode.json')
    require(mode['source_binding'] == envelope['payload']['source_binding'], 'Mode source differs')
    for key, value in row.items():
        require(mode[key] == value, 'Mode identity differs')
    for key, value in envelope['payload']['task']['settings'].items():
        require(mode['settings'][key] == value, 'Mode settings differ')
    expected_files = {str((directory/'mode.json').relative_to(raw)),
                      str(safe_path(directory, mode['geometry_file']).relative_to(raw))}
    if phase == 'replay':
        production = raw/'production'/row['label']
        original = load(production/'mode.json')
        geometry = safe_path(production, original['geometry_file'])
        expected_files.update((str((production/'mode.json').relative_to(raw)), str(geometry.relative_to(raw))))
        require(mode['production_mode_sha256'] == digest(production/'mode.json') and
                mode['production_geometry_sha256'] == digest(geometry), 'Replay production dependencies differ')
        settings = envelope['payload']['task']['settings']
        require(mode['degree'] == original['degree']+settings['replay_extra_degree'] and
                mode['settings']['bits'] == original['settings']['bits']+settings['replay_extra_bits'],
                'Replay increments differ')
    require(set(saved['files']) == expected_files, 'Incomplete generation checkpoint dependencies')
    return mode


def worker(args):
    arm_parent_death(args.expected_parent)
    unoptimized()
    raw = external(args.output)
    envelope = load(raw/'generation_binding.json')
    payload = envelope['payload']
    require(object_hash(payload) == envelope['fingerprint'], 'Corrupt generation binding')
    require(payload['package_sha256'] == code_binding(), 'Worker source/input package differs')
    require(payload['environment'] == environment(native=True), 'Worker environment differs')
    require(payload['source_binding'] == source_binding(), 'Worker mathematical source differs')
    require(payload['task'] == load(SOURCE/'full_task.json'), 'Worker task differs')
    inv = load(raw/'complete_inventory.json')
    require(object_hash(inv) == payload['inventory_sha256'], 'Worker inventory differs')
    row = next(r for r in inv['retained'] if r['label'] == args.label)
    settings, d = payload['task']['settings'], candidate()
    import resource
    resource.setrlimit(resource.RLIMIT_AS, (settings['worker_memory_bytes'],)*2)
    resource.setrlimit(resource.RLIMIT_FSIZE, (settings['maximum_file_bytes'],)*2)
    def alarm(signum, frame):
        raise TimeoutError('Per-mode time limit reached')
    signal.signal(signal.SIGALRM, alarm)
    signal.alarm(settings['mode_seconds'])
    from full_worker import replay, finite_width
    from validate_records import validate_mode
    started = time.monotonic()
    try:
        saved = checkpoint(raw, row, args.phase, envelope)
        directory = raw/args.phase/row['label']
        if saved is not None:
            validate_mode(saved, directory, d, payload['source_binding'])
            require(finite_width(saved['coefficients'], row['weight']) <= Q(settings['per_mode_width']),
                    'Resumed width exceeds its budget')
            require(abs(Q(row['weight']))*Q(saved['geometry']['scalar_tail_upper']) <= Q(settings['per_mode_tail']),
                    'Resumed tail exceeds its budget')
            checkpoint(raw, row, args.phase, envelope)
            return
        if directory.exists():
            directory.rename(directory.with_name('incomplete_'+row['label']+'_'+str(time.time_ns())))
        if args.phase == 'production':
            record = scheduled_production(row, d, settings, payload['source_binding'], directory)
        else:
            require(checkpoint(raw, row, 'production', envelope) is not None, 'Missing completed production')
            record = replay(row, d, settings, payload['source_binding'], raw/'production'/row['label'], directory)
        paths = [directory/'mode.json', safe_path(directory, record['geometry_file'])]
        if args.phase == 'replay':
            original_dir = raw/'production'/row['label']
            original = load(original_dir/'mode.json')
            paths += [original_dir/'mode.json', safe_path(original_dir, original['geometry_file'])]
            require(record['production_mode_sha256'] == digest(paths[-2]) and
                    record['production_geometry_sha256'] == digest(paths[-1]), 'Production changed during replay')
        require(payload['package_sha256'] == code_binding(), 'Package changed during generation')
        atomic_new(raw/'generation_checkpoints'/(args.phase+'_'+row['label']+'.json'),
                   dict(fingerprint=envelope['fingerprint'], row=row,
                        files={str(p.relative_to(raw)):digest(p) for p in paths},
                        seconds=time.monotonic()-started))
    finally:
        signal.alarm(0)


def run_phase(args, rows, phase, envelope, deadline):
    raw = args.output
    pending, active, completed, failures = list(rows), {}, [], []
    interruption = None
    def stop(signum, frame):
        nonlocal interruption
        interruption = signum
    previous = {s:signal.signal(s, stop) for s in (signal.SIGTERM, signal.SIGINT)}
    settings = envelope['payload']['task']['settings']
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONUNBUFFERED='1',
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    stopped = False
    try:
        while pending or active:
            if interruption is not None:
                break
            if time.monotonic()+settings['mode_seconds']+10 >= deadline:
                stopped = True
            while pending and len(active) < args.workers and not stopped:
                row = pending.pop(0)
                command = [sys.executable, '-B', str(ROOT/'run_generation.py'), '--output', str(raw),
                           '--label', row['label'], '--phase', phase, '--expected-parent', str(os.getpid())]
                logpath = raw/'logs'/(phase+'_'+row['label']+'_'+str(time.time_ns())+'.txt')
                logpath.parent.mkdir(exist_ok=True)
                log = logpath.open('x')
                def register(child):
                    active[child] = (row, time.monotonic(), log)
                try:
                    protected_spawn(command, register, env=env, stdout=log, stderr=log)
                except BaseException:
                    log.close()
                    raise
            for process, (row, began, log) in list(active.items()):
                code = process.poll()
                if code is None and time.monotonic()-began > settings['mode_seconds']+10:
                    stop_group(process)
                    code = process.returncode
                if code is None:
                    continue
                del active[process]
                log.close()
                if code == 0 and checkpoint(raw, row, phase, envelope) is not None:
                    completed.append(row['label'])
                else:
                    failures.append(dict(label=row['label'], returncode=code))
                if (len(completed)+len(failures)) % 128 == 0:
                    print(dict(phase=phase, completed=len(completed), failed=len(failures), required=len(rows)), flush=True)
                    if sum(p.stat().st_size for p in raw.rglob('*') if p.is_file()) > envelope['payload']['task']['maximum_output_bytes']:
                        stopped = True
            if stopped and not active:
                break
            if active:
                time.sleep(.05)
    finally:
        for process, (_, _, log) in active.items():
            stop_group(process)
            log.close()
        for s, handler in previous.items():
            signal.signal(s, handler)
        atomic_new(raw/'runs'/(phase+'_'+str(time.time_ns())+'.json'),
                   dict(fingerprint=envelope['fingerprint'], selected_labels=[r['label'] for r in rows],
                        completed=completed, failures=failures, interrupted=interruption,
                        stopped_for_limit=stopped, full_candidate_certified=False))
    require(interruption is None and not failures and len(completed) == len(rows),
            'Incomplete generation; retain the directory and resume')


def run(args):
    unoptimized()
    require(sys.platform.startswith('linux'), 'Full numerical execution requires Linux process limits')
    require(1 <= args.workers <= 8, 'Use one to eight workers')
    require(args.stage is not None, 'Choose conditioning or complete')
    task = load(SOURCE/'full_task.json')
    d, package = candidate(), code_binding()
    inv = inventory(d, Q(task['threshold']))
    require(inv == inventory(d, Q(task['threshold']), independently=True) == compressed('inventory.json.gz'),
            'Complete independent inventories or supplied inventory differ')
    payload = dict(package_sha256=package, source_binding=source_binding(), task=task,
                   candidate_sha256=object_hash(d), inventory_sha256=object_hash(inv), environment=environment(native=True))
    envelope = dict(payload=payload, fingerprint=object_hash(payload))
    raw = args.output = external(args.output)
    raw.mkdir(parents=True, exist_ok=True)
    with (raw/'generation.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        existing = raw/'generation_binding.json'
        initial = {'complete_inventory.json':inv,
                   'execution_inputs.json':dict(source_binding=payload['source_binding'], task=task,
                       environment=payload['environment'], generation_fingerprint=envelope['fingerprint']),
                   'generation_binding.json':envelope}
        names = {p.name for p in raw.iterdir()}
        require(args.resume or names == {'generation.lock'}, 'Existing output requires --resume')
        if not existing.exists():
            require(all(n == 'generation.lock' or n in initial or
                        any(n.startswith(k+'.tmp.') for k in initial) for n in names),
                    'Cannot resume numerical records without their source/input binding')
        for name, value in initial.items():
            path = raw/name
            if path.exists():
                require(load(path) == value, 'Resume source, inputs, inventory or environment differ: '+name)
            else:
                atomic_new(path, value)
        require(sum(p.stat().st_size for p in raw.rglob('*') if p.is_file()) <= task['maximum_output_bytes'],
                'Existing generation output exceeds the scheduling disk limit')
        rows = conditioning_rows(inv['retained']) if args.stage == 'conditioning' else inv['retained']
        if args.indices:
            selected = list(map(int, args.indices.split(',')))
            require(len(selected) == len(set(selected)) and all(0 <= i < len(inv['retained']) for i in selected),
                    'Unknown or duplicate retained index')
            rows = [inv['retained'][i] for i in sorted(selected)]
        gate = raw/'conditioning.json'
        if args.stage == 'complete' and not args.indices:
            saved = load(gate)
            expected = conditioning_rows(inv['retained'])
            require(saved['fingerprint'] == envelope['fingerprint'] and saved['labels'] == [r['label'] for r in expected],
                    'This source and input require completed numerical conditioning')
            for row in expected:
                for phase in ('production', 'replay'):
                    require(checkpoint(raw, row, phase, envelope) is not None, 'Conditioning evidence is missing')
        deadline = time.monotonic()+(task['conditioning_wall_seconds'] if args.stage == 'conditioning' else task['wall_seconds'])
        for phase in ('production', 'replay'):
            run_phase(args, rows, phase, envelope, deadline)
        require(package == code_binding(), 'Source or inputs changed during generation')
        if args.stage == 'conditioning' and not args.indices:
            value = dict(fingerprint=envelope['fingerprint'], labels=[r['label'] for r in rows])
            if gate.exists():
                require(load(gate) == value, 'Conditioning gate differs')
            else:
                atomic_new(gate, value)
        print(dict(status='COMPUTATION_COMPLETE', selected_modes=len(rows),
                   full_selection=args.stage == 'complete' and not bool(args.indices),
                   full_candidate_certified=False, next_step='Complete exact validation is required.'), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--stage', choices=('conditioning', 'complete'))
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--indices', default='', help='Partial retained indices for diagnosis only')
    parser.add_argument('--label', help=argparse.SUPPRESS)
    parser.add_argument('--phase', choices=('production', 'replay'), help=argparse.SUPPRESS)
    parser.add_argument('--expected-parent', type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.label:
        worker(args)
    else:
        run(args)
