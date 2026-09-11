"""Bounded, separately checkpointed exact validation of frozen coupled evidence."""
from linux_lifecycle import (arm_parent_death, guard_spawn_entry, guard_worker,
                             WORKER_PARENT_ENV)
if __name__ == '__mp_main__':
    guard_spawn_entry()

import argparse
import fcntl
from fractions import Fraction as Q
from hashlib import sha256
import importlib.metadata
import json
import multiprocessing as mp
import os
from pathlib import Path
import platform
import signal
import sys
import time
import uuid

HERE = Path(__file__).resolve().parent
CAP = 1 << 30
MEMORY = 3_000_000_000
MODE_SECONDS = 150
CONDITION_INDICES = (0, 1, 2, 1243, 2487, 3730, 4974, 6217, 7461, 8704, 9946, 9947)
CONTRACT = 'unchanged_validate_mode_with_full_replay_binding_and_reconstructed_dyadic_contribution_v1'
TRUST = ('The contour and primitive implications rely on successful execution of the specified '
         'exact validate_mode source on all hashed evidence. A checkpoint is an execution record, '
         'not a proof reconstructed from a Boolean. Resume rechecks every dependency hash and '
         'reconstructs its exact arithmetic summary and dyadic contribution from the raw record.')
G = {}
sys.path.insert(0, str(HERE.parent))
from package_support import verify_package, compressed, environment, external, atomic_new, object_hash as package_object_hash


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':')).encode()


def object_hash(value):
    return sha256(canonical(value)).hexdigest()


def file_hash(path):
    h = sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(path):
    def unique(items):
        value = {}
        for key, item in items:
            require(key not in value, 'Duplicate JSON key: '+key)
            value[key] = item
        return value
    return json.loads(path.read_text(), object_pairs_hook=unique)


def install_source(source):
    sys.path.insert(0, str(source))
    import coupled_exact
    require(coupled_exact.HERE.resolve() == source, 'Imported a different frozen source directory')
    from full_common import source_binding, inventory, DyadicAggregate
    from validate_records import validate_mode
    from winding_exact_core import binary_ball
    return coupled_exact, source_binding, inventory, DyadicAggregate, validate_mode, binary_ball


def reserve(size):
    with G['budget_lock']:
        require(G['used_bytes'].value+size <= CAP, 'Separate output reached its 1 GiB logical-byte cap')
        G['used_bytes'].value += size


def write_new(path, value):
    require(not path.exists(), 'Refusing to overwrite evidence: '+str(path))
    data = json.dumps(value, indent=2)+'\n'
    encoded = data.encode()
    reserve(len(encoded))
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name+'.tmp.'+str(os.getpid())+'.'+uuid.uuid4().hex[:8])
    with temporary.open('xb') as stream:
        stream.write(encoded)
    temporary.replace(path)


def append_progress(path, value):
    encoded = (json.dumps(value, separators=(',', ':'))+'\n').encode()
    reserve(len(encoded))
    with path.open('ab') as stream:
        stream.write(encoded)


def initializer(source, evidence, output, task, d, binding_sha, frozen_binding, own_hash,
                used_bytes, budget_lock, expected_parent, reference):
    guard_worker(expected_parent)
    import resource
    resource.setrlimit(resource.RLIMIT_AS, (MEMORY, MEMORY))
    resource.setrlimit(resource.RLIMIT_FSIZE, (8 << 20, 8 << 20))
    modules = install_source(Path(source))
    G.update(evidence=Path(evidence), output=Path(output), task=task, d=d,
             binding_sha=binding_sha, frozen_binding=frozen_binding, own_hash=own_hash,
             aggregate_class=modules[3], validate_mode=modules[4], binary_ball=modules[5],
             used_bytes=used_bytes, budget_lock=budget_lock,
             reference_hashes=compressed('reference_raw_sha256.json.gz') if reference else None,
             schedule={r['label']:r for r in compressed('reference_schedule.json.gz')})


def geometry_path(directory, record):
    path = (directory/record['geometry_file']).resolve()
    require(path.is_relative_to(directory.resolve()) and path.is_file(), 'Invalid contour path')
    return path


def prepare_mode(row):
    evidence, task = G['evidence'], G['task']
    replay_dir, production_dir = evidence/'replay'/row['label'], evidence/'production'/row['label']
    replay_file, production_file = replay_dir/'mode.json', production_dir/'mode.json'
    record, original = load(replay_file), load(production_file)
    paths = (replay_file, geometry_path(replay_dir, record),
             production_file, geometry_path(production_dir, original))
    hashes = {str(path.relative_to(evidence)): file_hash(path) for path in paths}
    if G['reference_hashes'] is not None:
        require(all(G['reference_hashes'].get(name) == value for name, value in hashes.items()),
                'Supplied original raw file differs from its reference hash')
    # Bind the production metadata as well as the fresh replay metadata.
    for item in (record, original):
        require(item['status'] == 'MODE_ARITHMETIC_PASS', 'Incomplete mode arithmetic')
        require(item['source_binding'] == G['frozen_binding'], 'Mode source binding differs')
        for name in ('label', 'key', 'weight', 'omission_cost'):
            require(item[name] == row[name], 'Mode identity differs: '+name)
        for name, value in task['settings'].items():
            require(item['settings'][name] == value, 'Mode settings differ: '+name)
        require(item['degree'] == item['settings']['degree'] and
                len(item['coefficients']) == item['degree']+1, 'Coefficient degree differs')
    require(record['production_mode_sha256'] == hashes[str(production_file.relative_to(evidence))] and
            record['production_geometry_sha256'] == hashes[str(paths[3].relative_to(evidence))],
            'Replay binds different production evidence')
    require(record['degree'] == original['degree']+task['settings']['replay_extra_degree'] and
            record['settings']['bits'] == original['settings']['bits']+task['settings']['replay_extra_bits'],
            'Replay precision or degree increment differs')
    schedule = G['schedule'][row['label']]
    require(record['degree'] == schedule['degree'] and record['settings']['bits'] == schedule['bits'] and
            record['geometry']['radius'] == original['geometry']['radius'] == schedule['radius'],
            'Prescribed complete reproduction schedule differs')
    intervals = [G['binary_ball'](value) for value in record['coefficients']]
    width = abs(Q(row['weight']))*sum(hi-lo for lo, hi in intervals)
    require(str(width) == record['weighted_finite_coefficient_interval_width'], 'Exact coefficient width differs')
    require(width <= Q(task['settings']['per_mode_width']), 'Finite coefficient width exceeds budget')
    proof = record['geometry']
    radius = Q(proof['radius'])
    require(radius > 1 and proof['radius'] in task['settings']['radii'], 'Invalid tail radius')
    tail = Q(proof['boundary_modulus_upper'])*radius**(-2*record['degree']-3)/(1-radius**-2)
    require(str(tail) == proof['scalar_tail_upper'] and tail >= 0, 'Exact Cauchy tail differs')
    require(abs(Q(row['weight']))*tail <= Q(task['settings']['per_mode_tail']), 'Tail exceeds budget')
    aggregate = G['aggregate_class'](task['aggregate_bits'])
    aggregate.add(record)
    contribution = dict(lower=aggregate.lower, upper=aggregate.upper, tail=aggregate.tail,
                        count=aggregate.count, bits=aggregate.bits)
    summary = dict(weight=row['weight'], degree=record['degree'], bits=record['settings']['bits'],
                   coefficient_count=len(intervals), exact_weighted_width=str(width),
                   exact_unweighted_tail=str(tail), exact_weighted_tail=str(abs(Q(row['weight']))*tail),
                   geometry=proof, encoded_primitives_sha256=object_hash(record['primitives']),
                   dyadic_contribution_sha256=object_hash(contribution),
                   dyadic_coefficient_count=len(aggregate.lower), dyadic_tail=aggregate.tail)
    return record, replay_dir, hashes, summary, contribution, paths


def validate_one(row):
    started = time.monotonic()
    def alarm(signum, frame):
        raise TimeoutError('Bounded per-mode exact validation exceeded its time limit')
    signal.signal(signal.SIGALRM, alarm)
    signal.alarm(MODE_SECONDS)
    path = G['output']/'checkpoints'/(row['label']+'.json')
    try:
        begin = time.monotonic()
        record, directory, hashes, summary, contribution, paths = prepare_mode(row)
        preparation_seconds = time.monotonic()-begin
        descriptor = dict(contract=CONTRACT, binding_sha256=G['binding_sha'], row=row,
                          dependency_sha256=hashes, exact_summary=summary,
                          trusted_execution_dependency=TRUST)
        validation_seconds = 0.0
        resumed = path.exists()
        if resumed:
            saved = load(path)
            require(saved['validated_descriptor'] == descriptor, 'Checkpoint dependency or reconstructed arithmetic differs')
            require(saved['validator_source_sha256'] == G['own_hash'], 'Checkpoint validator source differs')
        else:
            begin = time.monotonic()
            intervals, tail = G['validate_mode'](record, directory, G['d'], G['frozen_binding'])
            validation_seconds = time.monotonic()-begin
            require(len(intervals) == summary['coefficient_count'] and
                    str(tail) == summary['exact_unweighted_tail'], 'Frozen validator result differs')
            require(abs(Q(row['weight']))*sum(hi-lo for lo, hi in intervals) <=
                    Q(G['task']['settings']['per_mode_width']), 'Frozen width budget exceeded')
            require(abs(Q(row['weight']))*tail <= Q(G['task']['settings']['per_mode_tail']), 'Frozen tail budget exceeded')
        for item in paths:
            require(file_hash(item) == hashes[str(item.relative_to(G['evidence']))], 'Evidence changed during validation')
        if not resumed:
            write_new(path, dict(validated_descriptor=descriptor, validator_source_sha256=G['own_hash'],
                      timings=dict(preparation_seconds=preparation_seconds,
                                   exact_validation_seconds=validation_seconds,
                                   total_seconds=time.monotonic()-started)))
        return dict(label=row['label'], status='EXACT_MODE_VALIDATED', resumed=resumed,
                    checkpoint_sha256=file_hash(path), contribution=contribution,
                    degree=record['degree'], leaves=record['geometry']['leaves'],
                    preparation_seconds=preparation_seconds, exact_validation_seconds=validation_seconds,
                    seconds=time.monotonic()-started)
    except Exception as error:
        return dict(label=row['label'], status='UNRESOLVED', error=type(error).__name__+': '+str(error),
                    seconds=time.monotonic()-started)
    finally:
        signal.alarm(0)


def combine(aggregate, contribution):
    require(contribution['bits'] == aggregate.bits and contribution['count'] == 1, 'Wrong contribution normalization')
    require(len(contribution['lower']) == len(contribution['upper']), 'Malformed contribution')
    while len(aggregate.lower) < len(contribution['lower']):
        aggregate.lower.append(0)
        aggregate.upper.append(0)
    for i, (a, b) in enumerate(zip(contribution['lower'], contribution['upper'])):
        require(type(a) is int and type(b) is int and a <= b, 'Invalid integer interval')
        aggregate.lower[i] += a
        aggregate.upper[i] += b
    require(type(contribution['tail']) is int and contribution['tail'] >= 0, 'Invalid tail contribution')
    aggregate.tail += contribution['tail']
    aggregate.count += 1


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--frozen-source', required=True, type=Path)
    ap.add_argument('--evidence', required=True, type=Path)
    ap.add_argument('--output', required=True, type=Path)
    ap.add_argument('--stage', required=True, choices=('conditioning', 'complete'))
    ap.add_argument('--cpus', type=int, choices=(1, 2, 3, 4), default=1)
    ap.add_argument('--seconds', type=int, required=True)
    ap.add_argument('--expected-parent', type=int, required=True)
    ap.add_argument('--reference', action='store_true')
    args = ap.parse_args()
    arm_parent_death(args.expected_parent)
    require(sys.platform.startswith('linux') and __debug__,
            'Linux with assertions enabled is required')
    package_manifest = verify_package()
    limit = 280 if args.stage == 'conditioning' else 7180
    require(1 <= args.seconds <= limit, 'Internal stage time exceeds its bounded proposal')
    source, evidence, output = (p.resolve() for p in (args.frozen_source, args.evidence, args.output))
    external(output)
    require(source == (HERE.parent/'sources').resolve(), 'Use the supplied mathematical sources')
    require(not output.is_relative_to(evidence) and not evidence.is_relative_to(output), 'Output and raw evidence must be disjoint')
    started, deadline = time.monotonic(), time.monotonic()+args.seconds
    module, source_binding, inventory, aggregate_class, _, _ = install_source(source)
    frozen_binding, task, d = source_binding(), load(source/'full_task.json'), module.candidate()
    require(frozen_binding == package_manifest['sources/source_manifest.json'] and
            file_hash(source/'full_task.json') == package_manifest['sources/full_task.json'],
            'Source or task differs from the package manifest')
    record_binding = frozen_binding
    if args.reference:
        reference = load(HERE.parent/'data/reference_binding.json')
        record_binding = reference['source_manifest_sha256']
        task = load(HERE.parent/'data/reference_task.json')
        require(file_hash(HERE.parent/'data/reference_task.json') == reference['task_sha256'], 'Reference task differs')
        for name, expected in reference['evidence_dependency_sha256'].items():
            require(file_hash(evidence/name) == expected, 'Original input metadata differs: '+name)
    own_manifest = load(HERE/'manifest.json')
    for name, expected in own_manifest['files'].items():
        require(file_hash(HERE/name) == expected, 'New validator source manifest mismatch: '+name)
    own_manifest_hash, own_hash = file_hash(HERE/'manifest.json'), file_hash(Path(__file__))
    inv = inventory(d, Q(task['threshold']), independently=True)
    require(load(evidence/'complete_inventory.json') == inv, 'Independent complete inventory differs')
    execution = load(evidence/'execution_inputs.json')
    require(execution['source_binding'] == record_binding and execution['task'] == task, 'Frozen execution dependency differs')
    require(inv == compressed('inventory.json.gz'), 'Supplied complete inventory differs')
    dependency_names = ['complete_inventory.json', 'execution_inputs.json']
    if not args.reference:
        generation = load(evidence/'generation_binding.json')
        require(object_hash(generation['payload']) == generation['fingerprint'] == execution['generation_fingerprint'],
                'Generation execution binding differs')
        require(generation['payload']['source_binding'] == frozen_binding and
                generation['payload']['package_sha256'] == object_hash(package_manifest) and
                generation['payload']['inventory_sha256'] == object_hash(inv) and
                generation['payload']['candidate_sha256'] == object_hash(d) and
                generation['payload']['task'] == task and generation['payload']['environment'] == execution['environment'],
                'Generation source, candidate, settings, inventory or environment differs')
        dependency_names.append('generation_binding.json')
    labels = {row['label'] for row in inv['retained']}
    require({p.name for p in (evidence/'replay').iterdir() if p.is_dir() and p.name.startswith('mode_')} == labels,
            'Unexpected or missing replay mode directories')
    condition_rows = [inv['retained'][i] for i in CONDITION_INDICES]
    binding = dict(contract=CONTRACT, source_manifest_sha256=frozen_binding,
                   evidence_source_binding=record_binding, reference_records=args.reference,
                   effective_task_sha256=object_hash(task), validator_environment=environment(),
                   package_sha256=object_hash(package_manifest),
                   validator_manifest_sha256=own_manifest_hash, validator_source_sha256=own_hash,
                   task_sha256=file_hash(source/'full_task.json'), candidate_sha256=object_hash(d),
                   evidence_dependency_sha256={name:file_hash(evidence/name) for name in
                                              dependency_names},
                   complete_inventory_sha256=object_hash(inv), condition_indices=list(CONDITION_INDICES),
                   aggregate_bits=task['aggregate_bits'], maximum_output_bytes=CAP,
                   per_worker_memory_bytes=MEMORY, per_mode_seconds=MODE_SECONDS,
                   trusted_execution_dependency=TRUST)
    output.mkdir(parents=True, exist_ok=True)
    with (output/'process.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        context = mp.get_context('spawn')
        used = context.Value('Q', sum(p.stat().st_size for p in output.rglob('*') if p.is_file()))
        require(used.value <= CAP, 'Existing separate output already exceeds its byte cap')
        budget_lock = context.Lock()
        G.update(used_bytes=used, budget_lock=budget_lock)
        binding_path = output/'binding.json'
        if binding_path.exists():
            require(load(binding_path) == binding, 'Cannot resume under different dependencies')
        else:
            write_new(binding_path, binding)
        binding_sha = file_hash(binding_path)
        if args.stage == 'complete':
            gate = load(output/'conditioning_acceptance.json')
            require(gate['binding_sha256'] == binding_sha and
                    set(gate['checkpoint_sha256']) == {r['label'] for r in condition_rows}, 'Conditioning gate differs')
            for label, expected in gate['checkpoint_sha256'].items():
                require(file_hash(output/'checkpoints'/(label+'.json')) == expected, 'Conditioning checkpoint changed')
        rows = condition_rows if args.stage == 'conditioning' else inv['retained']
        run = output/'runs'/(time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())+'_'+uuid.uuid4().hex[:8])
        run.mkdir(parents=True, exist_ok=False)
        write_new(run/'invocation.json', dict(binding_sha256=binding_sha, stage=args.stage,
                  cpus=args.cpus, seconds=args.seconds, python=sys.version, platform=platform.platform(),
                  frozen_outward_environment=execution['environment'], epoch=time.time()))
        aggregate = aggregate_class(task['aggregate_bits'])
        completed, failures, pending = {}, [], {}
        iterator, exhausted, stop = iter(rows), False, False
        previous_worker_parent = os.environ.get(WORKER_PARENT_ENV)
        os.environ[WORKER_PARENT_ENV] = str(os.getpid())
        pool = context.Pool(args.cpus, initializer=initializer,
                 initargs=(str(source),str(evidence),str(output),task,d,binding_sha,record_binding,own_hash,
                           used,budget_lock,os.getpid(),args.reference))
        try:
            while pending or not exhausted:
                if time.monotonic()+5 >= deadline:
                    stop = True
                    break
                while not exhausted and len(pending) < args.cpus:
                    try:
                        row = next(iterator)
                    except StopIteration:
                        exhausted = True
                        break
                    pending[row['label']] = pool.apply_async(validate_one, (row,))
                ready = [label for label, value in pending.items() if value.ready()]
                if not ready:
                    time.sleep(0.05)
                    continue
                for label in ready:
                    result = pending.pop(label).get()
                    contribution = result.pop('contribution', None)
                    if result['status'] == 'EXACT_MODE_VALIDATED':
                        require(label not in completed and contribution is not None, 'Duplicate or missing exact contribution')
                        combine(aggregate, contribution)
                        completed[label] = result['checkpoint_sha256']
                    else:
                        failures.append(result)
                    append_progress(run/'progress.jsonl', result)
                    if (len(completed)+len(failures)) % 16 == 0:
                        print(json.dumps(dict(stage=args.stage,validated=len(completed),unresolved=len(failures),
                                              required=len(rows),elapsed=time.monotonic()-started)), flush=True)
        finally:
            pool.terminate()
            pool.join()
            if previous_worker_parent is None:
                os.environ.pop(WORKER_PARENT_ENV, None)
            else:
                os.environ[WORKER_PARENT_ENV] = previous_worker_parent
        require(source_binding() == frozen_binding, 'Frozen source changed during validation')
        require(file_hash(HERE/'manifest.json') == own_manifest_hash, 'Validator manifest changed during validation')
        for name, expected in own_manifest['files'].items():
            require(file_hash(HERE/name) == expected, 'Validator source changed during validation: '+name)
        for name, expected in binding['evidence_dependency_sha256'].items():
            require(file_hash(evidence/name) == expected, 'Global evidence dependency changed')
        complete = len(completed) == len(rows) and not failures and not stop
        report = dict(status='PARTIAL_EXACT_VALIDATION', binding_sha256=binding_sha, stage=args.stage,
                      validated=len(completed), required=len(rows), unresolved=failures,
                      unvisited_or_interrupted=len(rows)-len(completed)-len(failures),
                      checkpoint_sha256=completed, elapsed_seconds=time.monotonic()-started,
                      logical_output_bytes_before_report=used.value, full_candidate_certified=False,
                      trusted_execution_dependency=TRUST)
        if complete and args.stage == 'conditioning':
            report['status'] = 'CONDITIONING_EXACT_VALIDATION_PASS'
            gate_path = output/'conditioning_acceptance.json'
            gate = dict(binding_sha256=binding_sha, checkpoint_sha256=completed)
            if gate_path.exists():
                require(load(gate_path) == gate, 'Conditioning acceptance differs')
            else:
                write_new(gate_path, gate)
        elif complete:
            require(set(completed) == labels, 'Complete inventory coverage differs')
            require({p.stem for p in (output/'checkpoints').glob('mode_*.json')} == labels,
                    'Unexpected or missing checkpoint labels')
            implication = aggregate.finish(inv, d, Q(task['target']))
            report.update(status='COMPLETE_CRITERION_PASS' if implication['sufficient_criterion_passed']
                          else 'COMPLETE_CRITERION_INSUFFICIENT', implication=implication,
                          full_candidate_certified=implication['sufficient_criterion_passed'],
                          qualification='Same complete sufficient criterion; no exactness or significant-gap claim.')
        write_new(run/'report.json', report)
        if report['status'] == 'COMPLETE_CRITERION_PASS':
            atomic_new(output/'complete.json', report, replace=True)
        print(json.dumps({k:v for k,v in report.items() if k not in ('checkpoint_sha256','unresolved')},indent=2),flush=True)
        return 0 if report['status'] in ('CONDITIONING_EXACT_VALIDATION_PASS','COMPLETE_CRITERION_PASS') else 2


if __name__ == '__main__':
    raise SystemExit(main())
