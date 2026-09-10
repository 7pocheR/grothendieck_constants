"""Portable bounded orchestration of the unchanged Gaussian integration kernel."""
import os
import sys
if sys.flags.optimize or os.environ.get('PYTHONOPTIMIZE'):
    raise RuntimeError('Python optimization is unsupported; remove -O/-OO and PYTHONOPTIMIZE.')

import argparse
from concurrent.futures import ThreadPoolExecutor
import importlib.util
import json
from pathlib import Path
import signal
import subprocess
import threading
import time

import verify_certificate as check

ROOT = Path(__file__).resolve().parent
CHILDREN = set()
CHILD_LOCK = threading.RLock()
STOPPING = threading.Event()


def request_child_termination():
    STOPPING.set()
    with CHILD_LOCK:
        children = list(CHILDREN)
    for process in children:
        if process.poll() is None:
            process.terminate()


def terminate_and_reap(process):
    if process.poll() is None:
        process.terminate()
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait()


def install_signal_handlers():
    def interrupted(signum, frame):
        # Do not call wait() from a signal handler: the interrupted main
        # thread may already hold Popen's wait lock. Each invoke() reaps its
        # own child after unwinding; the executor waits for all invokes.
        request_child_termination()
        raise SystemExit(128+signum)
    signal.signal(signal.SIGINT, interrupted)
    signal.signal(signal.SIGTERM, interrupted)


def invoke(arguments, timeout=None):
    # Always create/wait for children outside the main signal-handling thread.
    # The registration lock then closes the creation/registration race without
    # passing a blocked signal mask into the new child process.
    if threading.current_thread() is threading.main_thread():
        with ThreadPoolExecutor(max_workers=1) as pool:
            return pool.submit(invoke, arguments, timeout).result()
    env = dict(os.environ)
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        env[key] = '1'
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    env['PYTHONUNBUFFERED'] = '1'
    process = None
    try:
        with CHILD_LOCK:
            check.require(not STOPPING.is_set(), 'Replay has been interrupted')
            process = subprocess.Popen([sys.executable, '-B', *map(str, arguments)], env=env)
            CHILDREN.add(process)
        deadline = None if timeout is None else time.monotonic()+timeout
        while True:
            if STOPPING.is_set():
                terminate_and_reap(process)
                return process.returncode
            remaining = None if deadline is None else deadline-time.monotonic()
            if remaining is not None and remaining <= 0:
                terminate_and_reap(process)
                return 124
            try:
                return process.wait(timeout=.2 if remaining is None else min(.2, remaining))
            except subprocess.TimeoutExpired:
                continue
    except BaseException:
        if process is not None:
            terminate_and_reap(process)
        raise
    finally:
        if process is not None:
            with CHILD_LOCK:
                CHILDREN.discard(process)


def fresh_domains(manifest, directory):
    # The source hash was checked before import. No code is modified or patched.
    spec = importlib.util.spec_from_file_location('quintic_outward_kernel', ROOT/'check_quintic_01.py')
    kernel = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(kernel)
    domains = {}
    for label, target in check.read(ROOT/'targets_01.json').items():
        _, matrices, compact = kernel.model(target, domains=True)
        domains[label] = {'matrices': matrices, 'compact_domain': compact}
    check.save(directory/'domains.json', {'checker_sha256': manifest['integration_checker_sha256'],
               'inventory_sha256': manifest['inventory_content_sha256'], 'domain_proofs': domains})


def main():
    install_signal_handlers()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['scalar', 'smoke', 'prepare', 'replay', 'finish', 'all'], default='all')
    parser.add_argument('--work', type=Path, default=ROOT/'verification/fresh')
    parser.add_argument('--workers', type=int, choices=range(1, 5), default=4)
    parser.add_argument('--seconds', type=int, default=240)
    parser.add_argument('--max-new', type=int, default=4000)
    parser.add_argument('--rounds', type=int, default=100)
    args = parser.parse_args()
    check.require(1 <= args.seconds <= 240 and 1 <= args.max_new <= 10000 and
                  1 <= args.rounds <= 100, 'Bounded replay limits exceeded')
    manifest = check.manifest_check()
    work = args.work.resolve()
    # Generated files cannot overwrite any distributed source or certificate.
    for relative in manifest['files']:
        check.require(not (ROOT/relative).resolve().is_relative_to(work), 'Work directory contains distributed files')
    certificates = (ROOT/'certificates').resolve()
    check.require(not certificates.is_relative_to(work) and not work.is_relative_to(certificates),
                  'Work directory overlaps distributed certificates')
    work.mkdir(parents=True, exist_ok=True)
    if args.stage == 'scalar':
        return invoke([ROOT/'check_scalar_01.py', '--sources', ROOT/'inputs/scalars', '--output', work/'scalars'])
    if args.stage == 'smoke':
        status = invoke([ROOT/'check_quintic_01.py', 'smoke', '--output', work/'smoke.json'])
        check.require(status == 0, 'Fresh smoke integration failed')
        check.save(work/'smoke_scope.json', {'complete_universal_replay': False,
                   'fresh_vertices': 6, 'fresh_compact_domain_columns': 18,
                   'checker_sha256': manifest['integration_checker_sha256']})
        return 0
    inventory = check.unpack_artifact(manifest, 'inventory.json', work/'inventory.json')
    binding = {'checker_sha256': manifest['integration_checker_sha256'],
               'inventory_sha256': manifest['inventory_content_sha256'], 'parts': 4,
               'precision': 256, 'grid_bits': 52, 'tolerance': manifest['integral_tolerance']}
    binding_path = work/'run_binding.json'
    if binding_path.exists():
        check.require(check.read(binding_path) == binding, 'Working directory belongs to another replay')
    else:
        check.save(binding_path, binding)
    replay = work/'replay'
    replay.mkdir(exist_ok=True)
    check.save(work/'run_status.json', {'complete_universal_replay': False,
               'reason': 'Replay or final validation is in progress; successful completion is written last.'})
    geometry_path = work/'geometry.json'
    # Geometry is checked independently, once per invocation that prepares.
    if args.stage in ('prepare', 'all') or not geometry_path.exists():
        inv = check.load_inventory(manifest, inventory)
        check.save(geometry_path, {'inventory_sha256': inv['content_sha256'],
                                  'geometry': check.geometry(inv)})
        del inv
    if args.stage == 'prepare':
        print(json.dumps({'prepared': True, 'complete_universal_replay': False}))
        return 0

    def complete(part):
        path = replay/f'progress_{part:02d}.json'
        if not path.exists():
            return False
        report = check.read(path)
        expected_binding = {**binding, 'part': part}
        check.require(report['binding'] == expected_binding, 'Progress source binding mismatch')
        return report['complete'] is True

    if args.stage in ('replay', 'all'):
        for round_number in range(args.rounds):
            todo = [p for p in range(4) if not complete(p)]
            if not todo:
                break

            def part_run(part):
                return invoke([ROOT/'check_quintic_01.py', 'replay', '--inventory', inventory,
                               '--output', replay, '--part', part, '--parts', 4,
                               '--seconds', args.seconds, '--max-new', args.max_new], args.seconds+90)

            with ThreadPoolExecutor(max_workers=args.workers) as pool:
                codes = list(pool.map(part_run, todo))
            check.require(all(code == 0 for code in codes), 'A bounded integration worker failed')
            print(json.dumps({'round': round_number+1, 'parts_complete': sum(complete(p) for p in range(4))}), flush=True)
    ready = all(complete(p) for p in range(4))
    if not ready:
        check.save(work/'run_status.json', {'complete_universal_replay': False,
                   'reason': 'Bounded invocation ended before all four parts completed.',
                   'resume': 'Repeat the same command with the same work directory.'})
        print('Incomplete bounded replay; repeat with the same work directory.', flush=True)
        return 2
    # Completion flags are not accepted as proof. Revalidate every actual record.
    fresh_domains(manifest, replay)
    inv = check.load_inventory(manifest, inventory)
    report, result = check.verify(manifest, inv, replay_dir=replay)
    check.save(replay/'verified.json', result)
    check.save(work/'final_validation.json', report)
    check.save(work/'run_status.json', {'complete_universal_replay': True,
               'vertices': 252021, 'leaves': 47399,
               'checker_sha256': manifest['integration_checker_sha256'],
               'inventory_sha256': manifest['inventory_content_sha256'],
               'fresh_domain_columns': 18,
               'interpretation': 'Complete source-bound outward execution and exact final assembly; analytic implications are in PROOF_01.md.'})
    print(json.dumps({'complete_universal_replay': True, 'vertices': 252021, 'leaves': 47399}), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
