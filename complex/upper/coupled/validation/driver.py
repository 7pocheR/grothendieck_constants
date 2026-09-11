"""Global time bound, Linux lifecycle conditioning, and exact-validator cleanup."""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import uuid

from linux_lifecycle import protected_spawn, stop_group

HERE = Path(__file__).resolve().parent
CHILD = None


def register(child):
    global CHILD
    CHILD = child


def stop():
    global CHILD
    stop_group(CHILD)
    CHILD = None


def terminate(signum, frame):
    # A second termination request cannot interrupt process-group cleanup.
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    stop()
    raise SystemExit(128+signum)


def run(command, env, deadline, maximum=None):
    left = deadline-time.monotonic()
    if maximum is not None:
        left = min(left, maximum)
    if left <= 0:
        return 124
    wait_deadline = time.monotonic()+left
    protected_spawn(command, register, env=env, stdin=subprocess.DEVNULL)
    try:
        return CHILD.wait(timeout=max(0.001, wait_deadline-time.monotonic()))
    except subprocess.TimeoutExpired:
        return 124
    finally:
        stop()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--frozen-source', required=True, type=Path)
    parser.add_argument('--evidence', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--stage', required=True, choices=('conditioning', 'complete'))
    parser.add_argument('--cpus', type=int, choices=(1, 2, 3, 4), default=4)
    parser.add_argument('--reference', action='store_true')
    args = parser.parse_args()
    if not (sys.platform.startswith('linux') and __debug__):
        raise ValueError('Linux with assertions enabled is required')
    started = time.monotonic()
    outer = 295 if args.stage == 'conditioning' else 7195
    deadline = started+outer
    output, evidence = args.output.resolve(), args.evidence.resolve()
    if output.is_relative_to(evidence) or evidence.is_relative_to(output):
        raise ValueError('Output and raw evidence must be disjoint')
    existing_bytes = sum(p.stat().st_size for p in output.rglob('*') if p.is_file())
    if existing_bytes+(1 << 20) > (1 << 30):
        raise ValueError('Insufficient room for bounded controls under the 1 GiB output cap')
    env = dict(os.environ, PYTHONUNBUFFERED='1', PYTHONDONTWRITEBYTECODE='1', PYTHONOPTIMIZE='0',
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    signal.signal(signal.SIGTERM, terminate)
    signal.signal(signal.SIGINT, terminate)
    controls = output/'linux_controls'/uuid.uuid4().hex
    controls_command = [sys.executable, '-B', '-u', str(HERE/'check_linux_controls.py'),
                        '--output', str(controls), '--expected-parent', str(os.getpid())]
    code = run(controls_command, env, deadline-5, maximum=30)
    if code != 0:
        print(json.dumps(dict(status='LINUX_CONTROLS_FAILED', returncode=code,
                              elapsed_seconds=time.monotonic()-started)), flush=True)
        return code
    control_report = json.loads((controls/'report.json').read_text())
    if control_report['status'] != 'LINUX_PROCESS_CONTROLS_PASS':
        raise ValueError('Linux process controls did not pass')
    for name in ('linux_lifecycle.py', 'check_linux_controls.py',
                 'driver.py', 'manifest.json'):
        if control_report['source_sha256'][name] != sha256((HERE/name).read_bytes()).hexdigest():
            raise ValueError('Linux control source binding differs: '+name)
    maximum = 280 if args.stage == 'conditioning' else 7180
    seconds = min(maximum, int(deadline-time.monotonic())-10)
    if seconds < 1:
        return 124
    command = [sys.executable, '-B', '-u', str(HERE/'checkpoint_validator.py'),
               '--frozen-source', str(args.frozen_source.resolve()), '--evidence', str(evidence),
               '--output', str(output), '--stage', args.stage, '--cpus', str(args.cpus),
               '--seconds', str(seconds), '--expected-parent', str(os.getpid())]
    if args.reference:
        command.append('--reference')
    code = run(command, env, deadline-3)
    print(json.dumps(dict(stage=args.stage, returncode=code,
                          linux_controls=str(controls/'report.json'),
                          elapsed_seconds=time.monotonic()-started)), flush=True)
    return code


if __name__ == '__main__':
    raise SystemExit(main())
