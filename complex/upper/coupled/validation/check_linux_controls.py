"""Bounded Linux-only synthetic lifecycle tests; no numerical evidence is read."""
from linux_lifecycle import (arm_parent_death, guard_spawn_entry, guard_worker,
                             linux_prctl, parent_death_signal, protected_spawn,
                             stop_group, WORKER_PARENT_ENV)
if __name__ == '__mp_main__':
    guard_spawn_entry()

import argparse
from hashlib import sha256
import json
import multiprocessing as mp
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
CHILD = None
ROOTS = []
KNOWN = {}


def write(path, value):
    encoded = (json.dumps(value, sort_keys=True)+'\n').encode()
    assert len(encoded) <= 4096
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as stream:
        stream.write(encoded)


def identity(pid):
    try:
        text = Path('/proc')/str(pid)/'stat'
        fields = text.read_text().rsplit(')', 1)[1].split()
        return dict(pid=pid, state=fields[0], ppid=int(fields[1]), starttime=int(fields[19]))
    except FileNotFoundError:
        return None


def descendants(pid):
    found = {}
    todo = [pid]
    while todo:
        parent = todo.pop()
        item = identity(parent)
        if item is None or parent in found:
            continue
        found[parent] = item
        try:
            children = (Path('/proc')/str(parent)/'task'/str(parent)/'children').read_text()
        except FileNotFoundError:
            children = ''
        todo.extend(map(int, children.split()))
    return found


def remember(pid):
    KNOWN.update(descendants(pid))


def same_process(item):
    current = identity(item['pid'])
    return current is not None and current['starttime'] == item['starttime']


def reap():
    while True:
        try:
            pid, _ = os.waitpid(-1, os.WNOHANG)
        except ChildProcessError:
            return
        if pid == 0:
            return


def wait_until(predicate, seconds=4):
    end = time.monotonic()+seconds
    while time.monotonic() < end:
        if predicate():
            return
        time.sleep(0.01)
    raise TimeoutError('Synthetic Linux process condition did not occur')


def absent(items):
    # Reaping is limited to the test harness, which is explicitly a subreaper.
    reap()
    return not any(same_process(item) for item in items.values())


def controls_state():
    return dict(pid=os.getpid(), ppid=os.getppid(), pdeathsig=parent_death_signal(),
                blocked_shutdown=sorted(int(s) for s in signal.pthread_sigmask(signal.SIG_BLOCK, ())
                                        if s in (signal.SIGTERM, signal.SIGINT)))


def command(role, output, expected):
    return [sys.executable, '-B', '-u', str(Path(__file__).resolve()), '--role', role,
            '--output', str(output), '--expected-parent', str(expected)]


def register(child):
    global CHILD
    CHILD = child


def child_stop(signum, frame):
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    stop_group(CHILD)
    raise SystemExit(128+signum)


def worker_init(expected):
    guard_worker(expected)


def worker_hold(output):
    write(Path(output)/('worker_'+str(os.getpid())+'.json'), controls_state())
    while True:
        time.sleep(0.05)


def coordinator(output):
    os.environ[WORKER_PARENT_ENV] = str(os.getpid())
    context = mp.get_context('spawn')
    pool = context.Pool(2, initializer=worker_init, initargs=(os.getpid(),))
    for _ in range(2):
        pool.apply_async(worker_hold, (str(output),))
    write(output/'coordinator.json', controls_state())
    wait_until(lambda: len(list(output.glob('worker_*.json'))) == 2)
    write(output/'ready.json', dict(pid=os.getpid()))
    while True:
        time.sleep(0.05)


def simulated_driver(output, window):
    signal.signal(signal.SIGTERM, child_stop)
    signal.signal(signal.SIGINT, child_stop)
    if window:
        def pending_then_register(child):
            global CHILD
            assert CHILD is None
            write(output/'window_child.json', dict(pid=child.pid, identity=identity(child.pid)))
            os.kill(os.getpid(), signal.SIGTERM)
            # A pending signal must not enter the handler before registration.
            assert CHILD is None
            CHILD = child
            write(output/'registered.json', dict(pid=child.pid, pending=int(signal.SIGTERM) in
                                                 [int(s) for s in signal.sigpending()]))
        protected_spawn(command('leaf', output, os.getpid()), pending_then_register,
                        stdin=subprocess.DEVNULL)
        raise AssertionError('The pending SIGTERM did not terminate the driver')
    protected_spawn(command('coordinator', output, os.getpid()), register, stdin=subprocess.DEVNULL)
    write(output/'driver.json', dict(controls_state(), child_pid=CHILD.pid))
    # Reap a terminated coordinator but never clean up its group here: worker
    # termination must follow its kernel guard in the coordinator-death test.
    while True:
        CHILD.poll()
        time.sleep(0.05)


def launch(role, output, expected=None):
    def registered(child):
        ROOTS.append(child)
        remember(child.pid)
    return protected_spawn(command(role, output, os.getpid() if expected is None else expected),
                           registered, stdin=subprocess.DEVNULL)


def cleanup():
    for child in ROOTS:
        remember(child.pid)
    # Kill descendants first, including Python's resource tracker. Normally it
    # has already exited on pipe EOF; this is failure-only cleanup.
    for item in reversed(list(KNOWN.values())):
        if same_process(item):
            try:
                os.kill(item['pid'], signal.SIGKILL)
            except ProcessLookupError:
                pass
    for child in ROOTS:
        try:
            child.wait(timeout=0.5)
        except subprocess.TimeoutExpired:
            pass
    for _ in range(50):
        reap()
        if not any(same_process(item) for item in KNOWN.values()):
            return
        time.sleep(0.01)


def harness(output):
    started = time.monotonic()
    # Adopt and reap intentionally orphaned synthetic descendants. This is a
    # property of the test harness only, never of the production validator.
    if linux_prctl()(36, 1, 0, 0, 0) != 0:  # PR_SET_CHILD_SUBREAPER
        raise RuntimeError('Could not install the synthetic test subreaper')
    def interrupted(signum, frame):
        raise TimeoutError('Linux controls exceeded 25 seconds or were interrupted')
    signal.signal(signal.SIGALRM, interrupted)
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    signal.alarm(25)
    results = {}
    try:
        path = output/'registration_window'
        proc = launch('window', path)
        assert proc.wait(timeout=5) == 128+signal.SIGTERM
        data = json.loads((path/'registered.json').read_text())
        assert data['pending'] is True
        child = json.loads((path/'window_child.json').read_text())['identity']
        assert child is not None
        KNOWN[child['pid']] = child
        wait_until(lambda: absent({child['pid']: child}))
        results['pending_signal_after_creation_before_registration'] = True
        for killed in ('driver', 'coordinator'):
            path = output/('kill_'+killed)
            proc = launch('driver', path)
            wait_until(lambda: (path/'ready.json').exists())
            tree = descendants(proc.pid)
            KNOWN.update(tree)
            driver = json.loads((path/'driver.json').read_text())
            coord = json.loads((path/'coordinator.json').read_text())
            workers = [json.loads(p.read_text()) for p in path.glob('worker_*.json')]
            assert len(workers) == 2 and len(tree) >= 5, tree
            for row in [driver, coord]+workers:
                assert row['pdeathsig'] == signal.SIGKILL and row['blocked_shutdown'] == [], row
            if killed == 'driver':
                os.kill(proc.pid, signal.SIGKILL)
                assert proc.wait(timeout=4) == -signal.SIGKILL
                wait_until(lambda: absent(tree))
            else:
                remaining = {pid: row for pid, row in tree.items() if pid != proc.pid}
                os.kill(coord['pid'], signal.SIGKILL)
                # The driver remains alive and never invokes stop_group here.
                wait_until(lambda: absent(remaining))
                assert identity(proc.pid)['state'] != 'Z'
                os.kill(proc.pid, signal.SIGKILL)
                assert proc.wait(timeout=4) == -signal.SIGKILL
            results[killed+'_sigkill_removes_descendants_including_resource_tracker'] = True
        path = output/'wrong_parent'
        proc = launch('leaf', path, expected=os.getpid()+1000000)
        assert proc.wait(timeout=4) == -signal.SIGKILL
        assert not (path/'leaf.json').exists()
        results['wrong_expected_parent_prevents_task_work'] = True
        results['coordinator_and_workers_have_unblocked_shutdown_signals'] = True
        results['coordinator_and_workers_have_verified_sigkill_parent_death'] = True
        cleanup()
        assert not any(same_process(row) for row in KNOWN.values())
        report = dict(status='LINUX_PROCESS_CONTROLS_PASS', checks=results,
                      physical_modes_evaluated=0, elapsed_seconds=time.monotonic()-started,
                      tracked_synthetic_processes=len(KNOWN),
                      source_sha256={name: sha256((HERE/name).read_bytes()).hexdigest() for name in
                                     ('linux_lifecycle.py', 'check_linux_controls.py',
                                      'driver.py', 'manifest.json')})
        write(output/'report.json', report)
        assert sum(p.stat().st_size for p in output.rglob('*') if p.is_file()) < (1 << 20)
        print(json.dumps(report), flush=True)
    finally:
        signal.alarm(0)
        cleanup()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--role', choices=('harness', 'window', 'driver', 'coordinator', 'leaf'), default='harness')
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--expected-parent', required=True, type=int)
    args = parser.parse_args()
    if not __debug__:
        raise ValueError('Assertions must be enabled for the synthetic controls')
    arm_parent_death(args.expected_parent)
    resource.setrlimit(resource.RLIMIT_AS, (256 << 20, 256 << 20))
    resource.setrlimit(resource.RLIMIT_FSIZE, (65536, 65536))
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.role == 'harness':
        harness(output)
    elif args.role in ('window', 'driver'):
        simulated_driver(output, args.role == 'window')
    elif args.role == 'coordinator':
        coordinator(output)
    else:
        write(output/'leaf.json', controls_state())
        while True:
            time.sleep(0.05)


if __name__ == '__main__':
    main()
