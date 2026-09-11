"""Linux process-parent protection and signal-safe child registration."""
import ctypes
import os
import signal
import subprocess
import sys
import threading

SHUTDOWN = (signal.SIGTERM, signal.SIGINT)
WORKER_PARENT_ENV = 'COUPLED_EXACT_VALIDATOR_WORKER_PARENT_PID'
PR_SET_PDEATHSIG = 1
PR_GET_PDEATHSIG = 2
_LIBC = None


def linux_prctl():
    global _LIBC
    if not sys.platform.startswith('linux'):
        raise RuntimeError('Linux parent-death protection is required')
    if _LIBC is None:
        _LIBC = ctypes.CDLL(None, use_errno=True)
        _LIBC.prctl.restype = ctypes.c_int
        _LIBC.prctl.argtypes = (ctypes.c_int, ctypes.c_ulong, ctypes.c_ulong,
                               ctypes.c_ulong, ctypes.c_ulong)
    return _LIBC.prctl


def parent_death_signal():
    value = ctypes.c_int()
    if linux_prctl()(PR_GET_PDEATHSIG, ctypes.addressof(value), 0, 0, 0) != 0:
        raise OSError(ctypes.get_errno(), 'PR_GET_PDEATHSIG failed')
    return value.value


def arm_parent_death(expected_parent):
    """No task work is allowed after a missing parent or failed kernel guard."""
    if type(expected_parent) is not int or expected_parent <= 1:
        raise RuntimeError('An explicit live expected parent PID is required')
    prctl = linux_prctl()
    if prctl(PR_SET_PDEATHSIG, int(signal.SIGKILL), 0, 0, 0) != 0:
        raise OSError(ctypes.get_errno(), 'PR_SET_PDEATHSIG failed')
    # PR_SET_PDEATHSIG alone misses a parent that died before this call.
    # Reparenting is permanent here: these processes never change credentials,
    # reparent themselves, or enter a new PID namespace.
    if os.getppid() != expected_parent or parent_death_signal() != signal.SIGKILL:
        os.kill(os.getpid(), signal.SIGKILL)
        os._exit(125)


def unblock_shutdown():
    signal.pthread_sigmask(signal.SIG_UNBLOCK, SHUTDOWN)
    current = signal.pthread_sigmask(signal.SIG_BLOCK, ())
    if any(s in current for s in SHUTDOWN):
        raise RuntimeError('A child inherited a blocked shutdown signal')


def guard_worker(expected_parent):
    arm_parent_death(expected_parent)
    unblock_shutdown()


def guard_spawn_entry():
    """Called before ordinary imports when multiprocessing reloads __main__."""
    guard_worker(int(os.environ[WORKER_PARENT_ENV]))


def protected_spawn(command, register, **kwargs):
    """Only a single-threaded supervisor may use the pre-exec guard.

    register(child) runs with SIGTERM/SIGINT blocked. The callback must retain
    the child before returning, so any pending handler can clean it up when
    the parent restores its previous mask. The child arms its kernel guard
    before unblocking the inherited mask and before exec.
    """
    if threading.active_count() != 1:
        raise RuntimeError('The pre-exec supervisor must be single-threaded')
    if 'preexec_fn' in kwargs or 'start_new_session' in kwargs:
        raise ValueError('Process creation controls cannot be overridden')
    linux_prctl()  # Resolve libc before fork, while the supervisor is single-threaded.
    expected_parent = os.getpid()
    def child_setup():
        arm_parent_death(expected_parent)
        unblock_shutdown()
    previous = signal.pthread_sigmask(signal.SIG_BLOCK, SHUTDOWN)
    try:
        child = subprocess.Popen(command, start_new_session=True,
                                 preexec_fn=child_setup, **kwargs)
        try:
            register(child)
        except BaseException:
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            child.wait()
            raise
    finally:
        signal.pthread_sigmask(signal.SIG_SETMASK, previous)
    return child


def stop_group(child):
    if child is None:
        return
    try:
        os.killpg(child.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        child.wait(timeout=3)
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(child.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    child.wait()
