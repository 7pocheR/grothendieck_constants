"""Compute every coefficient and complete circle enclosure with Arb.

Uses one independent subprocess per mode, up to four concurrently. Complete
checkpoints are reused only when all source, input, environment and record
hashes agree. Interrupted attempts remain incomplete and are never accepted.
"""
import argparse
from contextlib import contextmanager
import gzip
import importlib.metadata
import io
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time
import traceback
sys.path.insert(0,str(Path(__file__).resolve().parent/'sources'))
from exact_core import digest,load,require,specialize
from package_core import (HERE,SETTINGS,check_inputs,fingerprint,safe_path,source_hashes,unoptimized)


def write_new(path,value):
    with path.open('x') as stream:json.dump(value,stream,indent=2);stream.write('\n')

@contextmanager
def deterministic_gzip(path):
    with path.open('xb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=3) as zipped:
            with io.TextIOWrapper(zipped,encoding='utf8',newline='\n') as text:yield text


def validate_checkpoint(output,row,binding):
    folder=output/f"mode_{row['index']:04d}"
    marker=folder/'complete.json'
    if not marker.exists():return False
    completion=load(marker)
    require(completion['fingerprint']==binding,'Checkpoint fingerprint differs')
    attempt=safe_path(folder,completion['attempt'])
    require(digest(attempt/'mode.json')==completion['mode_sha256'] and digest(attempt/'geometry.jsonl.gz')==completion['geometry_sha256'],'Checkpoint bytes differ')
    record=load(attempt/'mode.json')
    require(record['fingerprint']==binding and all(record[k]==v for k,v in row.items()) and len(record['coefficients'])==row['degree']+1,'Checkpoint mathematical identity differs')
    return True


def worker(args):
    unoptimized()
    envelope=load(args.output/'inputs.json');payload=envelope['payload'];binding=fingerprint(payload)
    require(binding==envelope['fingerprint'] and payload['source_hashes']==source_hashes(),'Worker input/source binding differs')
    require(payload['settings']==SETTINGS,'Worker settings differ')
    row=next(r for r in payload['plan']['retained'] if r['index']==args.worker_index)
    if validate_checkpoint(args.output,row,binding):return
    folder=args.output/f"mode_{row['index']:04d}";folder.mkdir(exist_ok=True)
    number=1
    while (folder/f'attempt_{number:03d}').exists():number+=1
    attempt=folder/f'attempt_{number:03d}';attempt.mkdir()
    def interrupted(signum,frame):raise InterruptedError('Replay worker interrupted')
    for signum in (signal.SIGTERM,signal.SIGINT,signal.SIGALRM):signal.signal(signum,interrupted)
    started=time.monotonic()
    try:
        # Full replay is supported on Linux; the native mathematical functions
        # can also be used by smoke_test.py on other supported Python platforms.
        import resource
        resource.setrlimit(resource.RLIMIT_AS,(SETTINGS['worker_memory_bytes'],SETTINGS['worker_memory_bytes']))
        signal.alarm(SETTINGS['mode_seconds'])
        from native import formal_coefficients,geometry,ctx
        from fractions import Fraction as Q
        ctx.prec=SETTINGS['bits'];ctx.cap=row['degree']+1
        primitives=specialize(payload['formula'],payload['data'],list(map(Q,row['key'])))
        coefficients=formal_coefficients(payload['data'],primitives,row['degree'])
        with deterministic_gzip(attempt/'geometry.jsonl.gz') as stream:
            cover=geometry(payload['data'],primitives,row,SETTINGS,stream)
        record=dict(**row,fingerprint=binding,coefficients=coefficients,geometry=cover)
        write_new(attempt/'mode.json',record)
        completion=dict(fingerprint=binding,attempt=attempt.name,mode_sha256=digest(attempt/'mode.json'),geometry_sha256=digest(attempt/'geometry.jsonl.gz'))
        temporary=folder/f'complete_{number:03d}.tmp'
        write_new(temporary,completion)
        os.replace(temporary,folder/'complete.json')
        print(json.dumps(dict(index=row['index'],status='COMPUTED',coefficients=len(coefficients),leaves=cover['leaves'],seconds=time.monotonic()-started)),flush=True)
    except BaseException:
        with (attempt/'failure.txt').open('x') as stream:traceback.print_exc(file=stream)
        raise
    finally:signal.alarm(0)


def run(args):
    unoptimized()
    require(sys.platform.startswith('linux'),'Full replay requires Linux process limits and signals; smoke_test.py is portable')
    require(1<=args.workers<=4,'Use one to four workers')
    data,plan,formula,*_=check_inputs()
    environment=dict(python=sys.version,platform=platform.platform(),python_flint=importlib.metadata.version('python-flint'))
    require(environment['python_flint']=='0.8.0','Use python-flint==0.8.0')
    payload=dict(schema=1,data=data,plan=plan,formula=formula,settings=SETTINGS,source_hashes=source_hashes(),environment=environment)
    binding=fingerprint(payload);envelope=dict(payload=payload,fingerprint=binding)
    args.output=args.output.resolve()
    if args.output.exists():
        require(args.resume,'Existing output requires --resume')
        require(load(args.output/'inputs.json')==envelope,'Resume source, environment, or input changed')
    else:
        args.output.mkdir(parents=True)
        write_new(args.output/'inputs.json',envelope)
    # An exclusive lock prevents simultaneous drivers from mutating checkpoints.
    import fcntl
    lock=(args.output/'driver.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    selected=plan['retained']
    if args.indices:
        indices=list(map(int,args.indices.split(',')))
        require(len(indices)==len(set(indices)) and set(indices)<={r['index'] for r in selected},'Unknown or duplicate selected mode')
        selected=[r for r in selected if r['index'] in indices]
    pending=[];completed=[]
    for row in selected:
        if validate_checkpoint(args.output,row,binding):completed.append(row['index'])
        else:pending.append(row)
    active={};failures=[];interruption=None;started=time.monotonic()
    def request_stop(signum,frame):
        nonlocal interruption
        interruption=signum
    previous={s:signal.signal(s,request_stop) for s in (signal.SIGINT,signal.SIGTERM)}
    try:
        while pending or active:
            if interruption is not None:break
            while pending and len(active)<args.workers:
                row=pending.pop(0)
                process=subprocess.Popen([sys.executable,'-B',str(HERE/'run_replay.py'),'--output',str(args.output),'--worker-index',str(row['index'])],start_new_session=True)
                active[process]=(row,time.monotonic())
            for process,(row,began) in list(active.items()):
                code=process.poll()
                if code is None and time.monotonic()-began>SETTINGS['mode_seconds']+30:
                    try:os.killpg(process.pid,signal.SIGTERM)
                    except ProcessLookupError:pass
                    try:code=process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        try:os.killpg(process.pid,signal.SIGKILL)
                        except ProcessLookupError:pass
                        code=process.wait()
                if code is None:continue
                del active[process]
                if code==0 and validate_checkpoint(args.output,row,binding):completed.append(row['index'])
                else:failures.append(dict(index=row['index'],exit_code=code))
            if active:time.sleep(.1)
    finally:
        # Terminate every worker group, then reap it, before the driver exits.
        for process in active:
            if process.poll() is None:
                try:os.killpg(process.pid,signal.SIGTERM)
                except ProcessLookupError:pass
        deadline=time.monotonic()+5
        for process in active:
            try:process.wait(timeout=max(.01,deadline-time.monotonic()))
            except subprocess.TimeoutExpired:
                try:os.killpg(process.pid,signal.SIGKILL)
                except ProcessLookupError:pass
                process.wait()
        for s,handler in previous.items():signal.signal(s,handler)
        summary=dict(fingerprint=binding,selected_modes=len(selected),completed_modes=sorted(completed),failed_modes=failures,interrupted=interruption is not None,full_selection=not bool(args.indices),seconds=time.monotonic()-started,qualification='Computed checkpoints require complete exact validation; this file does not establish the theorem.')
        write_new(args.output/f'execution_{time.time_ns()}.json',summary)
        lock.close()
    require(interruption is None,'Replay interrupted; completed checkpoints and incomplete attempts preserved')
    require(not failures and len(completed)==len(selected),'Incomplete replay; inspect failed attempts and resume')
    print(json.dumps(dict(status='COMPUTATION_COMPLETE',full_selection=not bool(args.indices),selected_modes=len(selected),next_step='Run verify_certificate.py --replay on this output.')),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--indices',default='',help='Optional partial computation; cannot prove the theorem')
    parser.add_argument('--worker-index',type=int,help=argparse.SUPPRESS)
    args=parser.parse_args()
    if args.worker_index is None:run(args)
    else:worker(args)
