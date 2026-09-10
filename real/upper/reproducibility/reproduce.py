"""Checkpointable complete reproduction using the unchanged numerical engines.

This program does not install software, connect to a remote machine, or
submit jobs. Run it only on computational resources already assigned to you.
"""
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ThreadPoolExecutor
import argparse, hashlib, importlib.util, json, os, shutil, subprocess, sys
from flint import arb, ctx

if sys.flags.optimize: raise SystemExit('Do not use python -O.')
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def get(p): return json.loads(p.read_text())
def put(p,v): p.write_text(json.dumps(v,indent=2)+'\n')
def endpoint(v,side): return (F(v['m'])+side*F(v['r']))*F(10)**v['e']

def prepare(work):
    assert not work.is_relative_to(ROOT), 'Choose a separate output directory.'
    manifest=get(ROOT/'immutable_manifest.json')
    copied=[]
    for name in ('head','norm'):
        for folder in ('checks','sources','evidence','run','prior'):
            (work/name/folder).mkdir(parents=True,exist_ok=True)
        for folder in ('checks','sources'):
            for src in sorted((ROOT/name/folder).glob('*')):
                if not src.is_file(): continue
                rel=src.relative_to(ROOT)
                assert sha(src)==manifest[str(rel)]
                dst=work/rel
                if dst.exists(): assert sha(dst)==sha(src), str(dst)
                else: shutil.copyfile(src,dst)
                copied.append(str(rel))
    for filename in ('candidate_rational.json','certify_head.py','assemble_head.py'):
        src=ROOT/'head/run'/filename; dst=work/'head/run'/filename
        assert sha(src)==manifest[str(src.relative_to(ROOT))]
        if dst.exists(): assert sha(dst)==sha(src)
        else: shutil.copyfile(src,dst)
        copied.append(str(dst.relative_to(work)))
    src=ROOT/'head/prior/gauss_nodes_112.json'
    dst=work/'head/prior/gauss_nodes_112.json'
    assert sha(src)==manifest[str(src.relative_to(ROOT))]
    if dst.exists(): assert sha(dst)==sha(src)
    else: shutil.copyfile(src,dst)
    copied.append(str(dst.relative_to(work)))
    (work/'logs').mkdir(exist_ok=True)
    binding={p:sha(work/p) for p in copied}
    old=work/'reproduction_inputs.json'
    if old.exists(): assert get(old)==binding
    else: put(old,binding)
    return binding

def command(work,label,args):
    log=work/'logs'/f'{label}.log'
    # Number repeated logs so that checkpoint resumes preserve old evidence.
    if log.exists():
        k=1
        while log.with_name(f'{label}.{k}.log').exists(): k+=1
        log=log.with_name(f'{label}.{k}.log')
    argv=[sys.executable,'-B','-u',*map(str,args)]
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
    with log.open('w') as out:
        result=subprocess.run(argv,stdout=out,stderr=subprocess.STDOUT,env=env)
    put(log.with_suffix('.json'),{'argv':argv,'exit_code':result.returncode,
        'script_sha256':sha(Path(args[0])),'log_sha256':sha(log)})
    print(label,'exit',result.returncode,'log',log,flush=True)
    if result.returncode: raise RuntimeError(f'{label} failed; inspect {log}')

def finish(work):
    ctx.prec=384
    spec=importlib.util.spec_from_file_location('portable_verification',ROOT/'verify.py')
    verify=importlib.util.module_from_spec(spec);spec.loader.exec_module(verify)
    head=get(work/'head/evidence/independent_head.json')
    assert head['moments_sha256']==sha(work/'head/run/head_moments_501.json')
    assert head['checker_sha256']==sha(work/'head/checks/independent.py')
    assert head['input_sha256']==sha(work/'head/sources/candidate_rational.json')
    assert set(head['coefficients'])=={str(n) for n in range(1,502,2)}
    b1=endpoint(head['b1'],-1)
    mass=sum((max(abs(endpoint(v,-1)),abs(endpoint(v,1)))
              for k,v in head['coefficients'].items() if k!='1'),F(0))
    B=F('0.88182999999999823752');h=F('0.000000001054')
    assert b1>B and mass<h
    circle=get(work/'norm/evidence/independent_circle.json')
    integrals=get(work/'norm/evidence/independent_integrals.json')
    for item in (circle,integrals):
        assert item['checker_sha256']==sha(work/'norm/checks/independent.py')
        assert item['input_sha256']==sha(work/'head/sources/candidate_rational.json')
    caps,q=verify.check_circle_records(circle,integrals)
    C=F('24744.587741');T=F('0.000002349547147305');D=14*501**7
    g=F('0.8818276493988');K=F('1.781296297373')
    assert q<C*C<T*T*D and B-h-T>g>0
    def atan(x):
        a=sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(64)),F(0))
        return a,a+x**129/129
    a,b=atan(F(1,5));c,d=atan(F(1,239))
    pi_hi=16*b-4*c;pi_U=F('3.141592653589793238462643383280')
    assert pi_hi<pi_U<2*g*K
    out={'moment_count':63252,'odd_coefficient_count':251,'panel_count':512,
         'smoothing_integral_count':266,'b1_lower_endpoint':str(b1),
         'nonlinear_upper_endpoint':str(mass),'circle_caps_numerators':caps,
         'circle_caps_denominator':10**8,'C4_squared_upper':str(q),
         'C4_squared_margin':str(C*C-q),'tail_squared_margin':str(T*T*D-C*C),
         'gamma_lower':str(B-h-T),'gamma_simplified_lower':str(g),
         'strict_reciprocal_margin':str(2*g*K-pi_U),'K_upper':'1.781296297373',
         'head_sha256':sha(work/'head/evidence/independent_head.json'),
         'norm_sha256':sha(work/'norm/evidence/independent_circle.json'),
         'integrals_sha256':sha(work/'norm/evidence/independent_integrals.json')}
    put(work/'final.json',out)
    print('FINAL: K_G^R <= pi/(2*0.8818276493988) < 1.781296297373',flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('stage',choices=['prepare','init','block','merge','head','norm','final','full'])
    ap.add_argument('--workdir',type=Path,required=True)
    ap.add_argument('--block',type=int)
    ap.add_argument('--workers',type=int,default=1)
    ap.add_argument('--allow-long-local-computation',action='store_true',
                    help='Explicitly allow moment blocks without an inherited SLURM allocation.')
    args=ap.parse_args();work=args.workdir.resolve()
    assert 1<=args.workers<=32
    prepare(work)
    cp=work/'head/checks/head501_checkpoint.py'
    def block(j):
        assert 0<=j<32
        extra=[]
        if not os.environ.get('SLURM_JOB_ID'):
            if not args.allow_long_local_computation:
                raise RuntimeError('Moment blocks require an existing allocation or the explicit local-computation option.')
            extra=['--probe'] # Original engine flag; the numerical kernel is identical.
        command(work,f'block_{j:02d}',[cp,'block','--block',str(j),*extra])
    if args.stage in ('init','full'): command(work,'init',[cp,'init'])
    if args.stage=='block':
        if args.block is None: ap.error('--block is required')
        block(args.block)
    if args.stage=='full':
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(block,range(32)))
    if args.stage in ('merge','full'): command(work,'merge',[cp,'merge'])
    if args.stage in ('head','full'):
        command(work,'source_assembly',[work/'head/run/assemble_head.py'])
        command(work,'rational_assembly',[work/'head/checks/independent.py','head'])
    if args.stage in ('norm','full'):
        # A new process creates S_CACHE={} and never reads any integral cache.
        command(work,'norm',[work/'norm/checks/independent.py','tail','--panels','512'])
    if args.stage in ('final','full'): finish(work)
    print('Completed stage',args.stage,flush=True)

if __name__=='__main__': main()
