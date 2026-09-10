"""Checkpointed degree-501 Gauss-Legendre moment certificate.

The mathematical kernel uses the unchanged audited source helpers. Each
four-panel block is a separately hashed, resumable interval calculation.
No supplied degree-501 numerical output is read.
"""
from pathlib import Path
import argparse,gzip,hashlib,importlib.util,json,os,sys,time,types
from flint import arb,arb_mat,ctx

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'run'; PARTS=RUN/'blocks'; PARTS.mkdir(exist_ok=True)
N=501; NG=112; PANELS=128; BLOCK_PANELS=4
EXPECTED_CANDIDATE='207b7ed03d6bc1fd44922863e361544bb87442d13801870f0a8bc6705e8a3040'
EXPECTED_SOURCE='55544069a6ef88642f3f7c4b01fa86a8361b45c9678260d5f003d3e3d3e6271a'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(z):
    m,r,e=z.mid_rad_10exp();return {'m':str(m),'r':str(r),'e':int(e)}
def read(z):return arb(z['m']+'e'+str(z['e']),z['r']+'e'+str(z['e']))
def atomic_json(p,obj,compressed=False):
    tmp=p.with_name(p.name+f'.tmp-{os.getpid()}')
    if compressed:
        with gzip.open(tmp,'wt',compresslevel=1) as f:json.dump(obj,f,separators=(',',':'))
    else:tmp.write_text(json.dumps(obj,indent=2)+'\n')
    tmp.replace(p)

assert sha(RUN/'candidate_rational.json')==EXPECTED_CANDIDATE
assert sha(RUN/'certify_head.py')==EXPECTED_SOURCE
assert sha(ROOT/'sources/candidate_rational.json')==EXPECTED_CANDIDATE
assert sha(ROOT/'sources/certify_head.py')==EXPECTED_SOURCE

# SciPy is needed by the source only for untrusted initial root guesses.
# When absent, initialize those guesses from the preserved old node midpoints;
# the source's 112 disjoint sign brackets and bisections are still redone.
try:
    import scipy.special
    INITIALIZER='scipy.special.roots_legendre; all guesses rechecked'
except ModuleNotFoundError as exc:
    if exc.name not in ('scipy','scipy.special'):raise
    scipy=types.ModuleType('scipy');special=types.ModuleType('scipy.special')
    def roots_legendre(n):
        old=json.loads((ROOT/'prior/gauss_nodes_112.json').read_text())
        assert n==112 and len(old['records'])==112
        return [float(read(v['root']).mid()) for v in old['records']],None
    special.roots_legendre=roots_legendre;scipy.special=special
    sys.modules['scipy']=scipy;sys.modules['scipy.special']=special
    INITIALIZER='preserved midpoint float guesses; all guesses rechecked'
spec=importlib.util.spec_from_file_location('unchanged_head',RUN/'certify_head.py')
S=importlib.util.module_from_spec(spec);spec.loader.exec_module(S)
S.N=N;S.NG=NG;ctx.prec=384

def binding():
    return {'candidate_sha256':EXPECTED_CANDIDATE,'head_source_sha256':EXPECTED_SOURCE,
        'checkpoint_program_sha256':sha(Path(__file__)),
        'N':N,'NG':NG,'PANELS':PANELS,'block_panels':BLOCK_PANELS,'precision_bits':ctx.prec}

def initialize():
    target=RUN/'head501_config.json'
    if target.exists():
        config=json.loads(target.read_text())
        assert config['binding']==binding()
        assert config['node_sha256']==sha(RUN/'gauss_nodes_112.json')
        print('INITIALIZATION_REUSED_WITH_MATCHING_HASHES',flush=True)
        return
    start=time.time();S.gauss_nodes();err,bounds=S.ellipse_error()
    assert err<arb('4.121e-22')
    config={'binding':binding(),'node_sha256':sha(RUN/'gauss_nodes_112.json'),
            'uniform_error':save(err),'ellipse_M':bounds,'node_initializer':INITIALIZER,
            'seconds':time.time()-start}
    atomic_json(target,config)
    print('INITIALIZED N',N,'MOMENT_COUNT',63252,'ERROR',err,'seconds',round(time.time()-start,3),flush=True)

def config_read():
    c=json.loads((RUN/'head501_config.json').read_text());assert c['binding']==binding()
    assert c['node_sha256']==sha(RUN/'gauss_nodes_112.json')
    return c

EVENS=list(range(0,N+1,2));ODDS=list(range(1,N+1,2))
KEYS=[(a,b) for a in range(N+1) for b in range(N+1-a) if (a+b)%2]
assert len(KEYS)==63252 and len(ODDS)==251

def block(index,probe=False):
    if not probe and not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Full moment blocks are reserved for root-owned DSI execution. Use --probe for one short local block.')
    assert 0<=index<32
    config=config_read();target=PARTS/f'block_{index:02d}.json.gz'
    config_sha=sha(RUN/'head501_config.json')
    if target.exists():
        with gzip.open(target,'rt') as f:old=json.load(f)
        assert old['config_sha256']==config_sha and old['block']==index
        assert len(old['A_partial'])==63252
        print('BLOCK_REUSED',index,flush=True);return
    start=time.time();nodes=json.loads((RUN/'gauss_nodes_112.json').read_text())
    zz=[read(v['root']) for v in nodes['records']];ww=[read(v['weight']) for v in nodes['records']]
    count=NG*BLOCK_PANELS
    qe=arb_mat(len(EVENS),count);qo=arb_mat(len(ODDS),count)
    he=arb_mat(count,len(EVENS));ho=arb_mat(count,len(ODDS))
    half=arb(16)/(2*PANELS);col=0;begin=index*BLOCK_PANELS
    for panel in range(begin,begin+BLOCK_PANELS):
        center=(2*panel+1)*half
        for z,w in zip(zz,ww):
            x=center+half*z;H=S.hermites(x,N);u=S.Lder(x,0)
            Q=[(u/S.SQ[2]).erf(),S.NU*(-u*u/2).exp()]
            for a in range(1,N):Q.append((-u*S.SQ[a]*Q[a]-(a-1)*Q[a-1])/(S.SQ[a]*S.SQ[a+1]))
            weight=2*half*w*S.PHI*(-x*x/2).exp()
            for j,a in enumerate(EVENS):qe[j,col]=weight*Q[a];he[col,j]=H[a]
            for j,a in enumerate(ODDS):qo[j,col]=weight*Q[a];ho[col,j]=H[a]
            col+=1
    ae=qe*ho;ao=qo*he
    values={}
    for a,b in KEYS:
        v=ae[a//2,b//2] if a%2==0 else ao[a//2,b//2]
        assert v.is_finite();values[f'{a},{b}']=save(v)
    out={'config_sha256':config_sha,'block':index,'panels':list(range(begin,begin+BLOCK_PANELS)),
         'quadrature_error_added':False,'A_partial':values,'seconds':time.time()-start,
         'execution':{'slurm_job_id':os.environ.get('SLURM_JOB_ID'),'local_probe':bool(probe)}}
    atomic_json(target,out,compressed=True)
    print('BLOCK_COMPLETE',index,'panels',begin,'through',begin+BLOCK_PANELS-1,
          'moments',len(values),'seconds',round(time.time()-start,3),flush=True)

def merge():
    config=config_read();config_sha=sha(RUN/'head501_config.json')
    values={f'{a},{b}':arb(0) for a,b in KEYS};seen=[];records=[]
    for index in range(32):
        p=PARTS/f'block_{index:02d}.json.gz'
        if not p.exists():raise RuntimeError(f'Missing block {index}; the full degree-501 moment certificate is incomplete')
        with gzip.open(p,'rt') as f:b=json.load(f)
        assert b['config_sha256']==config_sha and b['block']==index
        assert b['panels']==list(range(index*4,index*4+4))
        assert not b['quadrature_error_added'] and set(b['A_partial'])==set(values)
        for key,v in b['A_partial'].items():values[key]+=read(v)
        seen+=b['panels'];records.append({'block':index,'sha256':sha(p),'execution':b['execution']})
    assert seen==list(range(128))
    err=read(config['uniform_error']);assert err<arb('4.121e-22')
    out={'N':N,'NG':NG,'PANELS':PANELS,'precision_bits':ctx.prec,
        'uniform_A_error':str(err),'uniform_A_error_ball':save(err),'ellipse_M':config['ellipse_M'],
        'A':{key:save(v+arb(0,err.upper())) for key,v in values.items()},
        'binding':binding(),'config_sha256':config_sha,'blocks':records,
        'all_128_panels_covered_once':True,'moment_count':63252}
    atomic_json(RUN/'head_moments_501.json',out)
    print('MERGE_COMPLETE: all 63252 moments certified; every panel used exactly once',flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['init','block','merge'])
    ap.add_argument('--block',type=int);ap.add_argument('--probe',action='store_true');args=ap.parse_args()
    if args.mode=='init':initialize()
    elif args.mode=='block':
        if args.block is None:ap.error('--block is required')
        block(args.block,args.probe)
    else:merge()
