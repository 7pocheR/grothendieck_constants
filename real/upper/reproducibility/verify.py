"""Verify the portable degree-501 certificate without recomputing moment blocks.

The default additionally recomputes every smoothing integral from an empty
cache, every circle enclosure, and both complete coefficient assemblies.
Inputs are read-only; all newly generated evidence goes to --output.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, gzip, hashlib, importlib.util, json, shutil, sys, time
from flint import arb, ctx

if sys.flags.optimize:
    raise SystemExit('Assertions must be enabled: do not use python -O.')
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
HEAD = ROOT/'head'
NORM = ROOT/'norm'
ctx.prec = 384

def get(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(v): return arb(v['m']+'e'+str(v['e']), v['r']+'e'+str(v['e']))
def save(v):
    m,r,e=v.mid_rad_10exp()
    return {'m':str(m),'r':str(r),'e':int(e)}
def endpoint(v,side): return (F(v['m'])+side*F(v['r']))*F(10)**v['e']
def upper(v): return endpoint(save(v),1)
def ceil(v): return -(-v.numerator//v.denominator)
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec)
    sys.modules[name]=obj
    spec.loader.exec_module(obj)
    return obj
def write(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2)+'\n')

def integrity():
    manifest=get(ROOT/'immutable_manifest.json')
    actual={str(p.relative_to(ROOT)) for d in (HEAD,NORM)
            for p in d.rglob('*') if p.is_file()}
    actual.add('rational_certificate.json')
    assert actual==set(manifest)
    for p,digest in manifest.items(): assert sha(ROOT/p)==digest, p
    c=get(HEAD/'run/head501_config.json')
    g=get(HEAD/'run/head_moments_501.json')
    assert c['binding']==g['binding']
    b=c['binding']
    assert (b['N'],b['NG'],b['PANELS'],b['block_panels'],b['precision_bits'])==(501,112,128,4,384)
    assert b['candidate_sha256']==sha(HEAD/'sources/candidate_rational.json')==sha(HEAD/'run/candidate_rational.json')==sha(NORM/'sources/candidate_rational.json')
    assert b['head_source_sha256']==sha(HEAD/'sources/certify_head.py')==sha(HEAD/'run/certify_head.py')
    assert b['checkpoint_program_sha256']==sha(HEAD/'checks/head501_checkpoint.py')
    assert c['node_sha256']==sha(HEAD/'run/gauss_nodes_112.json')
    assert g['config_sha256']==sha(HEAD/'run/head501_config.json')
    assert g['uniform_A_error_ball']==save(read(c['uniform_error']))
    h=get(HEAD/'evidence/independent_head.json')
    assert h['moments_sha256']==sha(HEAD/'run/head_moments_501.json')
    assert h['input_sha256']==b['candidate_sha256']
    assert h['checker_sha256']==sha(HEAD/'checks/independent.py')
    for filename in ('independent_circle.json','independent_integrals.json'):
        v=get(NORM/'evidence'/filename)
        assert v['input_sha256']==b['candidate_sha256']
        assert v['checker_sha256']==sha(NORM/'checks/independent.py')
    print('INTEGRITY:',len(manifest),'immutable files and all numerical input bindings',flush=True)
    return {'immutable_files':len(manifest),'candidate_sha256':b['candidate_sha256']}

def merge_check():
    ctx.prec=384
    c=get(HEAD/'run/head501_config.json')
    g=get(HEAD/'run/head_moments_501.json')
    cp=module('portable_checkpoint',HEAD/'checks/head501_checkpoint.py')
    error,bounds=cp.S.ellipse_error()
    saved_error=read(c['uniform_error'])
    assert saved_error.contains(error) and saved_error<arb('4.121e-22')
    assert len(bounds)==len(c['ellipse_M'])==128
    keys={f'{a},{b}' for a in range(502) for b in range(502-a) if (a+b)%2}
    assert len(keys)==63252 and set(g['A'])==keys
    totals={k:arb(0) for k in keys}
    assert len(g['blocks'])==32
    covered=[]
    for j in range(32):
        p=HEAD/'run/blocks'/f'block_{j:02d}.json.gz'
        with gzip.open(p,'rt') as f: part=json.load(f)
        assert g['blocks'][j]['block']==part['block']==j
        assert g['blocks'][j]['sha256']==sha(p)
        assert part['config_sha256']==sha(HEAD/'run/head501_config.json')
        assert part['panels']==list(range(4*j,4*j+4))
        assert part['quadrature_error_added'] is False
        assert set(part['A_partial'])==keys
        for k,v in part['A_partial'].items():
            z=read(v); assert z.is_finite(); totals[k]+=z
        covered.extend(part['panels'])
        if j%8==7: print('BLOCK_SUM:',j+1,'of 32',flush=True)
    assert covered==list(range(128))
    for k,v in totals.items():
        assert save(v+arb(0,saved_error.upper()))==g['A'][k],k
    print('MERGE: all 63252 moment balls reproduced exactly; uniform error replayed',flush=True)
    return {'moment_count':len(keys),'panel_count':128,'uniform_error':save(error)}

def prepare_head(work):
    dest=work/'head'
    for folder in ('sources','run','evidence'): (dest/folder).mkdir(parents=True,exist_ok=True)
    for name in ('candidate_rational.json','certify_head.py','assemble_head.py'):
        shutil.copyfile(HEAD/'sources'/name,dest/'sources'/name)
        shutil.copyfile(HEAD/'run'/name,dest/'run'/name)
    for name in ('gauss_nodes_112.json','head_moments_501.json','forward_head_501.json'):
        shutil.copyfile(HEAD/'run'/name,dest/'run'/name)
    return dest

def head_check(work):
    dest=prepare_head(work)
    source=module('portable_source_assembly',dest/'run/assemble_head.py')
    source.run()
    exactnodes=module('nodes_exact',HEAD/'checks/nodes_exact.py')
    exactnodes.ROOT=dest
    independent=module('portable_head',HEAD/'checks/independent.py')
    independent.ROOT=dest
    independent.algebra()
    independent.head()
    h=get(dest/'evidence/independent_head.json')
    old=get(HEAD/'evidence/independent_head.json')
    odds={str(n) for n in range(1,502,2)}
    assert set(h['coefficients'])==set(old['coefficients'])==odds
    for n in odds: assert read(h['coefficients'][n]).overlaps(read(old['coefficients'][n]))
    mass=sum((max(abs(endpoint(v,-1)),abs(endpoint(v,1)))
              for n,v in h['coefficients'].items() if n!='1'),F(0))
    assert endpoint(h['b1'],-1)>F('0.88182999999999823752')
    assert mass<F('0.000000001054')
    print('HEAD: two complete assemblies, 251 odd coefficients, exact node containment',flush=True)
    return {'odd_coefficients':251,'b1_lower':str(endpoint(h['b1'],-1)),
            'nonlinear_upper':str(mass),'head_sha256':sha(dest/'evidence/independent_head.json')}

PAIRS={f'{p},{q}' for p in range(5) for q in range(5-p) if p+q}
def check_circle_records(circle,integrals):
    assert circle['panels']==512 and circle['order']==4 and len(circle['rows'])==512
    records=integrals['records']
    expected={f'{r}/100,{p},{q}' for r in range(79,98)
              for p in range(5) for q in range(5-p) if p+q}
    assert set(records)==expected and len(records)==266
    for v in records.values():
        assert read(v['integral']).is_finite() and read(v['upper']).is_finite()
        assert read(v['tail'])>=0 and read(v['tail'])<arb(2)**-100
        assert read(v['monotonicity_lhs'])>v['degree']
        assert read(v['leading_coefficient_lower_test'])>0
        assert read(v['lower_terms_margin'])>0
        assert 0<endpoint(v['integral'],-1)<=endpoint(v['upper'],1)
    caps=[]
    for j,row in enumerate(circle['rows']):
        assert row['panel']==j
        assert row['theta_lo_multiple_pi']==f'{j}/1024'
        assert row['theta_hi_multiple_pi']==f'{j+1}/1024'
        assert set(row['terms'])==PAIRS
        major=arb(0)
        for key,term in row['terms'].items():
            if term['method']=='smoothing':
                assert read(row['abs_P_upper'])<arb(row['radial'])
                assert read(term['mixed_upper']).overlaps(read(records[row['radial']+','+key]['upper']))
            else:
                assert term['method']=='density'
                assert read(row['abs_Re_P_upper'])<1 and read(row['abs_Re_z_upper'])<1
            major+=read(term['chain_abs_upper'])*read(term['mixed_upper'])
        major*=arb.pi()/2
        cap=ceil(max(upper(major),endpoint(row['majorant'],1))*10**8)
        caps.append(cap)
    q=sum((F(k,10**8)**2 for k in caps),F(0))/512
    assert q<F('24744.587741')**2
    return caps,q

def norm_check(work):
    dest=work/'norm'
    (dest/'evidence').mkdir(parents=True,exist_ok=True)
    (dest/'sources').mkdir(parents=True,exist_ok=True)
    shutil.copyfile(NORM/'sources/candidate_rational.json',dest/'sources/candidate_rational.json')
    engine=module('portable_norm',NORM/'checks/independent.py')
    engine.ROOT=dest
    assert engine.S_CACHE=={} and engine.S_RECORDS=={}
    chain=engine.algebra()
    engine.circle(chain,512)
    new=get(dest/'evidence/independent_integrals.json')
    old=get(NORM/'evidence/independent_integrals.json')
    for k in new['records']:
        assert read(new['records'][k]['integral']).overlaps(read(old['records'][k]['integral']))
    circle=get(dest/'evidence/independent_circle.json')
    caps,q=check_circle_records(circle,new)
    write(dest/'evidence/rational_norm.json',{'panel_caps_numerators':caps,'denominator':10**8,
          'mean_square_upper':str(q),'C4_squared_margin':str(F('24744.587741')**2-q)})
    print('NORM: all 266 integrals recomputed from an empty cache; 512 whole intervals;',
          'C4 < 24744.587741',flush=True)
    return {'integral_count':266,'panel_count':512,'squared_rational_upper':str(q)}

def final_check():
    cert=get(ROOT/'rational_certificate.json')
    assert cert['head_sha256']==sha(HEAD/'evidence/independent_head.json')
    assert cert['norm_sha256']==sha(NORM/'evidence/independent_circle.json')
    assert cert['integrals_sha256']==sha(NORM/'evidence/independent_integrals.json')
    head=get(HEAD/'evidence/independent_head.json')
    assert set(head['coefficients'])=={str(n) for n in range(1,502,2)}
    assert head['b1']==head['coefficients']['1']
    b1=endpoint(head['b1'],-1)
    mass=sum((max(abs(endpoint(v,-1)),abs(endpoint(v,1)))
              for k,v in head['coefficients'].items() if k!='1'),F(0))
    B=F('0.88182999999999823752');h=F('0.000000001054')
    assert b1>B and mass<h
    caps,q=check_circle_records(get(NORM/'evidence/independent_circle.json'),
                              get(NORM/'evidence/independent_integrals.json'))
    assert len(cert['panel_caps'])==512
    for j,(a,b) in enumerate(zip(caps,cert['panel_caps'])):
        assert b=={'panel':j,'upper_numerator':a,'denominator':10**8}
    assert q==F(cert['C4_squared_rational_upper'])
    C=F('24744.587741');T=F('0.000002349547147305');D=14*501**7
    assert C*C-q==F(cert['C4_squared_rational_margin'])>0
    assert T*T*D-C*C==F(cert['tail_squared_margin'])>0
    gamma=B-h-T;g=F('0.8818276493988');K=F('1.781296297373')
    assert gamma==F(cert['gamma_rational_lower'])>g>0
    assert gamma-g==F(cert['gamma_rational_rounding_margin'])
    def atan(x):
        s=sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(64)),F(0))
        return s,s+x**129/129
    a,b=atan(F(1,5));c,d=atan(F(1,239))
    pi_lo=16*a-4*d;pi_hi=16*b-4*c
    pi_U=F('3.141592653589793238462643383280')
    assert pi_lo==F(cert['pi_lower_machin']) and pi_hi==F(cert['pi_upper_machin'])
    assert pi_hi<pi_U<2*g*K
    assert 2*g*K-pi_U==F(cert['strict_reciprocal_margin_exact'])
    assert K<F('1.78130')<F('1.7813319810625639')<F('1.78184132423347648')
    Ts=F('0.000002578800327634');gs=F('0.8818274201456')
    assert F(27159)**2<Ts**2*D and B-h-Ts>gs
    assert pi_U<2*gs*F('1.781296760466')
    print('EXACT CONCLUSION: K_G^R <= pi/(2*0.8818276493988) < 1.781296297373 < 1.78130',flush=True)
    return {'gamma_lower':str(gamma),'C4_squared_margin':str(C*C-q),
            'tail_squared_margin':str(T*T*D-C*C),'strict_reciprocal_margin':str(2*g*K-pi_U),
            'K_upper':'1.781296297373'}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stage',choices=['all','integrity','merge','head','norm','final'],default='all')
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    work=args.output.resolve()
    assert not work.is_relative_to(ROOT), 'Output must be outside the immutable package.'
    work.mkdir(parents=True,exist_ok=True)
    start=time.time();results={'integrity':integrity()}
    for name,fn in [('merge',merge_check),('head',lambda:head_check(work)),
                    ('norm',lambda:norm_check(work)),('final',final_check)]:
        if args.stage in ('all',name):
            results[name]=fn()
            write(work/'verification.json',results)
    results['seconds']=time.time()-start
    results['stages_completed']=list(results.keys())[:-1]
    write(work/'verification.json',results)
    print('VERIFICATION COMPLETE:',args.stage,'seconds',round(time.time()-start,3),flush=True)

if __name__=='__main__': main()
