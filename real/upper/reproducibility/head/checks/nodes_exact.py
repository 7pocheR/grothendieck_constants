"""Exact rational containment proof for all saved Legendre roots and weights."""
from pathlib import Path
import json,hashlib,sys,time
from flint import fmpq,fmpq_poly

ROOT=Path(__file__).resolve().parents[1]
if sys.flags.optimize:raise SystemExit('Python assertions must be enabled.')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def endpoints(v):
    scale=fmpq(10)**v['e'];m=fmpq(v['m']);r=fmpq(v['r'])
    return (m-r)*scale,(m+r)*scale
def square_interval(a,b):
    assert a<=b
    return (fmpq(0) if a<=0<=b else min(a*a,b*b)),max(a*a,b*b)
def certify():
    start=time.time();path=ROOT/'run/gauss_nodes_112.json'
    saved=json.loads(path.read_text());assert saved['degree']==112 and len(saved['records'])==112
    x=fmpq_poly([0,1]);polys=[fmpq_poly([1]),x]
    for k in range(1,112):polys.append(((2*k+1)*x*polys[-1]-k*polys[-2])/(k+1))
    p=polys[-1];dp=p.derivative();lip=fmpq(111*112*113*114,8)
    previous=None;rows=[];slo=fmpq(0);shi=fmpq(0)
    for j,record in enumerate(saved['records']):
        lo,hi=endpoints(record['root']);assert -1<lo<hi<1
        assert previous is None or previous<lo;previous=hi
        flo=p(lo);fhi=p(hi);assert flo*fhi<0
        signs=(-1 if flo<0 else 1,-1 if fhi<0 else 1)
        # Narrow the certified bracket so that an independently computed
        # weight enclosure is strictly contained in the saved source ball.
        for _ in range(8):
            mid=(lo+hi)/2;fm=p(mid)
            if fm*flo>0:lo=mid;flo=fm
            else:
                assert fm*fhi>0;hi=mid;fhi=fm
        mid=(lo+hi)/2;rad=(hi-lo)/2;dmid=dp(mid)
        dlo=dmid-lip*rad;dhi=dmid+lip*rad
        xslo,xshi=square_interval(lo,hi);dslo,dshi=square_interval(dlo,dhi)
        assert dslo>0 and xshi<1
        wlo=2/((1-xslo)*dshi);whi=2/((1-xshi)*dslo)
        savedlo,savedhi=endpoints(record['weight'])
        assert 0<savedlo<wlo<whi<savedhi,(j,'saved weight fails exact containment')
        slo+=savedlo;shi+=savedhi
        rows.append({'index':j,'original_left_sign':signs[0],'original_right_sign':signs[1],
                     'exact_refinement_steps':8,'weight_strict_containment':True})
    assert slo<2<shi
    result={'nodes_sha256':sha(path),'checker_sha256':sha(Path(__file__)),
            'degree':112,'all_112_disjoint_root_brackets_exact':True,
            'all_112_source_weight_balls_rigorously_contain_exact_weights':True,
            'weight_sum_contains_2':True,'arithmetic':'exact rational','records':rows}
    (ROOT/'evidence/exact_node_containment.json').write_text(json.dumps(result,indent=2)+'\n')
    print('EXACT_NODE_CONTAINMENT: 112 disjoint root brackets; 896 exact rational bisections; all 112 source weight enclosures contain independently proved intervals; seconds',round(time.time()-start,3),flush=True)
    return result

if __name__=='__main__':certify()
