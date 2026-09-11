"""Small native smoke test and independent four-term rational series check.

This tests three modes and three small arcs; it is not a full reproduction.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import sys
from flint import arb,ctx
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'sources'))
from exact_core import arc_witness,binary_ball,load,require,specialize
from package_core import check_inputs,compact_records,compact_interval,unoptimized
from native import formal_coefficients,map_polynomials,map_arc,compiled,evaluate,point,pack_geometry_complex
N=4

def add(p,q):
    return [p[i]+q[i] for i in range(N)]


def mul(p,q):
    return [sum(p[j]*q[i-j] for j in range(i+1)) for i in range(N)]


def power(p,k):
    out=[F(1)]+[F(0)]*(N-1)
    for _ in range(k):out=mul(out,p)
    return out


def inverse(p):
    out=[1/p[0]]+[F(0)]*(N-1)
    for k in range(1,N):out[k]=-sum(p[j]*out[k-j] for j in range(1,k+1))/p[0]
    return out


def coefficients(data,formula,key):
    r,u,s,v=map(F,key)
    values=list(map(F,[data['theta'],r,s,u*F(data['damping_w']),u*F(data['damping_x']),v*F(data['damping_w']),v*F(data['damping_x'])]))
    maps=[[F(data[name].get(str(2*i+1),0)) for i in range(N)] for name in ('P','B')]
    series={}
    for name,terms in formula['polynomials'].items():
        parity=int(name in ('P0','Q0'))
        out=[F(0)]*N
        for term in terms:
            j,k,*powers=term['powers']
            if (j+k)%2!=parity:raise ValueError('Parity')
            shift=(j+k-parity)//2
            if shift>=N:continue
            c=F(term['coefficient'])
            for value,exponent in zip(values,powers):c*=value**exponent
            product=mul(power(maps[0],j),power(maps[1],k))
            out=add(out,[F(0)]*shift+[c*x for x in product[:N-shift]])
        series[name]=out
    T=mul(series['N1'],series['N2']);T0=T[0]
    Z=[F(0)]+mul(mul(series['P0'],series['Q0']),inverse(T))[:N-1]
    W=[F(0)]+[x/T0 for x in T[1:]]
    inverse_root=[F(0)]*N;H=[F(0)]*N
    binomial=F(1)
    for j in range(N):
        if j:binomial*=F(1,2)-j;binomial/=j
        inverse_root=add(inverse_root,[binomial*x for x in power(W,j)])
        hyper=F(comb(2*j,j)**2,16**j*(j+1))
        H=add(H,[hyper*x for x in power(Z,j)])
    q=mul(mul(series['P0'],inverse(series['D'])),mul(H,inverse_root))
    return q,T0


def A(q):
    q=F(q);return arb(q.numerator)/q.denominator


def main():
    unoptimized();ctx.prec=768
    data,plan,formula,*_=check_inputs()
    base={**data,'theta':'0','P':{'1':'1'},'B':{'1':'1'}}
    q,T0=coefficients(base,formula,['0','0','0','0'])
    require(T0==1 and q==[F(1),F(1,8),F(3,64),F(25,1024)],'Complex Gaussian normalization')
    chosen=[plan['retained'][0],max(plan['retained'],key=lambda r:r['degree']),plan['retained'][len(plan['retained'])//2]]
    wanted={r['index'] for r in chosen};stored={r['index']:r for r in compact_records() if r['index'] in wanted}
    for row in chosen:
        primitive=specialize(formula,data,list(map(F,row['key'])))
        native=formal_coefficients(data,primitive,N-1)
        q,T0=coefficients(data,formula,row['key']);factor=arb.pi()/(4*A(T0).sqrt())
        for j,c in enumerate(q):
            result=A(c)*factor;lo,hi=binary_ball(native[j])
            lower=F(str(result.lower().fmpq()));upper=F(str(result.upper().fmpq()))
            require(max(lo,lower)<=min(hi,upper),'Independent exact series mismatch')
            packed_lo,packed_hi=compact_interval(stored[row['index']]['coefficients'][j])
            require(packed_lo<=lo<=hi<=packed_hi,'Native coefficient outside compact data')
        R=F(row['radius']);left,right=F(0),F(1,524288)
        maps=map_polynomials(data,R);polynomials=compiled(primitive)
        a,b=(map_arc(p,left,right,R) for p in maps)
        values=evaluate(polynomials,a,b);Z=values[3]*values[4]/(values[1]*values[2])
        endpoints=[]
        for t in (left,right):
            z=point(t,R);a,b=(z*P(z*z) for P,d,L in maps)
            endpoints.append(evaluate(polynomials,a,b))
        ratios=[b/a for a,b in zip(endpoints[0][:3],endpoints[1][:3])]
        arc_witness(list(map(pack_geometry_complex,values)),pack_geometry_complex(Z),list(map(pack_geometry_complex,ratios)))
    print(json.dumps(dict(status='PASS',modes=sorted(wanted),independently_derived_coefficients=12,native_arcs=3,scope='Small smoke test only; no complete circle or scalar degree replay.')))

if __name__=='__main__':main()
