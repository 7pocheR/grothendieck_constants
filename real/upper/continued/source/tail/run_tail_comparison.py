"""Complete circle comparison for fixed rational boundaries; no finite head."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from math import comb
from pathlib import Path
import time
from flint import arb,ctx


HERE=Path(__file__).resolve().parent


def require(condition,message):
    if not condition:raise ValueError(message)


def hash_file(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def upper_fraction(ball):
    return F(ball['m'])+F(ball['r']) if ball['e']==0 else (F(ball['m'])+F(ball['r']))*F(10)**ball['e']


def first_odd_cutoff(C,epsilon):
    require(C>0 and epsilon>0,'Positive tail inputs required')
    def good(k):return C*C<=14*(2*k+1)**7*epsilon*epsilon
    upper=1
    while not good(upper):upper*=2
    lower=0
    while lower<upper:
        middle=(lower+upper)//2
        if good(middle):upper=middle
        else:lower=middle+1
    N=2*lower+1
    require(good(lower) and (lower==0 or not good(lower-1)),'Cutoff comparison failed')
    return N


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--feature',choices=['original','sharp'],required=True)
    parser.add_argument('--panels',type=int,default=2048)
    parser.add_argument('--precision',type=int,default=384)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    require(args.panels>0 and args.precision>=128,'Invalid finite cover settings')
    require(not args.output.exists(),'Output already exists')
    ctx.prec=args.precision
    path=HERE/f'tail_engine_{args.feature}.py'
    spec=importlib.util.spec_from_file_location('tail_engine',path)
    base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
    raw=json.loads(args.input.read_text())
    data=dict(raw)
    require(type(data['dimension']) is int and data['dimension']>0,'Invalid dimension')
    require(F(data['collective_amplitude'])==0,'This comparison uses an additive inner boundary')
    local_bounds=base.bell_bounds
    metadata={}
    def profile_bounds(data,moments):
        inner,terms=local_bounds(data,moments)
        profile=[(F(item['amplitude']),F(item['frequency'])) for item in data.get('outer_profile',[])]
        ds=[F(0)]+[F(j==1)+sum(abs(a)*abs(v)**j for a,v in profile) for j in range(1,5)]
        ordinary={(0,0):F(1)}
        for n in range(1,5):
            for k in range(1,n+1):
                ordinary[n,k]=sum(F(comb(n-1,j-1))*ds[j]*ordinary.get((n-j,k-1),F(0))
                                  for j in range(1,n-k+2))
        result={((),0):arb(1)}
        for n in range(1,5):
            for alpha in base.partitions(n):
                if len(alpha)>data['dimension']:continue
                for k in range(1,n+1):
                    value=sum((arb(str(ordinary[r,k]))*inner.get((alpha,r),arb(0))
                               for r in range(k,n+1)),arb(0))
                    result[alpha,k]=arb(value.upper())
        metadata.update(derivative_suprema=list(map(str,ds[1:])),
                        ordinary_Bell_coefficients={str(k):str(v) for k,v in ordinary.items()})
        return result,terms
    base.bell_bounds=profile_bounds
    started=time.monotonic()
    result=base.run(data,args.panels)
    C=upper_fraction(result['C4_upper'])
    budgets=[F(1,1000),F(1,10000),F(1,20000),F(1,100000),F(1,1000000)]
    result.update(exact_C4_upper=str(C),feature_estimate=args.feature,
                  rational_cutoffs=[dict(tail_budget=str(e),first_odd_N=first_odd_cutoff(C,e)) for e in budgets],
                  profile_composition=metadata,
                  source_sha256={path.name:hash_file(path),Path(__file__).name:hash_file(__file__)},
                  input_sha256=hash_file(args.input),
                  finite_coefficient_head_certified=False,
                  Grothendieck_upper_bound_certified=False,
                  total_seconds=time.monotonic()-started)
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['intervals','local_moments','profile_composition','tails']},indent=2))


if __name__=='__main__':main()
