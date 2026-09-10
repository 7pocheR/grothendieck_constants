"""A second complete scalar-tail bound from polynomial norms and mode minima.

This uses no per-arc numerator or modulus bound and no numerical library.
The positive mode minima remain supplied by the verified full covers.
"""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from audit_exact import decimal,load_json,phase_modes,require
from rational_core import dyadic,pi_bounds,specialize,sqrt_bounds


def run(args):
    data=load_json(args.candidate)
    formulas=load_json(args.formulas)
    validation=load_json(args.validation)
    R=Q(103,100)
    A=sum(abs(Q(c))*R**int(n) for n,c in data["P"].items())
    B=sum(abs(Q(c))*R**int(n) for n,c in data["B"].items())
    require(pi_bounds()[1]<Q(22,7),"Rational pi bound")
    factor=R**-963/(1-R**-2)
    tail=Q(0)
    rows=[]
    for i,(key,weight) in enumerate(phase_modes(data,4)):
        if i in (0,165):continue
        record=load_json(Path(args.records)/f"mode_{i:04d}.json")
        D,N1,N2,delta=map(dyadic,record["cover"]["minimum_witnesses"])
        require(min(D,N1,N2,delta)>0,"Nonpositive verified mode minimum")
        poly=specialize(formulas,Q(data["theta"]),key)["P0"]
        numerator=sum(abs(c)*A**j*B**k for (j,k),c in poly.items())
        radicand=N1*N2*delta
        low,_=sqrt_bounds(radicand,digits=35)
        require(low>0 and low*low<=radicand,"Incorrect rational square-root lower bound")
        M=Q(22,7)*numerator/(4*D*low)
        cost=abs(weight)*M*factor
        tail+=cost
        rows.append(dict(index=i,numerator_upper=str(numerator),sqrt_lower=str(low),
                         boundary_modulus_upper=str(M),weighted_scalar_tail_upper=str(cost)))
    gamma=Q(validation["linear_lower"])-Q(validation["finite_nonlinear_upper"])-tail-Q(validation["complete_phase_tail"])-Q(validation["omission_cost"])
    require(len(rows)==208 and gamma>Q(71185999,100000000),"Alternative complete tail does not establish the claimed conservative coefficient bound")
    result=dict(primary_polynomial_modulus_upper=str(A),auxiliary_polynomial_modulus_upper=str(B),
                complete_scalar_tail_upper=str(tail),complete_scalar_tail_upper_decimal=decimal(tail,True),
                gamma_lower=str(gamma),gamma_lower_decimal=decimal(gamma),modes=rows,
                dependency="Exact determinant polynomials and the already verified positive mode minima; no per-arc numerator or modulus fields are used.")
    with open(args.output,"x") as f:
        json.dump(result,f,indent=2);f.write("\n")
    print(json.dumps({k:v for k,v in result.items() if k in ("complete_scalar_tail_upper_decimal","gamma_lower_decimal","dependency")},indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("candidate","formulas","records","validation","output"):
        p.add_argument("--"+name,required=True)
    run(p.parse_args())
