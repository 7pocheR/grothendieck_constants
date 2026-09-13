"""Complete fourth-derivative circle bound from certified mixed-energy endpoints."""
import argparse
from fractions import Fraction as F
import json
from math import factorial
from pathlib import Path
from flint import arb,acb,ctx
import tail_engine_sharp as engine
from conditional_fourier_terms import energy_terms
from fourier_dominant_energy_pilot import require,digest,rat,saveball,loadball,atomic_json
from run_tail_comparison import first_odd_cutoff

HERE=Path(__file__).resolve().parent
PAIRS=[(p,q) for p in range(5) for q in range(5-p) if p+q]


def endpoint_fraction(value):
    m,e=map(int,value.upper().man_exp())
    return F(m)*F(2)**e


def monotone_endpoints(values):
    result=list(values)
    for j in range(len(result)-2,-1,-1):
        result[j]=min(result[j],result[j+1])
    return result


def chord_upper(r,radii,endpoints):
    if r>radii[-1]:
        return None
    if r<=radii[0]:
        return endpoints[0]
    for j in range(1,len(radii)):
        if r<=radii[j]:
            t=(r-radii[j-1])/(radii[j]-radii[j-1])
            return (1-t)*endpoints[j-1]+t*endpoints[j]
    raise AssertionError('Radius was not covered')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--energy-result',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--panels',type=int,default=2048)
    args=parser.parse_args()
    require(args.panels>=64 and not args.output.exists(),'Invalid output or angular partition')
    manifest_path=HERE/'fourier_circle_manifest.json'
    manifest=json.loads(manifest_path.read_text())
    for name,value in manifest['files'].items():
        require(digest(HERE/name)==value,'Frozen source or input changed')
    ctx.prec=384
    data=json.loads((HERE/'inputs/latest_profile_rational.json').read_text())
    result=json.loads(args.energy_result.read_text())
    require(result['status']=='certified_mixed_energy_grid','Input is not a certified grid')
    binding=result['binding']
    require(binding['manifest_sha256']==digest(manifest_path)
            and binding['input_sha256']==digest(HERE/'inputs/latest_profile_rational.json')
            and binding['budget_sha256']==digest(HERE/'fourier_circle_error_budget_final.json'),
            'Mixed-energy source/input bindings differ')
    require(result['radii']==manifest['radii'] and result['pairs']==[list(pair) for pair in PAIRS],
            'Mixed-energy index sets differ')
    require(len(result['rows'])==len(manifest['radii'])*14,'Missing mixed-energy endpoint')
    radii=list(map(F,manifest['radii']))
    advertised={}
    for row in result['rows']:
        key=(F(row['radius']),row['p'],row['q'])
        require(key not in advertised,'Repeated mixed-energy endpoint')
        advertised[key]=row
    analytic=json.loads((HERE/'fourier_circle_error_budget_final.json').read_text())
    error_rows={(F(row['radius']),row['p'],row['q']):row for row in analytic['rows']}
    files={}
    for entry in result['scalar_contributions']:
        relative=Path(entry['path'])
        require(not relative.is_absolute() and '..' not in relative.parts,'Invalid scalar file path')
        require(entry['path'] not in files,'Repeated scalar file path')
        path=args.energy_result.parent/relative
        require(digest(path)==entry['sha256'],'A completed scalar contribution changed')
        files[entry['path']]=path
    require(len(files)==109*len(radii),'Incomplete scalar-function file inventory')
    used=set()
    indexed={}
    for radius in radii:
        groups={}
        totals={pair:arb(0) for pair in PAIRS}
        for p,q in PAIRS:
            for (h,outer,local),coefficient in energy_terms(4,radius,p,q).items():
                groups.setdefault((h,outer),[]).append((p,q,local,coefficient))
        require(len(groups)==109,'Wrong exact scalar-function inventory')
        for (h,outer),terms in groups.items():
            name='radius_'+str(radius.numerator)+'_'+str(radius.denominator)+'/'+str(h)+'_'+'_'.join(map(str,outer))+'.json'
            require(name in files,'Missing scalar function')
            used.add(name)
            stored=json.loads(files[name].read_text())
            require(stored['binding']==binding and stored['radius']==str(radius)
                    and stored['h']==h and stored['outer_exponents']==list(outer),
                    'Scalar function has different inputs')
            signed=[dict(p=p,q=q,local_exponents=[list(row) for row in local],coefficient=str(coefficient))
                    for p,q,local,coefficient in sorted(terms)]
            require(stored['signed_terms']==signed,'Signed scalar terms do not match the exact compiler')
            expected={str((p,q)) for p,q,local,coefficient in terms}
            require(set(stored['finite_contributions'])==expected,'Incomplete scalar contributions')
            for p,q in PAIRS:
                if str((p,q)) in expected:
                    totals[p,q]+=loadball(stored['finite_contributions'][str((p,q))])
        prefactor=2/(arb.pi()*(1-rat(radius)**2).sqrt())
        for p,q in PAIRS:
            key=radius,p,q
            require(key in advertised and key in error_rows,'Missing energy or error allowance')
            error=F(error_rows[key]['analytic_energy_error_upper'])
            require(error>0 and advertised[key]['analytic_error_upper']==str(error),'Omitted-error allowance differs')
            finite=prefactor*totals[p,q]
            value=finite+arb(0,rat(error).upper())
            require(value>0 and value.overlaps(loadball(advertised[key]['energy_enclosure'])),
                    'Recomputed energy is inconsistent with the generator result')
            indexed[key]=endpoint_fraction(value)
    require(used==set(files),'Unrecognized scalar file remains')
    require(set(indexed)=={(r,p,q) for r in radii for p,q in PAIRS},'Wrong endpoint set')
    endpoints={pair:monotone_endpoints([indexed[r,*pair] for r in radii]) for pair in PAIRS}
    # Any returned endpoint may be improved by monotonicity using a later one.
    # No sampled values are used in this interpolation.
    polynomials=[]
    for name in ['P','Q']:
        coefficients={int(n):F(value) for n,value in data[name].items()}
        require(all(n>=1 and n%2 for n in coefficients),'Preprocessing is not odd')
        require(sum(map(abs,coefficients.values()))<=1,'Preprocessing is inadmissible')
        if name=='P':
            require(sum(map(abs,coefficients.values()))<1,'Strict primary disk margin is missing')
        polynomials.append({n:rat(v) for n,v in coefficients.items() if v})
    total=arb(0)
    panels=[]
    count={'interpolated_energy':0,'density':0}
    for k in range(args.panels):
        theta=arb.union(arb.pi()*k/(2*args.panels),arb.pi()*(k+1)/(2*args.panels))
        trig={n:acb((n*theta).cos(),(n*theta).sin()) for n in set(polynomials[0])|set(polynomials[1])}
        ad,bd=[[sum((v*n**r*trig[n] for n,v in polynomial.items()),acb(0)) for r in range(5)]
               for polynomial in polynomials]
        pa,pb=ad[0],bd[0]
        radius=endpoint_fraction(abs(pa))
        da=db=None
        if abs(pa.real)<1 and abs(pb.real)<1:
            da=engine.density_all(pa)
            single=engine.density_all(pb)
            auxiliary=[arb(1)]+[arb(0)]*4
            for _ in range(data['dimension']):
                auxiliary=engine.conv(auxiliary,[single[q]/factorial(q) for q in range(5)])
            db=[auxiliary[q]*factorial(q) for q in range(5)]
        major=arb(0)
        bounds=[]
        for pair,coefficient in engine.jet(ad,bd).items():
            p,q=pair
            choices=[]
            chord=chord_upper(radius,radii,endpoints[pair])
            if chord is not None:
                choices.append(('interpolated_energy',rat(chord)))
            if da is not None:
                choices.append(('density',da[p]*db[q]))
            require(choices,'Uncovered angular interval or mixed derivative')
            name,value=min(choices,key=lambda item:endpoint_fraction(item[1]))
            value=arb(value.upper())
            count[name]+=1
            major+=abs(coefficient)*value
            bounds.append(dict(p=p,q=q,method=name,bound=saveball(value)))
        major*=arb.pi()/2
        total+=major*major
        panels.append(dict(index=k,primary_modulus_upper=str(radius),
                           fourth_derivative_upper=saveball(arb(major.upper())),mixed_bounds=bounds))
    C=endpoint_fraction((total/args.panels).sqrt())
    budgets=[F(1,1000),F(1,2000),F(1,5000),F(1,10000),F(1,20000),F(1,100000)]
    output=dict(status='certified_complete_fourth_derivative_tail_component',
                manifest_sha256=digest(manifest_path),energy_result_sha256=digest(args.energy_result),
                input_sha256=digest(HERE/'inputs/latest_profile_rational.json'),precision_bits=ctx.prec,
                panels=args.panels,exact_C4_upper=str(C),C4_decimal=str(rat(C)),
                radius_interpolation='nonnegative power series; exact convex secant bounds',
                method_counts=count,intervals=panels,
                tails={str(N):saveball(rat(C)/(14*arb(N)**7).sqrt()) for N in [81,251,501,1001]},
                rational_cutoffs=[dict(tail_budget=str(e),first_odd_N=first_odd_cutoff(C,e)) for e in budgets],
                finite_scalar_coefficients_certified=False,Grothendieck_upper_bound_certified=False)
    atomic_json(args.output,output)
    print(json.dumps({key:value for key,value in output.items() if key!='intervals'},indent=2))


if __name__=='__main__':
    main()
