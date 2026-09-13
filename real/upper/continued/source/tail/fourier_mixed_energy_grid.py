"""Enclose all fourteen mixed energies using previously certified local transforms."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
from flint import arb,acb,ctx
from conditional_fourier_terms import energy_terms,serialize_terms
from fourier_dominant_energy_pilot import (require,digest,rat,saveball,loadball,
                                          loadcomplex,atomic_json,negative_dft)

HERE=Path(__file__).resolve().parent
PAIRS=[(p,q) for p in range(5) for q in range(5-p) if p+q]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reuse-dir',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--max-seconds',type=int,default=1200)
    args=parser.parse_args()
    require(args.max_seconds>0,'Positive runtime limit required')
    manifest_path=HERE/'fourier_circle_manifest.json'
    manifest=json.loads(manifest_path.read_text())
    for name,expected in manifest['files'].items():
        require(digest(HERE/name)==expected,'Frozen source or input changed: '+name)
    require(manifest['precision_bits']==256 and manifest['pairs']==[list(pair) for pair in PAIRS],
            'Unexpected mixed-energy settings')
    ctx.prec=manifest['precision_bits']
    input_path=HERE/'inputs/latest_profile_rational.json'
    data=json.loads(input_path.read_text())
    require(data['dimension']==4 and F(data['collective_amplitude'])==0,'Wrong boundary class')
    budget_path=HERE/'fourier_circle_error_budget_final.json'
    budget=json.loads(budget_path.read_text())
    require(budget['input_sha256']==digest(input_path) and budget['status']=='analytic_error_budget_only',
            'Wrong analytic-error input or scope')
    for name,expected in budget['source_sha256'].items():
        require(expected==manifest['files'][name],'Analytic-error source changed')
    expected=dict(L='64',M=8192,N=3072,Delta='pi/32',T='96*pi',hx='1/512',Kx=7168,
                  hs='1/128',Ks=2048,yx='1/8',ys='2/5',y='2/5',beta='1/2',
                  bin_delta='1/250',Taylor_degree=31)
    require(budget['parameters']==expected,'Fourier integration parameters differ')
    radii=list(map(F,manifest['radii']))
    require(radii==[F(n,100) for n in range(81,98,2)],'Wrong radius cover')
    require(len(budget['rows'])==len(radii)*len(PAIRS),'Incomplete analytic error grid')
    budgets={}
    for row in budget['rows']:
        key=(F(row['radius']),row['p'],row['q'])
        require(key not in budgets,'Repeated analytic-error row')
        budgets[key]=row
    require(set(budgets)=={(radius,p,q) for radius in radii for p,q in PAIRS},
            'Wrong analytic error index set')
    reused=json.loads((args.reuse_dir/'result.json').read_text())
    reuse=manifest['reuse']
    require(digest(args.reuse_dir/'result.json')==reuse['result_sha256'],'Reused result changed')
    require(reused['status']=='certified_conditional_energy_component' and reused['binding']==reuse['binding'],
            'Wrong reused computation')
    require(reuse['binding']['input_sha256']==digest(input_path)
            and reuse['binding']['parameters']==expected
            and reuse['binding']['precision_bits']==ctx.prec,'Local transforms have different mathematical inputs')
    local_values={}
    for name,expected_hash in reuse['local_files'].items():
        path=args.reuse_dir/'local_weights'/name
        require(digest(path)==expected_hash,'A reused transform changed: '+name)
        stored=json.loads(path.read_text())
        multi=tuple(stored['exponents'])
        require(name=='_'.join(map(str,multi))+'.json' and multi not in local_values,'Wrong local transform index')
        require(stored['binding']==reuse['binding'] and stored['whole_bin_conditions_checked']==14337,
                'Reused transform has incomplete bin conditions or different bindings')
        require(len(stored['values'])==3073 and stored['finite_frequency_count']==3073,'Incomplete reused transform')
        local_values[multi]=[loadcomplex(value) for value in stored['values']]
    require(len(local_values)==25,'Not all local transforms are present')
    binding=dict(manifest_sha256=digest(manifest_path),input_sha256=digest(input_path),
                 budget_sha256=digest(budget_path),precision_bits=ctx.prec,parameters=expected,
                 reused_result_sha256=reuse['result_sha256'],reused_local_files=reuse['local_files'])
    args.output_dir.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    completed=[]
    energy_rows=[]
    def checkpoint(complete=False):
        atomic_json(args.output_dir/'checkpoint.json',
                    dict(scope='complete finite grid of conditional energies; scalar circle is separate',
                         complete=complete,binding=binding,scalar_files=completed,
                         energy_rows=energy_rows,seconds=time.monotonic()-started))
    def exhausted():
        if time.monotonic()-started>=args.max_seconds:
            checkpoint()
            print('Bounded runtime ended; complete units saved.',flush=True)
            return True
        return False
    N,M,Ks=3072,8192,2048
    hs,Delta=arb(1)/128,arb.pi()/32
    amplitudes=[rat(item['amplitude']) for item in data['outer_profile']]
    frequencies=[rat(item['frequency']) for item in data['outer_profile']]
    s_nodes=[]
    for k in range(Ks+1):
        s=k*hs
        U=s
        derivatives=[arb(1),arb(0),arb(0),arb(0)]
        for coefficient,w in zip(amplitudes,frequencies):
            sn,co=(w*s).sin(),(w*s).cos()
            U+=coefficient*sn
            derivatives[0]+=coefficient*w*co
            derivatives[1]-=coefficient*w*w*sn
            derivatives[2]-=coefficient*w**3*co
            derivatives[3]+=coefficient*w**4*sn
        s_nodes.append((U,derivatives))
    products={}
    def product(local):
        if local not in products:
            require(all(multi in local_values for multi in local),'Missing required local transform')
            values=[acb(1) for _ in range(N+1)]
            for multi in local:
                for n,value in enumerate(local_values[multi]):
                    values[n]*=value
            products[local]=values
        return products[local]
    for radius in radii:
        rdir=args.output_dir/('radius_'+str(radius.numerator)+'_'+str(radius.denominator))
        rdir.mkdir(exist_ok=True)
        groups={}
        for p,q in PAIRS:
            terms=energy_terms(4,radius,p,q)
            simple=[{key:row[key] for key in ['h','outer_exponents','local_exponents','coefficient']}
                    for row in budgets[radius,p,q]['terms']]
            require(simple==serialize_terms(terms),'Finite terms differ from analytic error terms')
            for (h,outer,local),coefficient in terms.items():
                groups.setdefault((h,outer),{}).setdefault((p,q),[]).append((local,coefficient))
        require(len(groups)==109,'Wrong number of scalar functions')
        totals={pair:arb(0) for pair in PAIRS}
        a=1+rat(radius)
        prefactor=2/(arb.pi()*(1-rat(radius)**2).sqrt())
        for (h,outer),pair_terms in sorted(groups.items()):
            label=str(h)+'_'+'_'.join(map(str,outer))
            path=rdir/(label+'.json')
            signed_terms=[dict(p=p,q=q,local_exponents=[list(row) for row in local],coefficient=str(coefficient))
                          for (p,q),terms in sorted(pair_terms.items()) for local,coefficient in sorted(terms)]
            if path.exists():
                stored=json.loads(path.read_text())
                require(stored['binding']==binding and stored['radius']==str(radius)
                        and stored['h']==h and stored['outer_exponents']==list(outer),
                        'Stored scalar function has different inputs')
                require(stored['signed_terms']==signed_terms,'Stored signed terms changed')
                require(set(stored['finite_contributions'])=={str(pair) for pair in pair_terms},
                        'Stored contribution indices are incomplete')
                contributions={pair:loadball(stored['finite_contributions'][str(pair)]) for pair in pair_terms}
            else:
                if exhausted():
                    return
                samples=[arb(0) for _ in range(M)]
                parity=(h+sum((j+2)*exponent for j,exponent in enumerate(outer)))%2
                for k,(U,derivatives) in enumerate(s_nodes):
                    value=hs/(2*arb.pi())*(-U*U/a).exp()*U**h
                    for exponent,derivative in zip(outer,derivatives):
                        value*=derivative**exponent
                    samples[k]=value
                    if k:
                        samples[-k]=(-1)**parity*value
                transformed=negative_dft(samples)
                contributions={}
                for pair,terms in pair_terms.items():
                    combined=[acb(0) for _ in range(N+1)]
                    for local,coefficient in terms:
                        weight=rat(coefficient)
                        for n,value in enumerate(product(local)):
                            combined[n]+=weight*value
                    value=transformed[0]*combined[0]
                    for n in range(1,N+1):
                        value+=2*transformed[n]*combined[n]
                    contributions[pair]=Delta*value.real
                atomic_json(path,dict(binding=binding,radius=str(radius),h=h,outer_exponents=list(outer),
                                      signed_terms=signed_terms,finite_contributions={str(pair):saveball(value)
                                      for pair,value in contributions.items()}))
            for pair,value in contributions.items():
                totals[pair]+=value
            completed.append(dict(path=str(path.relative_to(args.output_dir)),sha256=digest(path)))
            checkpoint()
        radius_rows=[]
        for p,q in PAIRS:
            finite=prefactor*totals[p,q]
            error=F(budgets[radius,p,q]['analytic_energy_error_upper'])
            require(error>0,'Nonpositive omitted-error allowance')
            enclosure=finite+arb(0,rat(error).upper())
            require(enclosure.is_finite() and enclosure>0,'Nonpositive or invalid mixed energy')
            row=dict(radius=str(radius),p=p,q=q,finite_energy=saveball(finite),
                     analytic_error_upper=str(error),energy_enclosure=saveball(enclosure),
                     energy_decimal=str(enclosure))
            energy_rows.append(row)
            radius_rows.append(row)
        atomic_json(rdir/'energy_table.json',dict(binding=binding,complete=True,rows=radius_rows))
        checkpoint()
        print('Completed all mixed energies at radius '+str(radius),flush=True)
    require(len(energy_rows)==126,'Incomplete mixed-energy grid')
    atomic_json(args.output_dir/'result.json',dict(status='certified_mixed_energy_grid',binding=binding,
                radii=list(map(str,radii)),pairs=[list(pair) for pair in PAIRS],rows=energy_rows,
                scalar_contributions=completed,full_fourth_derivative_circle_certified=False,
                finite_scalar_coefficients_certified=False,Grothendieck_upper_bound_certified=False,
                seconds=time.monotonic()-started))
    checkpoint(complete=True)
    print(json.dumps(dict(status='certified_mixed_energy_grid',energy_count=len(energy_rows),
                          scalar_functions=len(completed),seconds=time.monotonic()-started)),flush=True)


if __name__=='__main__':
    main()
