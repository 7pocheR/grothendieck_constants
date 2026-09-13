"""Certified finite arithmetic for one conditional energy, with reusable local transforms."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import time
from flint import arb,acb,ctx
from conditional_fourier_terms import energy_terms,serialize_terms

HERE=Path(__file__).resolve().parent


def require(condition,message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def rat(value):
    return arb(str(F(value)))


def saveball(value):
    require(value.is_finite(),'Nonfinite interval arithmetic')
    m,r,e=value.mid_rad_10exp()
    return dict(m=str(m),r=str(r),e=int(e))


def loadball(value):
    require(set(value)=={'m','r','e'} and type(value['e']) is int,'Malformed stored ball')
    m,r=F(value['m']),F(value['r'])
    require(r>=0,'Negative stored radius')
    return (rat(m)+arb(0,rat(r).upper()))*arb(10)**value['e']


def savecomplex(value):
    return dict(real=saveball(value.real),imag=saveball(value.imag))


def loadcomplex(value):
    require(set(value)=={'real','imag'},'Malformed stored complex ball')
    return acb(loadball(value['real']),loadball(value['imag']))


def atomic_json(path,value):
    temporary=path.with_suffix('.pending.json')
    temporary.write_text(json.dumps(value,separators=(',',':'))+'\n')
    temporary.replace(path)


def positive_dft(real_values):
    return [value.conjugate() for value in acb.dft(real_values)]


def negative_dft(real_values):
    return acb.dft(real_values)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--max-seconds',type=int,default=1200)
    args=parser.parse_args()
    require(args.max_seconds>0,'Positive runtime limit required')
    manifest=json.loads((HERE/'fourier_pilot_manifest.json').read_text())
    for name,expected in manifest['files'].items():
        require(digest(HERE/name)==expected,'Frozen file changed: '+name)
    require(manifest['radius']=='97/100' and manifest['precision_bits']==256,'Wrong frozen pilot')
    ctx.prec=manifest['precision_bits']
    source=HERE/'inputs/latest_profile_rational.json'
    data=json.loads(source.read_text())
    budget_path=HERE/'fourier_error_budget.json'
    budget=json.loads(budget_path.read_text())
    require(budget['status']=='analytic_error_budget_only' and budget['input_sha256']==digest(source),
            'Wrong analytic error budget')
    for name,value in budget['source_sha256'].items():
        require(value==manifest['files'][name],'Error budget source mismatch')
    expected=dict(L='64',M=8192,N=3072,Delta='pi/32',T='96*pi',hx='1/512',Kx=7168,
                  hs='1/128',Ks=2048,yx='1/8',ys='2/5',y='2/5',beta='1/2',
                  bin_delta='1/250',Taylor_degree=31)
    require(budget['parameters']==expected,'Different quadrature or error parameters')
    selected=[row for row in budget['rows'] if row['radius']==manifest['radius']]
    require(len(selected)==1 and selected[0]['p']==0 and selected[0]['q']==4,'Missing dominant-energy budget')
    row=selected[0]
    terms=energy_terms(4,F(manifest['radius']),0,4)
    simple=[{key:item[key] for key in ['h','outer_exponents','local_exponents','coefficient']}
            for item in row['terms']]
    require(simple==serialize_terms(terms),'Signed conditional terms differ from error budget')
    local_weights=sorted({multi for h,outer,local in terms for multi in local})
    require(len(local_weights)==25,'Wrong local-weight inventory')
    M,N,J=expected['M'],expected['N'],expected['Taylor_degree']
    Kx,Ks=expected['Kx'],expected['Ks']
    hx,hs,hb=rat(expected['hx']),rat(expected['hs']),rat(expected['hs'])
    Delta=arb.pi()/32
    delta=rat(expected['bin_delta'])
    frequencies=[Delta*n for n in range(N+1)]
    binding=dict(input_sha256=digest(source),budget_sha256=digest(budget_path),
                 manifest_sha256=digest(HERE/'fourier_pilot_manifest.json'),
                 precision_bits=ctx.prec,parameters=expected)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    local_dir=args.output_dir/'local_weights'
    local_dir.mkdir(exist_ok=True)
    scalar_dir=args.output_dir/'scalar_contributions'
    scalar_dir.mkdir(exist_ok=True)
    started=time.monotonic()
    completed_weights=[]
    completed_scalars=[]
    def checkpoint(complete=False):
        atomic_json(args.output_dir/'checkpoint.json',
                    dict(scope='one conditional Sobolev energy; no full scalar tail or finite head',
                         binding=binding,complete=complete,local_weights=completed_weights,
                         scalar_contributions=completed_scalars,seconds=time.monotonic()-started))
    def exhausted():
        if time.monotonic()-started>=args.max_seconds:
            checkpoint()
            print('Checkpoint saved; bounded runtime ended before completion.',flush=True)
            return True
        return False
    ell=rat(data['linear'])
    cs=[rat(value) for value in data['sine_coefficients']]
    omega=[rat(F(data['frequency_step'])*(j+1)) for j in range(len(cs))]
    amplitudes=[rat(item['amplitude']) for item in data['outer_profile']]
    nu=[rat(item['frequency']) for item in data['outer_profile']]
    require(data['dimension']==4 and F(data['collective_amplitude'])==0,'Wrong boundary class')
    x_nodes=None
    def make_x_nodes():
        output=[]
        for k in range(Kx+1):
            x=k*hx
            B=ell*x
            derivatives=[ell,arb(0),arb(0),arb(0)]
            for coefficient,w in zip(cs,omega):
                sn,co=(w*x).sin(),(w*x).cos()
                B+=coefficient*sn
                derivatives[0]+=coefficient*w*co
                derivatives[1]-=coefficient*w*w*sn
                derivatives[2]-=coefficient*w**3*co
                derivatives[3]+=coefficient*w**4*sn
            mid_m,mid_e=map(int,(B/hb).mid().man_exp())
            middle=F(mid_m)*F(2)**mid_e
            index=(2*middle.numerator+middle.denominator)//(2*middle.denominator)
            residual=B-index*hb
            require(abs(residual)<=delta,'Whole-interval phase-bin condition failed')
            powers=[arb(1)]
            for j in range(1,J+1):
                powers.append(powers[-1]*residual/j)
            gaussian=hx/(2*arb.pi()).sqrt()*(-x*x/2).exp()
            output.append((index,derivatives,powers,gaussian))
        return output
    local_values={}
    for multi in local_weights:
        label='_'.join(map(str,multi))
        path=local_dir/(label+'.json')
        if path.exists():
            stored=json.loads(path.read_text())
            require(stored['binding']==binding and stored['exponents']==list(multi),
                    'Stored local transform has different inputs')
            require(stored['whole_bin_conditions_checked']==2*Kx+1 and len(stored['values'])==N+1,
                    'Incomplete stored local transform')
            values=[loadcomplex(value) for value in stored['values']]
        else:
            if exhausted():
                return
            if x_nodes is None:
                x_nodes=make_x_nodes()
            moments=[[arb(0) for _ in range(M)] for _ in range(J+1)]
            parity=(multi[1]+multi[3])%2
            for k,(index,derivatives,powers,gaussian) in enumerate(x_nodes):
                weight=gaussian
                for exponent,derivative in zip(multi,derivatives):
                    weight*=derivative**exponent
                for j in range(J+1):
                    contribution=weight*powers[j]
                    moments[j][index%M]+=contribution
                    if k:
                        moments[j][(-index)%M]+=(-1)**(parity+j)*contribution
            values=[acb(0) for _ in range(N+1)]
            phase_powers=[acb(1) for _ in range(N+1)]
            for j in range(J+1):
                transformed=positive_dft(moments[j])
                for n in range(N+1):
                    values[n]+=phase_powers[n]*transformed[n]
                    phase_powers[n]*=acb(0,frequencies[n])
            atomic_json(path,dict(binding=binding,exponents=list(multi),
                                 whole_bin_conditions_checked=2*Kx+1,
                                 finite_frequency_count=N+1,values=[savecomplex(v) for v in values]))
        local_values[multi]=values
        completed_weights.append(dict(exponents=list(multi),sha256=digest(path)))
        checkpoint()
        print('Completed local weight '+label,flush=True)
    x_nodes=None
    products={}
    def product(local):
        if local not in products:
            values=[acb(1) for _ in range(N+1)]
            for multi in local:
                for n,value in enumerate(local_values[multi]):
                    values[n]*=value
            products[local]=values
        return products[local]
    groups={}
    for (h,outer,local),coefficient in terms.items():
        groups.setdefault((h,outer),[]).append((local,coefficient))
    require(len(groups)==80,'Wrong scalar-function inventory')
    radius=F(manifest['radius'])
    a=1+rat(radius)
    s_nodes=None
    def make_s_nodes():
        output=[]
        for k in range(Ks+1):
            s=k*hs
            U=s
            derivatives=[arb(1),arb(0),arb(0),arb(0)]
            for coefficient,w in zip(amplitudes,nu):
                sn,co=(w*s).sin(),(w*s).cos()
                U+=coefficient*sn
                derivatives[0]+=coefficient*w*co
                derivatives[1]-=coefficient*w*w*sn
                derivatives[2]-=coefficient*w**3*co
                derivatives[3]+=coefficient*w**4*sn
            output.append((U,derivatives,hs/(2*arb.pi())*(-U*U/a).exp()))
        return output
    contributions=[]
    for h,outer in sorted(groups):
        label=str(h)+'_'+'_'.join(map(str,outer))
        path=scalar_dir/(label+'.json')
        expected_terms=[dict(local_exponents=[list(row) for row in local],coefficient=str(coefficient))
                        for local,coefficient in sorted(groups[h,outer])]
        if path.exists():
            stored=json.loads(path.read_text())
            require(stored['binding']==binding and stored['radius']==str(radius),
                    'Stored scalar contribution has different inputs')
            require(stored['h']==h and stored['outer_exponents']==list(outer)
                    and stored['signed_terms']==expected_terms,'Stored scalar terms changed')
            contribution=loadball(stored['finite_contribution'])
        else:
            if exhausted():
                return
            if s_nodes is None:
                s_nodes=make_s_nodes()
            samples=[arb(0) for _ in range(M)]
            parity=(h+sum((j+2)*exponent for j,exponent in enumerate(outer)))%2
            for k,(U,derivatives,gaussian) in enumerate(s_nodes):
                value=gaussian*U**h
                for exponent,derivative in zip(outer,derivatives):
                    value*=derivative**exponent
                samples[k]=value
                if k:
                    samples[-k]=(-1)**parity*value
            transformed=negative_dft(samples)
            combined=[acb(0) for _ in range(N+1)]
            for local,coefficient in groups[h,outer]:
                weight=rat(coefficient)
                for n,value in enumerate(product(local)):
                    combined[n]+=weight*value
            finite=transformed[0]*combined[0]
            for n in range(1,N+1):
                finite+=2*transformed[n]*combined[n]
            contribution=Delta*finite.real
            atomic_json(path,dict(binding=binding,radius=str(radius),h=h,outer_exponents=list(outer),
                                 signed_terms=expected_terms,finite_contribution=saveball(contribution)))
        contributions.append(contribution)
        completed_scalars.append(dict(label=label,sha256=digest(path)))
        checkpoint()
        print('Completed scalar function '+label,flush=True)
    finite_energy=2/(arb.pi()*(1-rat(radius)**2).sqrt())*sum(contributions,arb(0))
    analytic_error=rat(row['analytic_energy_error_upper'])
    require(analytic_error>0,'Nonpositive omitted-error allowance')
    enclosure=finite_energy+arb(0,analytic_error.upper())
    require(enclosure.is_finite() and enclosure.upper()>0,'Invalid conditional energy enclosure')
    result=dict(status='certified_conditional_energy_component',binding=binding,radius=str(radius),p=0,q=4,
                finite_energy=saveball(finite_energy),analytic_error_upper=row['analytic_energy_error_upper'],
                conditional_energy_enclosure=saveball(enclosure),conditional_energy_decimal=str(enclosure),
                local_weights=completed_weights,scalar_contributions=completed_scalars,
                full_fourth_derivative_circle_certified=False,finite_scalar_coefficients_certified=False,
                Grothendieck_upper_bound_certified=False,seconds=time.monotonic()-started)
    atomic_json(args.output_dir/'result.json',result)
    checkpoint(complete=True)
    print(json.dumps({key:value for key,value in result.items()
                      if key not in ['local_weights','scalar_contributions','binding']},indent=2),flush=True)


if __name__=='__main__':
    main()
