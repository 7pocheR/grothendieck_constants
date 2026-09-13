"""Analytic integration errors for all fourth-order mixed energies on nine radii."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb,factorial
from pathlib import Path
from flint import arb,ctx
from conditional_fourier_terms import energy_terms,serialize_terms

HERE=Path(__file__).resolve().parent


def require(condition,message):
    if not condition:
        raise ValueError(message)


def rational(value):
    return arb(str(F(value)))


def upper_exact(value):
    require(value.is_finite() and value>=0,'Invalid error enclosure')
    mantissa,exponent=map(int,value.upper().man_exp())
    return F(mantissa)*F(2)**exponent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--radii',default='81/100,83/100,17/20,87/100,89/100,91/100,93/100,19/20,97/100')
    args=parser.parse_args()
    require(not args.output.exists(),'Output already exists')
    raw=args.input.read_bytes()
    data=json.loads(raw)
    require(data['dimension']==4 and F(data['collective_amplitude'])==0,'Wrong boundary class')
    ctx.prec=384
    radii=list(map(F,args.radii.split(',')))
    require(radii and sorted(set(radii))==radii and all(0<=r<1 for r in radii),'Invalid radii')
    ell=abs(rational(data['linear']))
    cs=[abs(rational(value)) for value in data['sine_coefficients']]
    omega=[rational(F(data['frequency_step'])*(j+1)) for j in range(len(cs))]
    amplitudes=[abs(rational(row['amplitude'])) for row in data['outer_profile']]
    nu=[abs(rational(row['frequency'])) for row in data['outer_profile']]
    dimension=data['dimension']
    L=arb(64);M=8192;N=3072
    Delta=2*arb.pi()/L;T=N*Delta
    hx=arb(1)/512;Ax=arb(14);hs=arb(1)/128;Z=arb(16)
    yx=arb(1)/8;ys=arb(2)/5;y=ys;beta=arb(1)/2
    delta=arb(1)/250;J=31
    A0=sum(amplitudes,arb(0));CB=sum(cs,arb(0))
    Bim=ell*yx+sum((c*(w*yx).sinh() for c,w in zip(cs,omega)),arb(0))
    gaussian_lattice_mass=1+2*(yx*yx/2).exp()/((2*arb.pi()*yx/hx).exp()-1)
    Vcache={}
    def vdata(multi):
        if multi in Vcache:
            return Vcache[multi]
        W=[]
        for yy in (arb(0),yx):
            value=arb(1)
            for j,exponent in enumerate(multi,1):
                bound=(ell if j==1 else arb(0))+sum((c*w**j*(w*yy).cosh()
                                                    for c,w in zip(cs,omega)),arb(0))
                value*=bound**exponent
            W.append(value)
        ex=2*W[1]*(yx*yx/2+T*Bim).exp()/((2*arb.pi()*yx/hx).exp()-1)
        ex+=W[0]*(Ax/arb(2).sqrt()).erfc()
        eb=W[0]*gaussian_lattice_mass*(T*delta)**(J+1)/factorial(J+1)
        Vcache[multi]=(W[0],ex,eb)
        return Vcache[multi]
    all_rows=[]
    for radius in radii:
        a=1+rational(radius)
        sigma2=dimension*ell*ell
        D=dimension*CB+A0
        c=a/(a+2*beta*sigma2)
        tau=(sigma2*c).sqrt()
        kappa=beta/(a+2*beta*sigma2)
        gcache={}
        def gdata(h,outer):
            key=(h,outer)
            if key in gcache:
                return gcache[key]
            def constant(yy):
                A=sum((aa*(nn*yy).cosh() for aa,nn in zip(amplitudes,nu)),arb(0))
                B=yy+sum((aa*(nn*yy).sinh() for aa,nn in zip(amplitudes,nu)),arb(0))
                Dstrip=A+B
                E=arb(1)
                for j,exponent in enumerate(outer,1):
                    bound=arb(j==1)+sum((aa*nn**j*(nn*yy).cosh()
                                         for aa,nn in zip(amplitudes,nu)),arb(0))
                    E*=bound**exponent
                total=((A+Dstrip)**(h+1)-Dstrip**(h+1))/(h+1)
                for k in range(h+1):
                    moment=a**(arb(k+1)/2)*(arb(k+1)/2).gamma()/2
                    total+=comb(h,k)*(A+Dstrip)**(h-k)*moment
                return E*(B*B/a).exp()*total/arb.pi(),E
            C0,E0=constant(arb(0))
            Cy,_=constant(y)
            require(Z>A0 and Z*Z>=A0*A0+a*h/2,'Scalar lattice tail is not monotone')
            R=Z-A0
            moments=[(arb.pi()*a).sqrt()*(R/a.sqrt()).erfc()/2,
                     a*(-R*R/a).exp()/2]
            for k in range(2,h+1):
                moments.append(a*R**(k-1)*(-R*R/a).exp()/2+a*(k-1)*moments[k-2]/2)
            scalar_tail=E0/arb.pi()*sum((comb(h,k)*(2*A0)**(h-k)*moments[k]
                                        for k in range(h+1)),arb(0))
            scalar_strip=2*Cy*(T*ys).exp()/((2*arb.pi()*ys/hs).exp()-1)
            es=scalar_tail+scalar_strip
            ef=2*Delta*Cy*(-y*(N+1)*Delta).exp()/(1-(-y*Delta).exp())
            rho=arb(2)**h*(-3*kappa*L*L).exp()
            require(rho<1,'Translate ratio is not below one')
            P=arb(0)
            for k in range(h+1):
                absolute_moment=arb(2)**(arb(k)/2)*(arb(k+1)/2).gamma()/arb.pi().sqrt()
                P+=comb(h,k)*(D+c*L)**(h-k)*tau**k*absolute_moment
            et=2*E0*c.sqrt()*(beta*D*D/(a*(1-beta))-kappa*L*L).exp()*P/(1-rho)
            gcache[key]=(C0,es,ef,et)
            return gcache[key]
        for p,q in [(p,q) for p in range(5) for q in range(5-p) if p+q]:
            terms=energy_terms(dimension,radius,p,q)
            total=arb(0)
            component_totals={name:arb(0) for name in ['frequency','translates','finite_factors']}
            term_rows=[]
            for (h,outer,local),coefficient in sorted(terms.items()):
                C0,es,ef,et=gdata(h,outer)
                base=C0
                propagated=es
                W=arb(1)
                for multi in local:
                    wi,ex,eb=vdata(multi)
                    error=ex+eb
                    propagated=propagated*(wi+error)+base*error
                    base*=wi
                    W*=wi
                central=Delta*(2*N+1)*propagated
                errors=dict(frequency=W*ef,translates=W*et,finite_factors=central)
                error=sum(errors.values(),arb(0))
                weight=rational(abs(coefficient))
                total+=weight*error
                for name,value in errors.items():
                    component_totals[name]+=weight*value
                term_rows.append(dict(h=h,outer_exponents=list(outer),
                                      local_exponents=[list(row) for row in local],
                                      coefficient=str(coefficient),
                                      expectation_error_upper=str(upper_exact(error))))
            prefactor=2/(arb.pi()*(1-rational(radius)**2).sqrt())
            total*=prefactor
            row=dict(radius=str(radius),p=p,q=q,term_count=len(terms),
                     unique_scalar_functions=len(gcache),
                     analytic_energy_error_upper=str(upper_exact(total)),
                     analytic_energy_error_decimal=str(total),
                     error_components={name:str(upper_exact(value*prefactor))
                                       for name,value in component_totals.items()},terms=term_rows)
            all_rows.append(row)
    result=dict(status='analytic_error_budget_only',input_sha256=sha256(raw).hexdigest(),
                source_sha256={name:sha256((HERE/name).read_bytes()).hexdigest()
                               for name in ['plan_fourier_circle_enclosure.py','conditional_fourier_terms.py']},
                parameters=dict(L='64',M=M,N=N,Delta='pi/32',T='96*pi',hx='1/512',Kx=7168,
                                hs='1/128',Ks=2048,yx='1/8',ys='2/5',y='2/5',beta='1/2',
                                bin_delta='1/250',Taylor_degree=J),
                rows=all_rows,unique_local_weights=len(Vcache),
                local_weights=[dict(exponents=list(key),supremum_upper=str(upper_exact(value[0])),
                                    trapezoidal_error_upper=str(upper_exact(value[1])),
                                    Taylor_error_upper=str(upper_exact(value[2])))
                               for key,value in sorted(Vcache.items())],
                finite_Fourier_sums_computed=False,nonlinear_energy_certified=False)
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print(json.dumps({key:value for key,value in result.items() if key not in ['rows','local_weights']},indent=2))
    print(json.dumps(dict(energy_rows=len(all_rows), maximum_analytic_error_upper=str(max(F(row['analytic_energy_error_upper']) for row in all_rows))),indent=2))


if __name__=='__main__':
    main()
