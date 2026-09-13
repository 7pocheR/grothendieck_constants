"""Independent positive-formula replay of the complete head error at higher precision."""
from fractions import Fraction as Q
from math import factorial
from flint import arb,ctx

def A(q):return arb(str(Q(q)))
def up(x):return Q(str(x.upper().fmpq()))

def replay(data,plan,F,bits=512):
    ctx.prec=bits
    parameters=plan['parameters'];N=parameters['scalar_degree'];M=parameters['phase_degree']
    assert (N,M,parameters['dimension'])==(403,31,4)
    row=next(r for r in plan['rows'] if r['frequency_cutoff']==F)
    c=(2/arb.pi()).sqrt();L=A(parameters['period']);Delta=2*arb.pi()/L
    amplitude=[abs(A(r['amplitude'])) for r in data['outer_profile']]
    frequencies=[A(r['frequency']) for r in data['outer_profile']]
    amp=sum(amplitude,arb(0));derivative=sum((a*w for a,w in zip(amplitude,frequencies)),arb(0))
    fac=[factorial(k) for k in range(N+M+1)]
    moments=[arb(1),c]
    for k in range(2,N+M+1):moments.append((k-1)*moments[k-2])
    phase=c*amp**(M+1)/fac[M+1]*(sum((moments[j+M]**2/fac[j] for j in range(N+1)),arb(0))).sqrt()

    def mehler_constants(t,a,da):
        coeff=[a**m/fac[m] for m in range(M+1)]
        hermite=[arb(fac[k]).sqrt()*t**(-A(k)/2) for k in range(N+M+1)]
        multiplier=c*(1-t*t)**(-A(Q(1,4)))
        zero=sum((coeff[m]*hermite[m] for m in range(M+1)),arb(0))
        zero+=da*sum((coeff[m-1]*hermite[m-1] for m in range(1,M+1)),arb(0))
        return [multiplier*zero]+[multiplier/arb(fac[j]).sqrt()*sum((coeff[m]*hermite[j+m-1] for m in range(M+1)),arb(0)) for j in range(1,N+1)]

    t=A(parameters['Mehler_t_real']);beta=A(parameters['alias_beta']);a=1+t
    base=mehler_constants(t,amp,derivative)
    linear=abs(A(data['linear']));sines=[abs(A(q)) for q in data['sine_coefficients']]
    omega=[abs(A(Q(data['frequency_step'])*(i+1))) for i in range(len(sines))]
    variance=4*linear*linear;D=4*sum(sines,arb(0));sigma=variance.sqrt()
    C=beta*D*D/(2*a*(1-beta));cL=beta*L*L/(2*(a+2*beta*variance));mu=(beta*variance/(2*a)).sqrt()
    assert 2*(-3*cL).exp()<1
    aliases=[2*base[0]*(D+sigma+mu*L)*(C-cL).exp()/(1-2*(-3*cL).exp())]
    aliases.extend(2*v*(C-cL).exp()/(1-(-3*cL).exp()) for v in base[1:])
    hx=A(parameters['x_spacing']);yx=A(parameters['strip_height_x']);cut=hx*parameters['x_cutoff']
    T=F*Delta;gx=2*(yx*yx/2).exp()/((2*arb.pi()*yx/hx).exp()-1)
    DN=sum((A(Q(factorial(N),factorial(k)*factorial(N-k)*factorial(k)))*yx**(2*k) for k in range(N+1)),arb(0))
    mass=((1+gx)*(1+gx*DN)).sqrt()
    Bim=linear*yx+sum((v*(w*yx).sinh() for v,w in zip(sines,omega)),arb(0))
    xstrip=2*DN.sqrt()*(yx*yx/2+T*Bim).exp()/((2*arb.pi()*yx/hx).exp()-1)
    xtail=(1-t*t)**(-A(Q(1,4)))*t**(-A(N)/2)*(1+t).sqrt()*(cut/(2*(1+t)).sqrt()).erfc()
    J=parameters['bin_Taylor_degree'];xbin=mass*(T*A(parameters['bin_radius']))**(J+1)/factorial(J+1)
    eta=arb(N+1).sqrt()*(xstrip+xtail+xbin)
    tensor=4*eta+6*eta**2+4*eta**3+eta**4
    harmonic=sum((Q(1,n) for n in range(1,F+1)),Q(0))
    inner=[c*(1+t).sqrt()*base[0]*A(harmonic)*tensor]
    inner.extend(c*Delta*F*(1+t).sqrt()*v*tensor for v in base[1:])
    cache={};freqs=[];choice_checks=[]
    assert [r['j'] for r in row['frequency_choices']]==list(range(N+1))
    for choice in row['frequency_choices']:
        tq,yq=Q(choice['t']),Q(choice['y']);key=tq,yq
        assert 0<tq<1 and yq>0
        if key not in cache:
            tt,y=A(tq),A(yq)
            aa=sum((v*(w*y).cosh() for v,w in zip(amplitude,frequencies)),arb(0))
            dd=sum((v*w*(w*y).cosh() for v,w in zip(amplitude,frequencies)),arb(0))
            factor=(1+tt).sqrt()/2*(y*y/(2*(1-tt))).exp()
            cache[key]=[factor*v for v in mehler_constants(tt,aa,dd)]
        j=choice['j'];y=A(yq)
        fresh=2*c*Delta*cache[key][j]*(-y*(F+1)*Delta).exp()/(1-(-y*Delta).exp())
        if j==0:fresh/=((F+1)*Delta)
        declared=Q(choice['error_upper']);assert up(fresh)<=declared
        freqs.append(A(declared));choice_checks.append({'j':j,'recomputed_upper':str(up(fresh)),'declared_upper':str(declared)})
    integration=sum(((aliases[j]+freqs[j]+inner[j])**2 for j in range(N+1)),arb(0)).sqrt()
    delta=phase+integration;E=arb.pi()/2*delta*(2+delta);claimed=Q(row['scalar_coefficient_error_upper'])
    assert up(E)<=claimed
    return {'status':'PASS_ANALYTIC_ERROR_REPLAY','precision_bits':bits,'frequency_cutoff':F,
            'phase_upper':str(up(phase)),'x_strip_upper':str(up(xstrip)),'x_tail_upper':str(up(xtail)),
            'bin_remainder_upper':str(up(xbin)),'basis_vector_error_upper':str(up(eta)),
            'total_coefficient_array_error_upper':str(up(delta)),
            'recomputed_scalar_error_upper':str(up(E)),'declared_scalar_error_upper':str(claimed),
            'frequency_choices':choice_checks,'passes_arbitrary_planning_target_is_not_required':True}
