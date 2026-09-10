"""Independent finite algebra, integral, circle, and coefficient verification.

Inputs are the exact candidate and freshly generated moments/nodes. The
independent tail calculation never reads sobolev_integrals.json. All computed
enclosures carry hashes of this program and the exact candidate.
"""
from pathlib import Path
from math import factorial, comb
from fractions import Fraction
import json, hashlib, sys, time, argparse
from flint import arb, acb, arb_poly, fmpq, fmpq_poly, ctx

ROOT = Path(__file__).resolve().parents[1]
ctx.prec = 384
DATA = json.loads((ROOT/'sources/candidate_rational.json').read_text())
N = 251
XQ = fmpq_poly([0, 1])
X = arb_poly([0, 1])

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def rational(s):
    f=Fraction(s)
    return fmpq(f.numerator,f.denominator)
def ball(v):
    m,r,e=v.mid_rad_10exp()
    return dict(m=str(m), r=str(r), e=int(e))
def read(v): return arb(v['m']+'e'+str(v['e']),v['r']+'e'+str(v['e']))
def endpoint(v, upper):
    return fmpq(int(v['m'])+(int(v['r']) if upper else -int(v['r'])))*fmpq(10)**v['e']
def export(name, data):
    out={'input_sha256':sha(ROOT/'sources/candidate_rational.json'),
         'checker_sha256':sha(Path(__file__)), 'precision_bits':ctx.prec, **data}
    (ROOT/'evidence'/name).write_text(json.dumps(out,indent=2)+'\n')

# Closed factorial formula, with rational coefficients, rather than the
# orthonormal three-term recurrence in the supplied programs.
def heprob(n):
    c=[fmpq(0)]*(n+1)
    for j in range(n//2+1):
        c[n-2*j]=fmpq((-1)**j*factorial(n),2**j*factorial(j)*factorial(n-2*j))
    return fmpq_poly(c)

HE=[heprob(n) for n in range(N+1)]
LP=arb_poly([])
for n,co in zip(DATA['boundary_Hermite_degrees'],DATA['boundary_coefficients']):
    LP += arb_poly([arb(c) for c in HE[n].coeffs()])*(arb(rational(co))/arb(factorial(n)).sqrt())
L1=LP.derivative()
Pq=fmpq_poly([0])
for j,v in enumerate(DATA['W_preprocessing_odd_coefficients']): Pq += rational(v)*XQ**(2*j+1)
P=arb_poly([arb(c) for c in Pq.coeffs()])

def algebra():
    for n in range(1,N): assert HE[n+1]==XQ*HE[n]-n*HE[n-1]
    for n in range(1,N+1): assert HE[n].derivative()==n*HE[n-1]
    mass=sum((abs(c) for c in Pq.coeffs()),fmpq(0))
    assert mass==fmpq(499995394004774337,500000000000000000) and mass<1
    # Exact recurrence for all fourth Euler chain-rule coefficients.
    c={(0,0):fmpq_poly([1])}; dp=XQ*Pq.derivative()
    for s in range(4):
        nxt={}
        for (p,q),v in c.items():
            for key,term in [((p,q),XQ*v.derivative()),((p+1,q),dp*v),((p,q+1),XQ*v)]:
                nxt[key]=nxt.get(key,fmpq_poly([]))+term
        c={key:v for key,v in nxt.items() if v}
    assert len(c)==14
    export('exact_algebra.json',{'hermite_closed_form_recurrence_checks':250,
        'hermite_derivative_checks':251,'P_degree':Pq.degree(),'P_l1':str(mass),
        'L_monomials':{str(k):ball(v) for k,v in enumerate(LP.coeffs())},
        'fourth_euler_coefficients':{f'{p},{q}':[str(t) for t in v.coeffs()] for (p,q),v in c.items()}})
    print('EXACT_ALGEBRA: 250 Hermite recurrences, 251 derivative identities; P l1 =',mass,flush=True)
    return c

# R is a polynomial in T with arb_poly coefficients in x. This constructs
# the derivatives by differential operators, without Bell polynomials or
# the supplied psi_ratio formula.
def tder(R): return [(j+1)*R[j+1] for j in range(len(R)-1)]
def radd(A,B):
    return [(A[j] if j<len(A) else arb_poly([]))+(B[j] if j<len(B) else arb_poly([])) for j in range(max(len(A),len(B)))]
def rt(R,v): return radd(tder(R),[arb_poly([])]+[a*(-1/v) for a in R])
def rx(R,v): return radd([a.derivative() for a in R],[L1*a for a in rt(R,v)])
def smoothing_poly(r,p,q):
    v=1-r
    if p:
        R=[arb_poly([1])]
        for _ in range(p-1): R=rt(R,v)
        for _ in range(q): R=rx(R,v)
    else:
        R=[L1]
        for _ in range(q-1): R=rx(R,v)
    mu=LP*(v/(1+r)); variance=r*v/(1+r)
    M=[]
    for d in range(2*len(R)-1):
        md=arb_poly([])
        for ell in range(d//2+1):
            # binom(d,2ell)*(2ell-1)!! = d!/(2^ell ell! (d-2ell)!).
            co=arb(fmpq(factorial(d),2**ell*factorial(ell)*factorial(d-2*ell)))
            md += mu**(d-2*ell)*(co*variance**ell)
        M.append(md)
    Q=arb_poly([])
    for j,A in enumerate(R):
        for k,B in enumerate(R): Q+=A*B*M[j+k]
    for k in range(1,Q.degree()+1,2): assert Q[k].contains(0)
    return Q

S_CACHE={}; S_RECORDS={}
def smooth(rtext,p,q):
    key=f'{rtext},{p},{q}'
    if key in S_CACHE: return S_CACHE[key]
    r=arb(fmpq(rtext)); Q=smoothing_poly(r,p,q)
    lead=abs(LP[9]); rem=sum((abs(LP[k])*arb(16)**(k-9) for k in range(9)),arb(0))
    assert lead>0 and rem<lead/2
    # A deliberately simple exact rational coefficient bound.
    cc=arb(fmpq(1,1000000))
    assert cc<lead/2
    degree=Q.degree(); monotone=18*cc**2*arb(16)**18/(1+r)
    assert monotone>degree
    pref=2/(arb.pi()*(1-r*r).sqrt())
    tail=pref*sum((abs(a)*arb(16)**k for k,a in enumerate(Q.coeffs())),arb(0))*(-cc**2*arb(16)**18/(1+r)).exp()
    assert tail<arb(2)**-100
    def integrand(z,analytic):
        # Both polynomial evaluations and the exponential are entire.
        ll=LP(z)
        return pref*Q(z)*(-ll*ll/(1+r)-z*z/2).exp()/(2*arb.pi()).sqrt()
    total=acb(0)
    # Explicit independent subdivision, in addition to adaptive certified
    # quadrature internal to each call. The full finite interval is covered.
    for j in range(16):
        total+=acb.integral(integrand,j,j+1,rel_tol=arb('1e-18'),abs_tol=arb('1e-20'),eval_limit=1000000)
    assert total.is_finite() and total.imag.contains(0)
    val=2*total.real+arb(0,tail.upper())
    assert val>0
    result=arb(val.upper())
    S_CACHE[key]=result
    S_RECORDS[key]={'integral':ball(val),'upper':ball(result),'tail':ball(tail),
       'degree':degree,'monotonicity_lhs':ball(monotone),
       'leading_coefficient_lower_test':ball(lead/2-cc),
       'lower_terms_margin':ball(lead/2-rem)}
    return result

# Sparse polynomials in U²,V², generated by logarithmic derivatives of the
# Gaussian density. The original code directly expands variance derivatives.
def padd(A,B):
    C=A.copy()
    for k,v in B.items(): C[k]=C.get(k,acb(0))+v
    return C
def pmul(A,B):
    C={}
    for (a,b),v in A.items():
        for (c,d),w in B.items(): C[a+c,b+d]=C.get((a+c,b+d),acb(0))+v*w
    return C
def pscale(A,s): return {k:v*s for k,v in A.items()}
def dpoly(z,n):
    ell={}
    for k in range(1,n+1):
        vp=1+z; vm=1-z
        # Invert the nonzero variance before taking powers; taking a broad
        # rectangle's power first may produce an interval containing zero.
        ell[k]={(0,0):acb((-1)**k*factorial(k-1))/2*vp**(-k)+acb(factorial(k-1))/2*vm**(-k),
                (1,0):acb((-1)**(k+1)*factorial(k))/2*vp**(-k-1),
                (0,1):-acb(factorial(k))/2*vm**(-k-1)}
    T=[{(0,0):acb(1)}]
    for m in range(1,n+1):
        tm={}
        for k in range(1,m+1): tm=padd(tm,pscale(pmul(ell[k],T[m-k]),comb(m-1,k-1)))
        T.append(tm)
    return T[n]

def density(z,n):
    re=z.real
    assert abs(re)<1
    # Factorized mass formula instead of abs(1-z*z).
    d0=(abs(1+z)*abs(1-z)/(1-re*re)).sqrt()
    T=dpoly(z,n)
    vu=abs(1+z)**2/(1+re); vv=abs(1-z)**2/(1-re)
    # Orthogonal Hermite expansion makes the square norm a sum of squares.
    hc={}
    for (k,m),co in T.items():
        scaled=co*vu**k*vv**m
        for i in range(k+1):
            ci=arb(fmpq(factorial(2*k),2**(k-i)*factorial(k-i)*factorial(2*i)))
            for j in range(m+1):
                cj=arb(fmpq(factorial(2*m),2**(m-j)*factorial(m-j)*factorial(2*j)))
                hc[i,j]=hc.get((i,j),acb(0))+scaled*ci*cj
    norm=sum((abs(v)**2*factorial(2*i)*factorial(2*j) for (i,j),v in hc.items()),arb(0))
    # The exact norm is nonnegative. Ball radii may extend below zero;
    # its positive upper endpoint still bounds the exact squared norm.
    assert norm.is_finite() and norm.upper()>0, (z,n,norm,vu,vv)
    return d0*arb(norm.upper()).sqrt()

def trig_poly(coeff,trig):
    return sum((arb(c)*trig[k] for k,c in enumerate(coeff) if c),acb(0))

def circle(c,panels):
    start=time.time(); rows=[]; square=arb(0)
    pcoef=Pq.coeffs(); ccoef={key:v.coeffs() for key,v in c.items()}
    deg=max(v.degree() for v in c.values())
    for panel in range(panels):
        angle=arb.union(arb.pi()*panel/(2*panels),arb.pi()*(panel+1)/(2*panels))
        trig={k:acb((k*angle).cos(),(k*angle).sin()) for k in range(deg+1)}
        a=trig_poly(pcoef,trig); z=trig[1]
        radial=None
        for ri in range(78,98):
            if abs(a)<arb(fmpq(ri,100)): radial=f'{ri}/100'; break
        da=db=None
        if abs(a.real)<1 and abs(z.real)<1:
            da=[density(a,p) for p in range(5)]
            db=[density(z,q) for q in range(5)]
        assert radial is not None or da is not None
        major=arb(0); terms={}
        for (p,q),cf in ccoef.items():
            bounds=[]
            if radial is not None: bounds.append(('smoothing',smooth(radial,p,q)))
            if da is not None: bounds.append(('density',da[p]*db[q]))
            method,ub=min(bounds,key=lambda t:t[1].upper())
            chain=abs(trig_poly(cf,trig)); term=chain*arb(ub.upper())
            major+=term
            terms[f'{p},{q}']={'method':method,'mixed_upper':ball(arb(ub.upper())),
                              'chain_abs_upper':ball(arb(chain.upper()))}
        major*=arb.pi()/2
        # Each row stores an exact upper endpoint and the accumulation also
        # uses upper endpoints, avoiding unnecessary interval correlations.
        major=arb(major.upper()); square+=major*major
        rows.append({'panel':panel,'theta_lo_multiple_pi':f'{panel}/{2*panels}',
            'theta_hi_multiple_pi':f'{panel+1}/{2*panels}',
            'radial':radial,'abs_P_upper':ball(arb(abs(a).upper())),
            'abs_Re_P_upper':ball(arb(abs(a.real).upper())),
            'abs_Re_z_upper':ball(arb(abs(z.real).upper())),
            'majorant':ball(major),'terms':terms})
        if (panel+1)%64==0:
            export('independent_integrals.json',{'records':S_RECORDS})
            print('INDEPENDENT_CIRCLE',panel+1,'/',panels,'elapsed',round(time.time()-start,2),flush=True)
    norm=arb((square/panels).sqrt().upper())
    tail=norm/(arb(14)*arb(251)**7).sqrt()
    export('independent_circle.json',{'panels':panels,'order':4,'norm_upper':ball(norm),
           'tail251':ball(tail),'rows':rows,'integral_count':len(S_RECORDS)})
    export('independent_integrals.json',{'records':S_RECORDS})
    print('INDEPENDENT_C4',norm,flush=True)
    print('INDEPENDENT_TAIL251',tail,flush=True)

def head():
    g=json.loads((ROOT/'run/gauss_nodes_112.json').read_text())
    le=[fmpq_poly([1]),XQ]
    for n in range(1,112): le.append(((2*n+1)*XQ*le[-1]-n*le[-2])/(n+1))
    L=le[-1]; Ld=L.derivative(); signs=[]; previous=None
    for j,v in enumerate(g['records']):
        lo,hi=endpoint(v['root'],False),endpoint(v['root'],True)
        assert -1<lo<hi<1 and (previous is None or previous<lo)
        previous=hi
        flo=L(lo); fhi=L(hi)
        assert flo*fhi<0
        # Direct derivative polynomial evaluation on the narrow root ball.
        root=read(v['root']); deriv=arb_poly([arb(t) for t in Ld.coeffs()])(root)
        weight=2/((1-root*root)*deriv*deriv)
        assert weight>0 and weight.overlaps(read(v['weight']))
        signs.append({'index':j,'left_sign':-1 if flo<0 else 1,
                     'right_sign':-1 if fhi<0 else 1,'positive_weight':True})
    export('exact_legendre_check.json',{'degree':112,'all_disjoint':True,
                                      'all_in_minus1_1':True,'records':signs})
    print('EXACT_LEGENDRE: 112 disjoint rational root brackets, 224 exact signs, 112 weights checked',flush=True)
    grid=json.loads((ROOT/'run/head_moments_251.json').read_text())
    A={tuple(map(int,k.split(','))):read(v) for k,v in grid['A'].items()}
    assert len(A)==16002
    assert set(A)=={(a,b) for a in range(252) for b in range(252-a) if (a+b)%2}
    # Powers are exact rational polynomials, independently of Arb series.
    powers=[]; power=fmpq_poly([1])
    for a in range(252):
        powers.append([arb(v) for v in power.coeffs()])
        power=(power*Pq).truncate(252)
    coeff=[arb(0) for _ in range(252)]
    for (a,b),ab in sorted(A.items(),key=lambda kv:(kv[0][1],kv[0][0])):
        fac=(-1)**b*ab**2
        for d in range(a,min(len(powers[a]),252-b)):
            if powers[a][d]: coeff[b+d]+=fac*powers[a][d]
    coeff=[v*arb.pi()/2 for v in coeff]
    for n in range(0,252,2): assert coeff[n].is_zero()
    saved=json.loads((ROOT/'run/forward_head_251.json').read_text())
    for n in range(1,252,2): assert coeff[n].overlaps(read(saved['coefficients'][str(n)]))
    nonlinear=sum((abs(coeff[n]) for n in range(3,252,2)),arb(0))
    export('independent_head.json',{'b1':ball(coeff[1]),'nonlinear_head':ball(nonlinear),
        'coefficients':{str(n):ball(coeff[n]) for n in range(1,252,2)},
        'all_126_odd_coefficients_agree':True,'all_even_coefficients_exactly_zero':True})
    print('EXACT_POWER_HEAD b1',coeff[1],'nonlinear',nonlinear,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('mode',choices=['algebra','tail','head']);
    ap.add_argument('--panels',type=int,default=512); args=ap.parse_args()
    chain=algebra()
    if args.mode=='tail': circle(chain,args.panels)
    if args.mode=='head': head()
