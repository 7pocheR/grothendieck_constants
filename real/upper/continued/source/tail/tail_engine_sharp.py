"""A complete fourth-derivative circle certificate for exact rational data.

Boundary: sum_i B(X_i)+a*sin(omega*sum_i X_i), with identical local B.
Only the tail component is certified. No finite coefficient table is assumed.
"""
from fractions import Fraction
from math import factorial,comb
from pathlib import Path
from itertools import product
from flint import arb,acb,fmpq,fmpq_poly,ctx
import argparse,json,hashlib,time
HERE=Path(__file__).resolve().parent

def saveball(x):
    m,r,e=x.mid_rad_10exp();return {'m':str(m),'r':str(r),'e':int(e)}
def partitions(n,least=1):
    if n==0:yield ()
    for k in range(least,n+1):
        for rest in partitions(n-k,k):yield (k,)+rest
def set_partitions(n):
    if n==0:yield [];return
    for pi in set_partitions(n-1):
        yield pi+[[n-1]]
        for j in range(len(pi)):
            copy=[v[:] for v in pi];copy[j].append(n-1);yield copy

def local_moments(data):
    step=Fraction(data['frequency_step']);cs=[Fraction(v) for v in data['sine_coefficients']];M=len(cs);ell=Fraction(data['linear']);derivs={}
    for r in range(1,5):
        coeff=[Fraction(0)]*(2*M+1)
        for k,c in enumerate(cs,1):
            coeff[M+k]+=c*(step*k)**r/2
            coeff[M-k]+=(-1)**(r+1)*c*(step*k)**r/2
        if r==1:coeff[M]+=ell
        derivs[r]=fmpq_poly([fmpq(str(v)) for v in coeff])
    moments={};exp_cache={};metadata={}
    for total in range(5):
        for part in partitions(total):
            multi=tuple(part.count(r) for r in range(1,5));p=fmpq_poly([1]);count=0;phase=1
            for r,m in enumerate(multi,1):
                p*=derivs[r]**(2*m);count+=2*m;phase*=(-1)**(m*(r-1))
            shift=M*count;value=arb(0);terms=0
            for k,c in enumerate(p):
                if not c:continue
                n=abs(k-shift)
                if n not in exp_cache:exp_cache[n]=(-arb(str(step*step*Fraction(n*n,2)))).exp()
                value+=phase*arb(c)*exp_cache[n];terms+=1
            assert value.is_finite() and value.upper()>=0
            moments[multi]=arb(value.upper());metadata[','.join(map(str,multi))]={'moment':saveball(value),'laurent_terms':terms}
    assert len(moments)==12
    return moments,metadata,len(exp_cache)

def bell_bounds(data,moments):
    a=abs(arb(data['collective_amplitude']));omega=abs(arb(data['collective_frequency']));B={((),0):arb(1)};termcount=0
    for q in range(1,5):
        for alpha in partitions(q):
            if len(alpha)>data['dimension']:continue
            labels=[i for i,n in enumerate(alpha) for _ in range(n)]
            values={k:arb(0) for k in range(1,q+1)}
            for pi in set_partitions(q):
                choices=[]
                for block in pi:
                    r=len(block);coords={labels[j] for j in block};options=[]
                    if len(coords)==1:options.append(('local',next(iter(coords)),r))
                    if a!=0:options.append(('collective',0,r))
                    choices.append(options)
                for selection in product(*choices):
                    coefficient=arb(1);mm=[[0]*4 for _ in alpha]
                    for kind,i,r in selection:
                        if kind=='collective':coefficient*=a*omega**r
                        else:mm[i][r-1]+=1
                    for multi in mm:coefficient*=moments[tuple(multi)].sqrt()
                    values[len(pi)]+=coefficient;termcount+=1
            for k,value in values.items():B[alpha,k]=arb(value.upper())
    return B,termcount

def feature(r,m):
    """Sharp uniform threshold energy on r >= (m-1)/m; monotone extension."""
    assert 0 <= m <= 4 and r >= 0 and r < 1
    if m == 0:return arb(1)
    threshold = (arb(m-1)/m).upper()
    rho = arb(max(r.upper(),threshold))
    den = 1-rho*rho
    assert den > 0
    if m == 1:return 2/(arb.pi()*den.sqrt())
    if m == 2:return 2*rho/(arb.pi()*den*den.sqrt())
    if m == 3:return 2*(1+2*rho*rho)/(arb.pi()*den**2*den.sqrt())
    return 6*rho*(3+2*rho*rho)/(arb.pi()*den**3*den.sqrt())

def polyadd(a,b):
    out=a.copy()
    for key,v in b.items():out[key]=out.get(key,arb(0))+v
    return out
def polymul(a,b):
    out={}
    for (i,j),x in a.items():
        for (k,l),y in b.items():out[i+k,j+l]=out.get((i+k,j+l),arb(0))+x*y
    return out

def density_all(z):
    re=z.real;assert abs(re)<1
    d0=(abs(1-z*z)/(1-re*re)).sqrt();sp=abs(1+z)**2/(1+re);sm=abs(1-z)**2/(1-re)
    L=[{}]
    for j in range(1,5):
        L.append({(0,0):arb(factorial(j-1))/2*(abs(1+z)**(-j)+abs(1-z)**(-j)),(1,0):arb(factorial(j))/2*abs(1+z)**(-j-1),(0,1):arb(factorial(j))/2*abs(1-z)**(-j-1)})
    T=[{(0,0):arb(1)}];D=[d0]
    def moment(s,j):return arb(factorial(2*j))/(2**j*factorial(j))*s**j
    for j in range(1,5):
        poly={}
        for k in range(1,j+1):poly=polyadd(poly,{e:comb(j-1,k-1)*v for e,v in polymul(T[j-k],L[k]).items()})
        T.append(poly);D.append(d0*sum((c*moment(sp,r)*moment(sm,s) for (r,s),c in poly.items()),arb(0)))
    return D

def conv(a,b,N=4):
    out=[a[0]*0]*(N+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b[:N+1-i]):out[i+j]+=x*y
    return out
def jet(ad,bd):
    aa=[acb(0)]+[ad[r]/factorial(r) for r in range(1,5)];bb=[acb(0)]+[bd[r]/factorial(r) for r in range(1,5)]
    ap=[[acb(1)]+[acb(0)]*4];bp=[ap[0]]
    for j in range(4):ap.append(conv(ap[-1],aa));bp.append(conv(bp[-1],bb))
    return {(p,q):arb(24)/factorial(p)/factorial(q)*conv(ap[p],bp[q])[4] for p in range(5) for q in range(5-p) if p+q}

def run(data,panels):
    start=time.time();moments,mm,exp_count=local_moments(data);B,termcount=bell_bounds(data,moments);dimension=data['dimension'];cache={}
    def sbound(ri,p,q):
        key=(ri,p,q)
        if key in cache:return cache[key]
        r=arb(ri)/200;value=arb(0)
        for alpha in partitions(q):
            length=len(alpha)
            if length>dimension:continue
            multiplicity=factorial(dimension)//factorial(dimension-length)
            for part in set(alpha):multiplicity//=factorial(alpha.count(part))
            weight=Fraction(multiplicity*factorial(q))
            for n in alpha:weight/=factorial(n)
            norm=sum((feature(r,p+k).sqrt()*B.get((alpha,k),arb(0)) for k in range(q+1)),arb(0))
            value+=arb(str(weight))*norm*norm
        cache[key]=arb(value.upper());return cache[key]
    polynomials=[]
    for name in ['P','Q']:
        coeff={int(n):Fraction(v) for n,v in data[name].items()}
        assert all(n>=1 and n%2==1 for n in coeff) and sum(map(abs,coeff.values()))<=1
        polynomials.append({n:arb(str(v)) for n,v in coeff.items() if v})
    rows=[];total=arb(0)
    for k in range(panels):
        theta=arb.union(arb.pi()*k/(2*panels),arb.pi()*(k+1)/(2*panels))
        trig={n:acb((n*theta).cos(),(n*theta).sin()) for n in set(polynomials[0])|set(polynomials[1])}
        derivatives=[[sum((v*n**r*trig[n] for n,v in coeff.items()),acb(0)) for r in range(5)] for coeff in polynomials]
        ad,bd=derivatives;pa,pb=ad[0],bd[0];radial=None
        for ri in range(100,195):
            if abs(pa)<arb(ri)/200:radial=ri;break
        da=db=None
        if abs(pa.real)<1 and abs(pb.real)<1:
            da=density_all(pa);single=density_all(pb);aux=[arb(1)]+[arb(0)]*4
            for _ in range(dimension):aux=conv(aux,[single[q]/factorial(q) for q in range(5)])
            db=[aux[q]*factorial(q) for q in range(5)]
        assert radial is not None or da is not None
        major=arb(0)
        for (p,q),coefficient in jet(ad,bd).items():
            choices=[]
            if radial is not None:choices.append(sbound(radial,p,q))
            if da is not None:choices.append(da[p]*db[q])
            major+=abs(coefficient)*arb(min(v.upper() for v in choices))
        major*=arb.pi()/2;total+=major*major
        rows.append({'interval':k,'primary_radius_numerator_over_200':radial,'fourth_derivative_upper':saveball(arb(major.upper()))})
    C4=arb((total/panels).sqrt().upper())
    return {'status':'certified_tail_component_only','precision_bits':ctx.prec,'panels':panels,'local_moments':mm,'distinct_gaussian_exponentials':exp_count,'bell_expansion_terms':termcount,'C4_upper':saveball(C4),'C4_decimal':str(C4),'tails':{str(N):saveball(arb((C4/(arb(14)*arb(N)**7).sqrt()).upper())) for N in [251,501,1001]},'intervals':rows,'seconds':time.time()-start}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('--panels',type=int,default=1024);ap.add_argument('--precision',type=int,default=320);ap.add_argument('--output',default='general_trigonometric_tail.json');args=ap.parse_args();ctx.prec=args.precision;source=HERE/args.input;data=json.loads(source.read_text());result=run(data,args.panels);result['input_sha256']=hashlib.sha256(source.read_bytes()).hexdigest();result['code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();(HERE/args.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['intervals','local_moments']},indent=2))
