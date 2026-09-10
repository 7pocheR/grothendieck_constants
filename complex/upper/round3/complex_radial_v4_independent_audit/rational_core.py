"""Exact arithmetic for matrix-derived primitives and short formal series."""
from collections import defaultdict
from fractions import Fraction as Q
from math import comb, isqrt

from audit_exact import require


def specialize(formula, theta, key):
    """Collect a,b monomials only after exact substitution of all parameters."""
    r, alpha, s, beta = key
    parameters = list(map(Q, (theta, r, s, alpha, beta)))
    out = {}
    require(formula["variables"] == ["a","b","h","r","s","alpha","beta"], "Wrong primitive variable order")
    for name, terms in formula["polynomials"].items():
        poly = defaultdict(Q)
        for term in terms:
            powers = term["powers"]
            c = Q(term["coefficient"])
            for x, n in zip(parameters, powers[2:]):
                c *= x**n
            poly[tuple(powers[:2])] += c
        out[name] = {key:c for key,c in poly.items() if c}
    return out


def add(a,b,n):
    return [(a[i] if i<len(a) else Q(0))+(b[i] if i<len(b) else Q(0)) for i in range(n)]


def scale(a,c,n):
    return [x*c for x in a[:n]]+[Q(0)]*max(0,n-len(a))


def mul(a,b,n):
    out=[Q(0)]*n
    for i,x in enumerate(a[:n]):
        for j,y in enumerate(b[:n-i]):
            out[i+j]+=x*y
    return out


def power(a,k,n):
    out=[Q(1)]+[Q(0)]*(n-1)
    for _ in range(k):
        out=mul(out,a,n)
    return out


def inverse(a,n):
    require(a[0]!=0,"Zero constant in reciprocal")
    out=[1/a[0]]
    for k in range(1,n):
        out.append(-sum(a[j]*out[k-j] for j in range(1,min(k+1,len(a))))/a[0])
    return out


def square_root_one(a,n):
    require(a[0]==1,"Square root is not normalized at one")
    out=[Q(1)]
    for k in range(1,n):
        out.append(((a[k] if k<len(a) else 0)-sum(out[j]*out[k-j] for j in range(1,k)))/2)
    return out


def correlation_polynomial(data,key,n):
    out=[Q(0)]*n
    for degree,value in data[key].items():
        degree=int(degree)
        require(degree>0 and degree%2,"Invalid correlation polynomial")
        j=(degree-1)//2
        if j<n:
            out[j]=Q(value)
    return out


def substitute_even_odd(poly,A,B,odd,n):
    out=[Q(0)]*n
    for (i,j),c in poly.items():
        require((i+j)%2==odd,"Unexpected primitive parity")
        shift=(i+j-odd)//2
        term=mul(power(A,i,n),power(B,j,n),n)
        out=add(out,[Q(0)]*shift+scale(term,c,n),n)
    return out


def normalized_coefficients(data,primitives,degree):
    """F(z)/z = pi/(4 sqrt(T0)) times an exact rational series."""
    n=degree+1
    A,B=(correlation_polynomial(data,k,n) for k in ("P","B"))
    D,N1,N2,L,Q0=(substitute_even_odd(primitives[k],A,B,int(k in ("P0","Q0")),n)
                  for k in ("D","N1","N2","P0","Q0"))
    T=mul(N1,N2,n)
    T0=T[0]
    require(T0>0 and D[0]>0,"Invalid origin branch")
    Z=[Q(0)]+mul(mul(L,Q0,n),inverse(T,n),n)[:n-1]
    HG=[Q(1)]+[Q(0)]*(n-1)
    for j in range(1,n):
        c=Q(comb(2*j,j)**2,16**j*(j+1))
        HG=add(HG,scale(power(Z,j,n),c,n),n)
    root=square_root_one(scale(T,1/T0,n),n)
    rational=mul(mul(mul(L,inverse(D,n),n),inverse(root,n),n),HG,n)
    return T0,rational


def atan_bounds(q,terms=120):
    q=Q(q)
    require(0<q<1 and terms%2==0,"Expected positive argument and even term count")
    value=sum((-1)**j*q**(2*j+1)/(2*j+1) for j in range(terms))
    return value,value+q**(2*terms+1)/(2*terms+1)


def pi_bounds():
    a,b=atan_bounds(Q(1,5)),atan_bounds(Q(1,239))
    return 16*a[0]-4*b[1],16*a[1]-4*b[0]


def sqrt_bounds(q,digits=120):
    require(q>0,"Nonpositive square root")
    scale=10**digits
    m=isqrt(q.numerator*scale*scale//q.denominator)
    return Q(m,scale),Q(m+1,scale)


def dyadic(pair):
    require(isinstance(pair,list) and len(pair)==2 and all(type(x) is int for x in pair),"Invalid dyadic endpoint")
    m,e=pair
    require(abs(e)<100000,"Unreasonable dyadic exponent")
    return Q(m)*Q(2)**e


def binary_ball(record):
    require(set(record)=={"mid","rad"},"Invalid binary ball")
    m,r=dyadic(record["mid"]),dyadic(record["rad"])
    require(r>=0,"Negative binary radius")
    return m-r,m+r
