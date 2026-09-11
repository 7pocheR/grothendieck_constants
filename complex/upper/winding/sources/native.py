"""Arb coefficient and complete-circle computation for the Gaussian primitives."""
from fractions import Fraction as Q
import json
from math import comb
from flint import arb,acb,arb_series,acb_poly,ctx
from exact_core import (arc_witness,ceil_dyadic,dyad,floor_dyadic,require)


def rat(value):
    value = Q(value)
    return arb(value.numerator)/value.denominator


def pack(value):
    require(value.is_finite(),"Nonfinite real enclosure")
    return dict(mid=list(map(int,value.mid().man_exp())),
                rad=list(map(int,value.rad().man_exp())))


def pack_complex(value):
    return dict(real=pack(value.real),imag=pack(value.imag))


def pack_geometry_real(value):
    """Outward absolute dyadic endpoint rounding limits certificate size."""
    require(value.is_finite(),"Nonfinite geometric enclosure")
    lo=dyad(list(map(int,value.lower().man_exp())))
    hi=dyad(list(map(int,value.upper().man_exp())))
    left=floor_dyadic(lo,80)[0]
    right=ceil_dyadic(hi,80)[0]
    return dict(mid=[left+right,-81],rad=[right-left,-81])


def pack_geometry_complex(value):
    return dict(real=pack_geometry_real(value.real),imag=pack_geometry_real(value.imag))


def formal_coefficients(data,primitives,degree):
    n = degree+1
    ctx.cap = n
    maps = []
    for name in ("P","B"):
        values = [arb(0)]*n
        for exponent,c in data[name].items():
            j = (int(exponent)-1)//2
            if j<n:
                values[j] = rat(c)
        maps.append(arb_series(values))
    A,B = maps
    y = arb_series([0,1])
    series = {}
    for name,polynomial in primitives.items():
        parity = int(name in ("P0","Q0"))
        value = arb_series([0])
        for (j,k),c in polynomial.items():
            require((j+k)%2==parity,"Incorrect primitive parity")
            value += rat(c)*y**((j+k-parity)//2)*A**j*B**k
        series[name] = value
    D,N1,N2,L,V = (series[name] for name in ("D","N1","N2","P0","Q0"))
    T = N1*N2
    T0 = primitives["N1"][(0,0)]*primitives["N2"][(0,0)]
    require(T0>0 and D[0]>0,"Invalid origin branch")
    Z = y*L*V/T
    require(Z[0]==0,"Hypergeometric argument has nonzero constant")
    H = arb_series([1])
    power = arb_series([1])
    for j in range(1,n):
        power *= Z
        H += rat(Q(comb(2*j,j)**2,16**j*(j+1)))*power
    result = arb.pi()/(4*rat(T0).sqrt())*L/D*H/(T/rat(T0)).sqrt()
    return [pack(result[j]) for j in range(n)]


def compiled(primitives):
    result = []
    for name in ("D","N1","N2","P0","Q0"):
        polynomial = primitives[name]
        n = max(j for j,k in polynomial)+1
        m = max(k for j,k in polynomial)+1
        result.append([[rat(polynomial.get((j,k),0)) for k in range(m)]
                       for j in range(n)])
    return result


def evaluate(polynomials,a,b):
    result = []
    for polynomial in polynomials:
        value = acb(0)
        for row in reversed(polynomial):
            sub = acb(0)
            for c in reversed(row):
                sub = sub*b+c
            value = value*a+sub
        result.append(value)
    return result


def map_polynomials(data,R):
    result = []
    for name in ("P","B"):
        degree = max(map(int,data[name]))
        values = [arb(0)]*((degree+1)//2)
        derivative = values.copy()
        L2 = Q(0)
        for exponent,c in data[name].items():
            exponent,c = int(exponent),Q(c)
            values[(exponent-1)//2] = rat(c)
            derivative[(exponent-1)//2] = rat(exponent*c)
            L2 += exponent**2*abs(c)*R**exponent
        result.append((acb_poly(values),acb_poly(derivative),L2))
    return result


def point(angle,R):
    if angle in (Q(0),Q(2)):
        return acb(rat(R))
    return rat(R)*acb(0,arb.pi()*rat(angle)).exp()


def map_arc(polynomial,left,right,R):
    P,derivative,L2 = polynomial
    z = point((left+right)/2,R)
    center = z*P(z*z)
    first_derivative = acb(0,1)*z*derivative(z*z)
    # Angular displacement is in [-pi*(right-left)/2, +pi*(right-left)/2].
    half_width = Q(22,7)*(right-left)/2
    displacement = arb(0,rat(half_width).upper())
    # Taylor remainder <= (angular displacement)^2 L2 / 2.
    error = half_width**2*L2/2
    remainder = arb(0,rat(error).upper())
    return center+first_derivative*displacement+acb(remainder,remainder)


def geometry(data,primitives,row,settings,stream):
    R = Q(row["radius"])
    maps = map_polynomials(data,R)
    polynomials = compiled(primitives)
    for name in ("D","N1","N2"):
        require(primitives[name][(0,0)]>0,"A primitive is not positive at zero")
    endpoints = {}
    def endpoint(angle):
        if angle==2:
            angle=Q(0)
        if angle not in endpoints:
            z=point(angle,R)
            a,b=(z*P(z*z) for P,derivative,L2 in maps)
            endpoints[angle]=evaluate(polynomials,a,b)
        return endpoints[angle]
    maximum=Q(0)
    sums=[[Q(0),Q(0)] for _ in range(3)]
    minima=[None]*4
    leaves=refinements=0
    panels=settings["panels"]
    for panel in range(panels):
        stack=[(Q(2*panel,panels),Q(2*(panel+1),panels),[])]
        while stack:
            left,right,path=stack.pop()
            a,b=(map_arc(P,left,right,R) for P in maps)
            values=evaluate(polynomials,a,b)
            try:
                require(all(v.is_finite() for v in values),"Nonfinite primitive arc")
                Z=values[3]*values[4]/(values[1]*values[2])
                first,last=endpoint(left),endpoint(right)
                ratios=[b/a for a,b in zip(first[:3],last[:3])]
                encoded=list(map(pack_geometry_complex,values))
                encoded_Z=pack_geometry_complex(Z)
                encoded_ratios=list(map(pack_geometry_complex,ratios))
                witness=arc_witness(encoded,encoded_Z,encoded_ratios)
            except (ValueError,ZeroDivisionError):
                require(len(path)<settings["max_depth"],
                        f"Unresolved complete arc at mode {row['index']}, panel {panel}, path {path}")
                middle=(left+right)/2
                stack.append((middle,right,path+[1]))
                stack.append((left,middle,path+[0]))
                refinements+=1
                continue
            record=dict(panel=panel,path=path,left=str(left),right=str(right),
                        values=encoded,Z=encoded_Z,endpoint_ratios=encoded_ratios,
                        exact_witness=witness)
            stream.write(json.dumps(record,separators=(",",":"))+"\n")
            leaves+=1
            require(leaves<=settings["max_leaves"],"Mode exceeded complete-arc budget")
            maximum=max(maximum,dyad(witness["modulus_upper"]))
            for j,pair in enumerate(witness["lower"]):
                value=dyad(pair)
                minima[j]=value if minima[j] is None else min(minima[j],value)
            for j,(lo,hi) in enumerate(witness["argument_increments"]):
                sums[j][0]+=dyad(lo)
                sums[j][1]+=dyad(hi)
    # pi>3, so this strictly excludes any nonzero multiple of 2*pi.
    require(all(-3<lo<=hi<3 for lo,hi in sums),"Winding number not proved zero")
    require(leaves==panels+refinements,"Incorrect subdivision accounting")
    tail=maximum*R**(-2*row["degree"]-3)/(1-R**(-2))
    return dict(boundary_modulus_upper=str(maximum),complete_scalar_tail_upper=str(tail),
                minima=list(map(str,minima)),argument_sums=[[str(lo),str(hi)] for lo,hi in sums],
                leaves=leaves,refinements=refinements,
                map_second_angular_derivative_bounds=[str(P[2]) for P in maps])


