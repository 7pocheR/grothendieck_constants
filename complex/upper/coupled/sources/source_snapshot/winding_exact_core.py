"""Rational construction and interval implications; Python standard library only."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from math import factorial, isqrt
from pathlib import Path
import re


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(path):
    def unique(items):
        result = {}
        for key,value in items:
            require(key not in result, "Duplicate JSON key")
            result[key] = value
        return result
    return json.loads(Path(path).read_text(), object_pairs_hook=unique)


def digest(path):
    result = sha256()
    with Path(path).open("rb") as stream:
        while block := stream.read(1024*1024):
            result.update(block)
    return result.hexdigest()


def modes(data):
    """Expand the combined row/column exponential by the multinomial theorem."""
    generators = []
    for side in range(2):
        for epsilon,frequency,damping in data["profiles"]:
            epsilon,frequency,damping = map(Q,(epsilon,frequency,damping))
            require(damping >= 0, "Negative profile damping")
            for sign in (-1,1):
                generators.append((side,sign*frequency,damping,sign*epsilon/2))
    coefficients = defaultdict(Q)
    k = Q(data["k"])
    for order in range(int(data["phase_order"])+1):
        for choice in combinations_with_replacement(range(len(generators)),order):
            frequencies = [Q(0),Q(0)]
            dampings = [Q(0),Q(0)]
            weight = Q(1)
            for index,multiplicity in Counter(choice).items():
                side,frequency,damping,c = generators[index]
                frequencies[side] += multiplicity*frequency
                dampings[side] += multiplicity*damping
                weight *= c**multiplicity/factorial(multiplicity)
            first,second = sorted((k-frequencies[i],dampings[i]) for i in range(2))
            coefficients[(*first,*second)] += weight
    return sorted((key,value) for key,value in coefficients.items() if value)


def inputs(data):
    for name in ("theta","k","damping_w","damping_x","damping_geometric_lower_bound"):
        Q(data[name])
    dw,dx,lower = (Q(data[n]) for n in
                   ("damping_w","damping_x","damping_geometric_lower_bound"))
    require(dw > 0 and dx > 0 and lower > 0 and dw*dx >= lower*lower,
            "Invalid damping hypotheses")
    norms = {}
    for name in ("P","B"):
        require(all(int(n)>0 and int(n)%2 for n in data[name]), "Nonodd input map")
        norms[name] = sum(abs(Q(c)) for c in data[name].values())
        require(norms[name] <= 1, "Input map exceeds the admissible norm")
    require(data["phase_order"] == 4, "Unexpected phase order")
    amplitude = Q(0)
    for e,l,d in data["profiles"]:
        e,l,d = Q(e),abs(Q(l)),Q(d)
        require(d >= 0, "Negative profile damping")
        amplitude += abs(e)*(Q(1) if d == 0 else min(Q(1),2*l/(5*d*lower)))
    return norms,amplitude,(2*amplitude)**5/factorial(5)


def omission_cost(key,weight,data):
    lower = Q(data["damping_geometric_lower_bound"])
    return abs(weight)/((1+2*key[1]*lower)*(1+2*key[3]*lower))


def specialize(formula,data,key):
    require(formula["variables"] == ["a","b","h","r","s","aw","ax","bw","bx"],
            "Wrong determinant variable order")
    r,u,s,v = key
    dw,dx = Q(data["damping_w"]),Q(data["damping_x"])
    values = [Q(data["theta"]),r,s,u*dw,u*dx,v*dw,v*dx]
    result = {}
    for name,terms in formula["polynomials"].items():
        collected = defaultdict(Q)
        for term in terms:
            c = Q(term["coefficient"])
            for value,power in zip(values,term["powers"][2:]):
                c *= value**power
            collected[tuple(term["powers"][:2])] += c
        result[name] = {powers:c for powers,c in collected.items() if c}
    return result


def dyad(pair):
    require(isinstance(pair,list) and len(pair)==2 and
            all(type(x) is int for x in pair), "Invalid dyadic encoding")
    m,e = pair
    return Q(m*2**e) if e>=0 else Q(m,2**(-e))


def binary_ball(value):
    require(set(value)=={"mid","rad"}, "Invalid binary ball")
    m,r = dyad(value["mid"]),dyad(value["rad"])
    require(r>=0, "Negative ball radius")
    return m-r,m+r


def decimal_ball(text):
    text = text.strip()
    if text.startswith("["):
        match = re.fullmatch(r"\[(.*?)\s*\+/-\s*(.*?)\]",text)
        require(match is not None, "Unsupported decimal interval")
        middle,radius = match.groups()
        m,r = Q(middle or "0"),Q(radius)
        require(r>=0, "Negative decimal radius")
        return m-r,m+r
    q = Q(text)
    return q,q


def floor_dyadic(q,bits=80):
    q = Q(q)
    return [q.numerator*2**bits//q.denominator,-bits]


def ceil_dyadic(q,bits=80):
    q = Q(q)
    return [-((-q.numerator*2**bits)//q.denominator),-bits]


def sqrt_lower(q,bits=96):
    q = Q(q)
    require(q>=0, "Negative square-root input")
    return Q(isqrt(q.numerator*2**(2*bits)//q.denominator),2**bits)


def sqrt_upper(q,bits=96):
    q = Q(q)
    root = sqrt_lower(q,bits)
    return root if root*root==q else root+Q(1,2**bits)


def abs_interval(interval):
    lo,hi = interval
    require(lo<=hi, "Reversed interval")
    return max(Q(0),lo,-hi),max(abs(lo),abs(hi))


def complex_rectangle(value):
    require(set(value)=={"real","imag"}, "Invalid complex rectangle")
    return binary_ball(value["real"]),binary_ball(value["imag"])


def modulus_squared_bounds(value):
    real,imag = complex_rectangle(value)
    ra,rb = abs_interval(real)
    ia,ib = abs_interval(imag)
    return ra*ra+ia*ia,rb*rb+ib*ib


def euler_lower(Z):
    real,imag = complex_rectangle(Z)
    imaginary_lower = abs_interval(imag)[0]
    upper = sqrt_upper(modulus_squared_bounds(Z)[1])
    return max(min(Q(1),1-real[1]),
               imaginary_lower/upper if upper else Q(0))


def multiply_intervals(first,second):
    products = [x*y for x in first for y in second]
    return min(products),max(products)


def quotient_intervals(first,positive):
    lo,hi = positive
    require(lo>0, "Denominator not positive")
    return multiply_intervals(first,(1/hi,1/lo))


def atan_point_bounds(x):
    """Twelve Taylor terms; explicit remainder for |x|<=1/2."""
    x = Q(x)
    require(abs(x)<=Q(1,2), "Argument increment requires refinement")
    value = sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(12))
    error = abs(x)**25/25
    return value-error,value+error


def argument_increment(ratio):
    real,imag = complex_rectangle(ratio)
    lo,hi = quotient_intervals(imag,real)
    lo,hi = dyad(floor_dyadic(lo,40)),dyad(ceil_dyadic(hi,40))
    require(max(abs(lo),abs(hi))<=Q(1,2), "Large endpoint argument increment")
    lower = atan_point_bounds(lo)[0]
    upper = atan_point_bounds(hi)[1]
    return [floor_dyadic(lower,48),ceil_dyadic(upper,48)]


def arc_witness(values,Z,ratios):
    """Exact implications of freshly enclosed rectangles; no displayed strings."""
    ds = [sqrt_lower(modulus_squared_bounds(v)[0]) for v in values[:3]]
    require(all(v>0 for v in ds), "An arc rectangle contains zero")
    delta = dyad(floor_dyadic(euler_lower(Z),80))
    require(delta>0, "No strictly positive Euler distance")
    p = sqrt_upper(modulus_squared_bounds(values[3])[1])
    # The classical elementary inequality pi < 22/7 is sufficient.
    square = (Q(22,7)*p/4)**2/(ds[0]**2*ds[1]*ds[2]*delta)
    M = sqrt_upper(square)
    increments = [argument_increment(ratio) for ratio in ratios]
    return dict(lower=[floor_dyadic(v,96) for v in ds]+[floor_dyadic(delta,80)],
                numerator_upper=ceil_dyadic(p,96),modulus_upper=ceil_dyadic(M,96),
                argument_increments=increments)
