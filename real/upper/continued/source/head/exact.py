"""Exact algebra shared by the staged replay and seconds-long synthetic checks."""
from fractions import Fraction as Q
from math import comb, factorial, lcm


def require(ok, message):
    if not ok:
        raise ValueError(message)


def mul(a, b, size=None):
    size = len(a) + len(b) - 1 if size is None else size
    c = [0] * max(0, size)
    for i, x in enumerate(a):
        for j, y in enumerate(b[:max(0, size-i)]):
            c[i+j] += x*y
    return c


def power(a, n, size=None):
    r = [1]
    while n:
        if n & 1:
            r = mul(r, a, size)
        n //= 2
        if n:
            a = mul(a, a, size)
    return r if size is None else (r+[0]*size)[:size]


def hermite_coefficients(k):
    c = [0] * (k+1)
    for a in range(k//2+1):
        c[k-2*a] = (-1)**a * factorial(k)//(2**a*factorial(a)*factorial(k-2*a))
    return c


def horner(c, x):
    out = 0
    for v in reversed(c):
        out = out*x+v
    return out


def laurent(data, degree, polynomial_type=None, tick=lambda: None):
    """Clear denominators; independently form each integer polynomial power.

    z^s R(z)=A(z)/D. Then r[m,k]=[z^(k+s*m)]A(z)^m/(D^m m!).
    The optional FLINT integer polynomial backend changes only exact arithmetic.
    No producer's repeated R/m recurrence is used.
    """
    terms = {}
    for row in data['outer_profile']:
        k, c = Q(row['frequency'])/2, Q(row['amplitude'])/2
        require(k.denominator == 1 and k > 0, 'Invalid profile frequency')
        k = int(k)
        terms[k] = terms.get(k, Q(0))+c
        terms[-k] = terms.get(-k, Q(0))-c
    s = max(map(abs, terms), default=0)
    denominator = lcm(*(c.denominator for c in terms.values())) if terms else 1
    base = [int(denominator*terms.get(k-s, 0)) for k in range(2*s+1)]
    rows = []
    for m in range(degree+1):
        tick()
        if polynomial_type is None:
            values = power(base, m)
        else:
            poly = polynomial_type(base)**m
            values = [int(poly[i]) for i in range(2*s*m+1)]
        divisor = denominator**m*factorial(m)
        row = {i-s*m: Q(v, divisor) for i, v in enumerate(values) if v}
        require(all(row.get(-k, 0) == (-1)**m*v for k, v in row.items()),
                'Laurent parity failed')
        require(sum(row.values(), Q(0)) == (1 if m == 0 else 0),
                'R(1)=0 identity failed')
        rows.append(row)
    return rows


def factored_laurent(data, degree):
    """Small-check oracle from the product of harmonic exponential factors."""
    terms = {(0, 0): Q(1)}
    for row in data['outer_profile']:
        k, c = int(Q(row['frequency'])/2), Q(row['amplitude'])/2
        factor = {(m, k*(m-2*b)): c**m*(-1)**b/Q(factorial(b)*factorial(m-b))
                  for m in range(degree+1) for b in range(m+1)}
        new = {}
        for (m, s), a in terms.items():
            for (n, t), b in factor.items():
                if m+n <= degree:
                    new[m+n, s+t] = new.get((m+n, s+t), Q(0))+a*b
        terms = {key: v for key, v in new.items() if v}
    return [{k: v for (n, k), v in terms.items() if n == m} for m in range(degree+1)]


def pairs(nodes):
    return [(u, v) for i, u in enumerate(nodes) for v in nodes[i:]]


def prefix_count(row, next_v, cutoff):
    require(1 <= row <= cutoff and row <= next_v <= cutoff+1, 'Invalid pair cursor')
    return next_v-row
