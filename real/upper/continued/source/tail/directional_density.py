"""Complete directional Gaussian density derivative norm via Laguerre polynomials."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import comb, factorial
from flint import arb, acb


def aq(x):
    return arb(str(Q(x)))


def convolution(a, b, N):
    result = [acb(0) for _ in range(N+1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b[:N+1-i]):
            result[i+j] += x*y
    return result


def inverse(a, N):
    assert not a[0].contains(0)
    result = [1/a[0]]
    for n in range(1, N+1):
        result.append(-sum((a[j]*result[n-j] for j in range(1, n+1)), acb(0))/a[0])
    return result


def logarithm_positive_degrees(a, N):
    assert not a[0].contains(0)
    y = [acb(0)]+[v/a[0] for v in a[1:]]
    power = [acb(1)]+[acb(0)]*N
    result = [acb(0)]*(N+1)
    for k in range(1, N+1):
        power = convolution(power, y, N)
        for n in range(1, N+1):
            result[n] += Q((-1)**(k+1), k).numerator*power[n]/Q((-1)**(k+1), k).denominator
    return result


def multiply(a, b):
    result = {}
    for u, x in a.items():
        for v, y in b.items():
            key = tuple(i+j for i, j in zip(u, v))
            result[key] = result.get(key, acb(0))+x*y
    return result


def rising(x, n):
    result = Q(1)
    for j in range(n):
        result *= x+j
    return result


@lru_cache(None)
def laguerre_monomial(power, shape):
    # Y=2T, T~Gamma(shape,1), and generalized Laguerre parameter shape-1.
    return tuple(Q(2**power*factorial(power)*(-1)**n, factorial(power-n))*rising(shape+n, power-n)
                 for n in range(power+1))


@lru_cache(None)
def laguerre_tensor(power, dimension):
    shapes = (Q(1, 2), Q(1, 2), Q(dimension, 2), Q(dimension, 2))
    local = [laguerre_monomial(k, a) for k, a in zip(power, shapes)]
    output = []
    for index in product(*(range(k+1) for k in power)):
        coefficient = Q(1)
        for row, n in zip(local, index):
            coefficient *= row[n]
        output.append((index, coefficient))
    return tuple(output)


@lru_cache(None)
def laguerre_norm(index, dimension):
    shapes = (Q(1, 2), Q(1, 2), Q(dimension, 2), Q(dimension, 2))
    result = Q(1)
    for n, a in zip(index, shapes):
        result *= rising(a, n)/factorial(n)
    return result


def density_derivative_norm(ad, bd, dimension, N):
    """Return L1 upper bounds for D^m of the full product density, m<=N.

    ad and bd are D-jets including degree zero, not Taylor coefficients.
    No pi/2 normalization is included. Reflection can be supplied as bd->-bd.
    """
    assert dimension >= 1 and len(ad) >= N+1 and len(bd) >= N+1
    assert abs(ad[0].real) < 1 and abs(bd[0].real) < 1
    correlations = [ad, bd]
    variances, logs, invs = [], [], []
    mass = arb(1)
    for jet, multiplicity in zip(correlations, (1, dimension)):
        z = jet[0]
        mass *= (abs(1-z*z)/(1-z.real*z.real)).sqrt()**multiplicity
        for sign in (1, -1):
            coefficients = [1+sign*z]+[sign*jet[j]/factorial(j) for j in range(1, N+1)]
            logs.append(logarithm_positive_degrees(coefficients, N))
            invs.append(inverse(coefficients, N))
            variances.append(abs(1+sign*z)**2/(1+sign*z.real))
    zero = (0, 0, 0, 0)
    linear = [None]
    for n in range(1, N+1):
        value = {zero: -factorial(n)*(logs[0][n]+logs[1][n]
                                         +dimension*(logs[2][n]+logs[3][n]))/2}
        for j in range(4):
            index = tuple(int(i == j) for i in range(4))
            value[index] = -factorial(n)*invs[j][n]*variances[j]/2
        linear.append(value)
    bell = [{zero: acb(1)}]
    output = [mass]
    for n in range(1, N+1):
        value = {}
        for k in range(1, n+1):
            for index, c in multiply(bell[n-k], linear[k]).items():
                value[index] = value.get(index, acb(0))+comb(n-1, k-1)*c
        bell.append(value)
        orthogonal = {}
        for monomial, c in value.items():
            assert c.is_finite(), (n, monomial)
            for index, coefficient in laguerre_tensor(monomial, dimension):
                orthogonal[index] = orthogonal.get(index, acb(0))+aq(coefficient)*c
        terms = []
        for index, c in orthogonal.items():
            # Ball abs may have a tiny negative lower endpoint when its exact
            # value can vanish. Square its nonnegative upper point directly.
            magnitude = abs(c).upper()
            term = aq(laguerre_norm(index, dimension))*magnitude*magnitude
            assert term.is_finite(), (n, index)
            terms.append(term)
        square = sum(terms, arb(0))
        assert square.is_finite() and square.lower() >= 0, (n, square)
        output.append(mass*square.sqrt())
    return output
