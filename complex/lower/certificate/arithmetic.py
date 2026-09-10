"""Directed rational intervals and exponential moments for exact verification."""
from fractions import Fraction as F
from pathlib import Path
from math import factorial
import hashlib
import json

PRECISION = 90
SCALE = 10**PRECISION


class IV:
    """Closed rational interval with endpoints on the grid 10^-PRECISION."""
    __slots__ = ('lo', 'hi')

    def __init__(self, lower=0, upper=None):
        lower = F(lower)
        upper = lower if upper is None else F(upper)
        assert lower <= upper
        self.lo = lower.numerator*SCALE//lower.denominator
        self.hi = -((-upper.numerator*SCALE)//upper.denominator)

    @classmethod
    def raw(cls, lo, hi):
        assert lo <= hi
        result = object.__new__(cls)
        result.lo, result.hi = lo, hi
        return result

    @staticmethod
    def convert(other):
        return other if isinstance(other, IV) else IV(other)

    def __add__(self, other):
        other = self.convert(other)
        return self.raw(self.lo+other.lo, self.hi+other.hi)

    __radd__ = __add__

    def __neg__(self):
        return self.raw(-self.hi, -self.lo)

    def __sub__(self, other):
        return self+(-self.convert(other))

    def __rsub__(self, other):
        return self.convert(other)+(-self)

    def __mul__(self, other):
        other = self.convert(other)
        products = (self.lo*other.lo, self.lo*other.hi,
                    self.hi*other.lo, self.hi*other.hi)
        return self.raw(min(products)//SCALE,
                        -((-max(products))//SCALE))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.convert(other)
        assert other.lo > 0 or other.hi < 0
        reciprocals = [F(SCALE*SCALE,other.lo),
                       F(SCALE*SCALE,other.hi)]
        low, high = min(reciprocals), max(reciprocals)
        reciprocal = self.raw(low.numerator//low.denominator,
                             -((-high.numerator)//high.denominator))
        return self*reciprocal

    def __rtruediv__(self, other):
        return self.convert(other)/self

    def below(self, value):
        return F(self.hi,SCALE) < F(value)

    def above(self, value):
        return F(self.lo,SCALE) > F(value)

    def text(self, places=18):
        unit = 10**(PRECISION-places)
        low = self.lo//unit
        high = -((-self.hi)//unit)
        def decimal(integer):
            sign='-' if integer < 0 else ''
            integer=abs(integer)
            return sign+str(integer//10**places)+'.'+str(integer%10**places).zfill(places)
        return '['+decimal(low)+', '+decimal(high)+']'


def exp_negative_step(step):
    """Alternating Taylor bounds; 0 <= step <= 1 guarantees decreasing terms."""
    assert 0 < step <= 1
    term, partial = F(1), F(1)
    for n in range(1,61):
        term *= -step/n
        partial += term
    upper = partial                 # even partial sum
    lower = partial+term*(-step/61) # odd partial sum
    return IV(lower,upper)


def tail_polynomials(s):
    """Integral_s^infinity t^j exp(-t) dt = exp(-s)*A_j(s)."""
    return [sum(F(factorial(j),factorial(k))*s**k for k in range(j+1))
            for j in range(4)]

