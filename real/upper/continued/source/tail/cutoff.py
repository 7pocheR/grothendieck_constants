"""Independent exact power-basis restriction and Bernstein derivative bounds."""
from fractions import Fraction as Q
from math import comb,factorial

def derivative_polynomials():
    p=[Q(1)]+[Q(0)]*13
    for j in range(7):p[7+j]=-Q(12012*(-1)**j*comb(6,j),7+j)
    out=[p]
    for _ in range(6):p=[(i+1)*p[i+1] for i in range(len(p)-1)];out.append(p)
    return out

def evaluate(p,x):
    value=Q(0)
    for c in reversed(p):value=value*x+c
    return value

def restricted_controls(p,lo,hi):
    n=len(p)-1;w=hi-lo
    assert 0<=lo<hi<=1
    # Direct affine substitution in powers, followed by the exact power-to-
    # Bernstein identity; no de Casteljau primitive is reused.
    shifted=[sum((p[i]*comb(i,j)*lo**(i-j)*w**j for i in range(j,n+1)),Q(0)) for j in range(n+1)]
    return [sum((shifted[j]*Q(comb(k,j),comb(n,j)) for j in range(k+1)),Q(0)) for k in range(n+1)]

def bounds(index,a,b,polynomials=None):
    assert type(index) is int and type(a) is int and type(b) is int and 0<=index<2048 and 400<=a<b<2048
    if index+1<=a:return [Q(1)]+[Q(0)]*6,Q(0)
    if index>=b:return [Q(0)]*7,Q(1)
    polynomials=derivative_polynomials() if polynomials is None else polynomials
    lo,hi=Q(index-a,b-a),Q(index+1-a,b-a)
    result=[evaluate(polynomials[0],lo)]
    complement=1-evaluate(polynomials[0],hi)
    assert 0<=result[0]<=1 and 0<=complement<=1
    result.extend(max(map(abs,restricted_controls(p,lo,hi))) for p in polynomials[1:])
    return result,complement
