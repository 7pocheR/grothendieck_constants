"""Exact rational expansion of the conditional Gaussian energy integrands."""
from fractions import Fraction as F
from math import comb,factorial


def partitions(n,least=1):
    if n==0:
        yield ()
    for k in range(least,n+1):
        for tail in partitions(n-k,k):
            yield (k,)+tail


def set_partitions(n):
    if n==0:
        yield []
        return
    for previous in set_partitions(n-1):
        yield previous+[[n-1]]
        for j in range(len(previous)):
            result=[list(block) for block in previous]
            result[j].append(n-1)
            yield result


def addpoly(destination,source,factor=F(1)):
    for key,value in source.items():
        destination[key]=destination.get(key,F(0))+factor*value
        if destination[key]==0:
            del destination[key]


def multiply_univariate(left,right):
    result={}
    for i,a in left.items():
        for j,b in right.items():
            result[i+j]=result.get(i+j,F(0))+a*b
    return {key:value for key,value in result.items() if value}


def conditional_gram(radius):
    r=F(radius)
    if not 0<=r<1:
        raise ValueError('Radius outside [0,1)')
    v,a=1-r,1+r
    mu={1:v/a}
    variance=r*v/a
    moments=[{0:F(1)},mu]
    for j in range(1,6):
        value=multiply_univariate(mu,moments[-1])
        addpoly(value,moments[-2],j*variance)
        moments.append(value)
    J={1:{0:F(1)},2:{1:-1/v},3:{2:1/v**2,0:-1/v},
       4:{3:-1/v**3,1:3/v**2}}
    gram={}
    for m in range(1,5):
        for n in range(1,5):
            value={}
            for power,coefficient in multiply_univariate(J[m],J[n]).items():
                addpoly(value,moments[power],coefficient)
            gram[m,n]=value
    return gram


def ordinary_bell():
    zero=(0,0,0,0)
    result={(0,0):{zero:F(1)}}
    for n in range(1,5):
        for k in range(1,n+1):
            value={}
            for j in range(1,n-k+2):
                for key,coefficient in result.get((n-j,k-1),{}).items():
                    new=list(key)
                    new[j-1]+=1
                    new=tuple(new)
                    value[new]=value.get(new,F(0))+comb(n-1,j-1)*coefficient
            result[n,k]=value
    return result


def boundary_bell(alpha,dimension,ordinary):
    zero=(0,0,0,0)
    q=sum(alpha)
    if q==0:
        return {0:{(zero,(zero,)*dimension):F(1)}}
    labels=[coordinate for coordinate,count in enumerate(alpha) for _ in range(count)]
    inner={}
    for blocks in set_partitions(q):
        local=[[0]*4 for _ in range(dimension)]
        for block in blocks:
            coordinates={labels[index] for index in block}
            if len(coordinates)!=1:
                break
            local[next(iter(coordinates))][len(block)-1]+=1
        else:
            key=tuple(tuple(row) for row in local)
            index=len(blocks)
            inner.setdefault(index,{})[key]=inner.setdefault(index,{}).get(key,F(0))+1
    output={}
    for k in range(1,q+1):
        value={}
        for j in range(k,q+1):
            for local,c1 in inner.get(j,{}).items():
                for outer,c2 in ordinary.get((j,k),{}).items():
                    key=(outer,local)
                    value[key]=value.get(key,F(0))+c1*c2
        output[k]=value
    return output


def energy_terms(dimension,radius,p,q):
    if not (type(dimension) is int and dimension>0 and p>=0 and q>=0 and 1<=p+q<=4):
        raise ValueError('Invalid conditional energy')
    ordinary=ordinary_bell()
    gram=conditional_gram(radius)
    result={}
    for alpha in partitions(q):
        length=len(alpha)
        if length>dimension:
            continue
        weight=F(factorial(dimension),factorial(dimension-length))*factorial(q)
        for h in set(alpha):
            weight/=factorial(alpha.count(h))
        for h in alpha:
            weight/=factorial(h)
        bell=boundary_bell(alpha,dimension,ordinary)
        for k,left in bell.items():
            for l,right in bell.items():
                for (outer1,local1),c1 in left.items():
                    for (outer2,local2),c2 in right.items():
                        outer=tuple(x+y for x,y in zip(outer1,outer2))
                        local=tuple(sorted(tuple(x+y for x,y in zip(row1,row2))
                                           for row1,row2 in zip(local1,local2)))
                        for h,c3 in gram[p+k,p+l].items():
                            key=(h,outer,local)
                            result[key]=result.get(key,F(0))+weight*c1*c2*c3
    return {key:value for key,value in result.items() if value}


def serialize_terms(terms):
    return [dict(h=h,outer_exponents=list(outer),local_exponents=[list(row) for row in local],
                 coefficient=str(coefficient))
            for (h,outer,local),coefficient in sorted(terms.items())]
