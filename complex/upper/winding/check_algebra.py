"""Derive the five Gaussian polynomials and the complete finite phase expansion.

All arithmetic uses the Python standard library. No numerical enclosures are
assumed by these algebraic checks.
"""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import permutations
from pathlib import Path
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent/'sources'))
from exact_core import load,require,modes

ZERO=(0,)*9

def add(*polynomials):
    result=defaultdict(Q)
    for p in polynomials:
        for powers,c in p.items():result[powers]+=c
    return {p:c for p,c in result.items() if c}

def mul(a,b):
    result=defaultdict(Q)
    for p,c in a.items():
        for q,d in b.items():result[tuple(x+y for x,y in zip(p,q))]+=c*d
    return {p:c for p,c in result.items() if c}

def scalar(c):return {ZERO:Q(c)} if c else {}
def neg(a):return {p:-c for p,c in a.items()}
def transpose(a):return list(map(list,zip(*a)))
def matmul(a,b):
    return [[add(*(mul(x,y) for x,y in zip(row,col))) for col in transpose(b)] for row in a]

def det(a):
    total={}
    for p in permutations(range(len(a))):
        term=scalar((-1)**sum(p[i]>p[j] for i in range(len(a)) for j in range(i+1,len(a))))
        for i,j in enumerate(p):term=mul(term,a[i][j])
        total=add(total,term)
    return total

def check_polynomials(formula):
    require(formula['variables']==['a','b','h','r','s','aw','ax','bw','bx'],'Polynomial variable order')
    variables=[{tuple(int(j==i) for j in range(9)):Q(1)} for i in range(9)]
    a,b,h,r,s,aw,ax,bw,bx=variables
    z,o=scalar(0),scalar(1)
    S=[[o,z,a,z],[z,o,z,b],[a,z,o,z],[z,b,z,o]]
    J=[[neg(aw),neg(r),z,z],[r,neg(ax),z,z],[z,z,neg(bw),neg(s)],[z,z,s,neg(bx)]]
    U=[[o,h,z,z],[z,z,o,neg(h)]]
    SJ=matmul(S,J)
    matrix=[[add(scalar(int(i==j)),neg(SJ[i][j])) for j in range(4)] for i in range(4)]
    adj=[]
    for i in range(4):
        row=[]
        for j in range(4):
            minor=[[matrix[k][l] for l in range(4) if l!=i] for k in range(4) if k!=j]
            row.append(mul(scalar((-1)**(i+j)),det(minor)))
        adj.append(row)
    compressed=matmul(matmul(matmul(U,adj),S),transpose(U))
    actual=dict(zip(('D','N1','N2','P0','Q0'),(det(matrix),compressed[0][0],compressed[1][1],compressed[0][1],compressed[1][0])))
    require(set(formula['polynomials'])==set(actual),'Polynomial names')
    for name,p in actual.items():
        terms=formula['polynomials'][name]
        require(len({tuple(t['powers']) for t in terms})==len(terms),'Duplicate polynomial term')
        expected={tuple(t['powers']):Q(t['coefficient']) for t in terms}
        require(p==expected,'Gaussian determinant identity: '+name)
    permutation=[0,1,2,4,3,7,8,5,6]
    for name,target in [('D','D'),('N1','N2'),('N2','N1'),('P0','P0'),('Q0','Q0')]:
        swapped={tuple(p[i] for i in permutation):c for p,c in actual[name].items()}
        require(swapped==actual[target],'Row/column exchange identity')
    # 22/7-pi = integral_0^1 x^4(1-x)^4/(1+x^2) dx > 0.
    dividend={4:Q(1),5:Q(-4),6:Q(6),7:Q(-4),8:Q(1)}
    quotient={}
    while dividend and max(dividend)>=2:
        n=max(dividend);c=dividend[n];quotient[n-2]=c
        for k in (n,n-2):
            dividend[k]=dividend.get(k,0)-c
            if not dividend[k]:del dividend[k]
    require(dividend=={0:Q(-4)} and sum(c/Q(n+1) for n,c in quotient.items())==Q(22,7),'Elementary pi upper bound')
    require(4*sum(Q((-1)**j,2*j+1) for j in range(8))>3,'Elementary pi lower bound')
    return {name:len(p) for name,p in actual.items()}

def expand_modes(data):
    generators=[]
    for side in (0,1):
        for e,l,u in data['profiles']:
            for sign in (-1,1):
                key=[Q(0)]*4
                key[2*side],key[2*side+1]=sign*Q(l),Q(u)
                generators.append((tuple(key),sign*Q(e)/2))
    order={(Q(0),)*4:Q(1)}
    total=defaultdict(Q)
    for n in range(data['phase_order']+1):
        for key,weight in order.items():
            first,second=sorted(((Q(data['k'])-key[0],key[1]),(Q(data['k'])-key[2],key[3])))
            total[first+second]+=weight
        if n<data['phase_order']:
            next_order=defaultdict(Q)
            for key,weight in order.items():
                for generator,c in generators:
                    next_order[tuple(a+b for a,b in zip(key,generator))]+=weight*c/(n+1)
            order={key:c for key,c in next_order.items() if c}
    return sorted((key,c) for key,c in total.items() if c)

def main():
    require(__debug__ and not sys.flags.optimize,'Run Python without -O or PYTHONOPTIMIZE')
    here=Path(__file__).resolve().parent
    counts=check_polynomials(load(here/'data/primitive_polynomials.json'))
    data=load(here/'data/candidate.json')
    independent=expand_modes(data)
    require(independent==modes(data) and len(independent)==2962,'Complete phase expansion')
    require(sum(weight for key,weight in independent)==1,'Merged weight sum')
    print(json.dumps(dict(status='PASS',polynomial_terms=counts,complete_modes=len(independent))))

if __name__=='__main__':main()
