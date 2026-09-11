"""Public input, checkpoint, and exact aggregation utilities."""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import gzip
import json
import sys
from exact_core import digest,inputs,load,modes,omission_cost,require,binary_ball

HERE=Path(__file__).resolve().parents[1]
SETTINGS=dict(bits=768,panels=1024,max_depth=18,max_leaves=100000,mode_seconds=600,worker_memory_bytes=2147483648)
SOURCE_FILES=['run_replay.py','verify_certificate.py','check_algebra.py','sources/package_core.py','sources/native.py','sources/exact_core.py','sources/geometry_check.py','data/candidate.json','data/plan.json','data/primitive_polynomials.json']
BITS=96
AGGREGATION_BITS=256

def unoptimized():
    require(__debug__ and not sys.flags.optimize,'Run Python without -O or PYTHONOPTIMIZE')

def canonical_json(value):return json.dumps(value,sort_keys=True,separators=(',',':'))
def fingerprint(payload):return sha256(canonical_json(payload).encode()).hexdigest()
def source_hashes():return {name:digest(HERE/name) for name in SOURCE_FILES}
def safe_path(root,name):
    require(type(name) is str and name and not Path(name).is_absolute(),'Invalid relative record path')
    p=(root/name).resolve()
    require(p.is_relative_to(root.resolve()),'Record path escapes its directory')
    return p

def check_inputs():
    data=load(HERE/'data/candidate.json');plan=load(HERE/'data/plan.json');formula=load(HERE/'data/primitive_polynomials.json')
    norms,amplitude,phase=inputs(data)
    all_modes=modes(data)
    require(len(all_modes)==2962,'Incorrect complete phase-mode count')
    indices=[r['index'] for r in plan['retained']]
    omitted=plan['omitted_indices']
    require(all(type(i) is int for i in indices+omitted),'Noninteger mode index')
    require(len(indices)==len(set(indices))==2204 and len(omitted)==len(set(omitted))==758,'Missing or duplicate phase modes')
    require(set(indices).isdisjoint(omitted) and set(indices)|set(omitted)==set(range(len(all_modes))),'Incomplete finite Taylor partition')
    require(indices==sorted(indices),'Unordered retained modes')
    for row in plan['retained']:
        key,weight=all_modes[row['index']]
        require(row['key']==list(map(str,key)) and row['weight']==str(weight),'Incorrect exact mode key or weight')
        require(type(row['degree']) is int and 0<=row['degree']<=768 and Q(row['radius'])>1,'Invalid scalar truncation or radius')
    require(sum(r['degree']+1 for r in plan['retained'])==792700,'Incomplete coefficient schedule')
    require(Q(plan['target_gamma'])==Q(71188883,100000000),'Wrong theorem target')
    omission=sum(omission_cost(all_modes[i][0],all_modes[i][1],data) for i in omitted)
    return data,plan,formula,norms,amplitude,phase,omission

def compact_records():
    with gzip.open(HERE/'data/coefficients.jsonl.gz','rt',encoding='utf8') as stream:
        for line in stream:yield json.loads(line)

def compact_interval(value):
    require(isinstance(value,list) and len(value)==2 and all(type(i) is int for i in value),'Invalid compact dyadic interval')
    lo,width=value
    require(width>=0,'Negative interval width')
    return Q(lo,2**BITS),Q(lo+width,2**BITS)

def coefficient_intervals(record,compact):
    convert=compact_interval if compact else binary_ball
    return (convert(c) for c in record['coefficients'])

def geometry_summary(record,row):
    cover=record['geometry']
    require(len(cover['minima'])==4 and all(Q(v)>0 for v in cover['minima']),'Nonpositive geometric summary')
    require(len(cover['argument_sums'])==3 and all(-3<Q(lo)<=Q(hi)<3 for lo,hi in cover['argument_sums']),'Zero winding not established')
    leaves,refinements=cover['leaves'],cover['refinements']
    require(type(leaves) is int and type(refinements) is int and refinements>=0 and leaves==SETTINGS['panels']+refinements and leaves<=SETTINGS['max_leaves'],'Invalid circle subdivision count')
    bound=Q(cover['boundary_modulus_upper'])
    require(bound>0,'Nonpositive circle maximum')
    R=Q(row['radius'])
    tail=bound*R**(-2*row['degree']-3)/(1-R**(-2))
    if 'complete_scalar_tail_upper' in cover:
        require(Q(cover['complete_scalar_tail_upper'])==tail,'Incorrect Cauchy tail')
    return tail,leaves

class Aggregate:
    def __init__(self,plan,phase,omission):
        self.total=[[0,0] for _ in range(1+max(r['degree'] for r in plan['retained']))]
        self.tail=0;self.coefficients=0;self.leaves=0
        self.scale=2**AGGREGATION_BITS;self.phase=phase;self.omission=omission
    def floor(self,q):return q.numerator*self.scale//q.denominator
    def ceil(self,q):return -((-q.numerator*self.scale)//q.denominator)
    def add(self,row,record,compact):
        require(record['index']==row['index'] and len(record['coefficients'])==row['degree']+1,'Missing or reordered coefficient record')
        weight=Q(row['weight'])
        for j,(lo,hi) in enumerate(coefficient_intervals(record,compact)):
            require(lo<=hi,'Reversed coefficient interval')
            lo,hi=sorted((lo*weight,hi*weight))
            self.total[j][0]+=self.floor(lo);self.total[j][1]+=self.ceil(hi)
            self.coefficients+=1
        tail,leaves=geometry_summary(record,row)
        self.tail+=self.ceil(abs(weight)*tail);self.leaves+=leaves
    def result(self,target):
        nonlinear=sum(max(abs(lo),abs(hi)) for lo,hi in self.total[1:])
        gamma=Q(self.total[0][0]-nonlinear-self.tail-self.ceil(self.phase)-self.ceil(self.omission),self.scale)
        require(gamma>Q(target),'Complete coefficient margin does not prove the target')
        return dict(linear_lower=str(Q(self.total[0][0],self.scale)),finite_nonlinear_upper=str(Q(nonlinear,self.scale)),complete_scalar_tail_upper=str(Q(self.tail,self.scale)),complete_phase_tail=str(self.phase),complete_omission_cost=str(self.omission),gamma_lower=str(gamma),target_gamma=target,upper_constant='100000000/71188883',coefficients=self.coefficients,contour_leaves=self.leaves,aggregation_bits=AGGREGATION_BITS)
