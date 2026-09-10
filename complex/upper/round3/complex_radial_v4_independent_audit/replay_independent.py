"""Full interval replay intended only for the existing cluster resident.

Primitive coefficients come from determinant-generated rational polynomials.
Taylor weights are reconstructed by repeated multiplication of eight terms.
Every boundary/annulus leaf is saved, with exact rational coordinates and
dyadic witnesses. No checkpoint is reused. This is distinct from the supplied
construction implementation, but still uses the Arb arithmetic library.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Q
import gzip
from hashlib import sha256
import importlib.metadata
import json
from math import comb,factorial
from pathlib import Path
import platform
import sys
import time
import traceback

from flint import arb,acb,arb_series,acb_poly,ctx
from audit_exact import load_json,phase_modes,require
from rational_core import specialize


def rat(q):
    q=Q(q)
    return arb(q.numerator)/q.denominator


def endpoint(x,upper=False,bits=48):
    """Round an outward endpoint again, in the same direction, to a dyadic."""
    y=x.upper() if upper else x.lower()
    m,e=map(int,y.man_exp())
    shift=max(0,abs(m).bit_length()-bits)
    if shift:
        m=-((-m)//(1<<shift)) if upper else m//(1<<shift)
        e+=shift
    return [m,e]


def exact(pair):
    return arb(tuple(pair))


def pack(x):
    require(x.is_finite(),"Nonfinite coefficient ball")
    return {"mid":list(map(int,x.mid().man_exp())),"rad":list(map(int,x.rad().man_exp()))}


def formal_coefficients(data,primitives,degree):
    n=degree+1
    arrays=[]
    for name in ("P","B"):
        values=[arb(0)]*n
        for exponent,c in data[name].items():
            j=(int(exponent)-1)//2
            if j<n:
                values[j]=rat(c)
        arrays.append(arb_series(values))
    A,B=arrays
    y=arb_series([0,1])
    values={}
    for name,poly in primitives.items():
        parity=int(name in ("P0","Q0"))
        value=arb_series([0])
        for (i,j),c in poly.items():
            require((i+j)%2==parity,"Primitive has incorrect parity")
            value+=rat(c)*y**((i+j-parity)//2)*A**i*B**j
        values[name]=value
    D,N1,N2,L,Q0=(values[name] for name in ("D","N1","N2","P0","Q0"))
    T=N1*N2
    T0=primitives["N1"][(0,0)]*primitives["N2"][(0,0)]
    require(T0>0 and D[0]>0,"Invalid origin branch")
    Z=y*L*Q0/T
    require(Z[0]==0,"Nonzero hypergeometric constant")
    HG=arb_series([1])
    power=arb_series([1])
    for j in range(1,n):
        power*=Z
        HG+=rat(Q(comb(2*j,j)**2,16**j*(j+1)))*power
    value=arb.pi()/(4*rat(T0).sqrt())*L/D*HG/(T/rat(T0)).sqrt()
    return [pack(value[j]) for j in range(n)]


def compiled_primitives(primitives):
    """Nested Horner evaluation of exact, already-collected a,b coefficients."""
    out={}
    for name,poly in primitives.items():
        imax=max(i for i,j in poly)
        jmax=max(j for i,j in poly)
        out[name]=[[rat(poly.get((i,j),0)) for j in range(jmax+1)] for i in range(imax+1)]
    return out


def evaluate(compiled,a,b):
    out=[]
    for name in ("D","N1","N2","P0","Q0"):
        value=acb(0)
        for row in reversed(compiled[name]):
            sub=acb(0)
            for c in reversed(row):
                sub=sub*b+c
            value=value*a+sub
        out.append(value)
    return out


def box_ball(lo,hi):
    require(lo<=hi,"Reversed interval")
    return rat((lo+hi)/2)+arb(0,rat((hi-lo)/2))


def check_box(box,compiled,P,B,pi_upper,boundary):
    rlo,rhi,alo,ahi=box
    z=box_ball(rlo,rhi)*acb(0,arb.pi()*box_ball(alo,ahi)).exp()
    a,b=z*P(z*z),z*B(z*z)
    D,N1,N2,p,q=evaluate(compiled,a,b)
    ds=[endpoint(abs(v)) for v in (D,N1,N2)]
    if not all(exact(v)>0 for v in ds):
        return None
    Z=p*q/(N1*N2)
    real=min(arb(1),(1-Z.real).lower())
    modulus=abs(Z).upper()
    imaginary=arb(0) if not modulus>0 else (abs(Z.imag).lower()/modulus).lower()
    delta=max(real,imaginary)
    delta_pair=endpoint(delta)
    if not exact(delta_pair)>0:
        return None
    witnesses=dict(lower=ds+[delta_pair])
    if boundary:
        p_upper=endpoint(abs(p),True)
        denominator=4*exact(ds[0])*(exact(ds[1])*exact(ds[2])*exact(delta_pair)).sqrt()
        M=endpoint(exact(pi_upper)*exact(p_upper)/denominator,True)
        witnesses.update(numerator_upper=p_upper,modulus_upper=M)
    return witnesses


def children(box,axis):
    rlo,rhi,alo,ahi=box
    if axis=="r":
        mid=(rlo+rhi)/2
        return (rlo,mid,alo,ahi),(mid,rhi,alo,ahi)
    mid=(alo+ahi)/2
    return (rlo,rhi,alo,mid),(rlo,rhi,mid,ahi)


def full_cover(data,primitives,settings,stream):
    R=Q(settings["radius"])
    panels=settings["panels"]
    pi_upper=endpoint(arb.pi(),True)
    polynomials=[]
    for name in ("P","B"):
        coefficients=[arb(0)]*((max(map(int,data[name]))+1)//2)
        for n,c in data[name].items():
            coefficients[(int(n)-1)//2]=rat(c)
        polynomials.append(acb_poly(coefficients))
    P,B=polynomials
    compiled=compiled_primitives(primitives)
    maximum=arb(0)
    minima=[None]*4
    counts={"boundary":0,"annulus":0,"subdivisions":0}
    for panel in range(panels):
        left,right=Q(2*panel,panels),Q(2*(panel+1),panels)
        for kind,box in (("boundary",(R,R,left,right)),("annulus",(Q(1),R,left,right))):
            stack=[(box,[])]
            while stack:
                box,path=stack.pop()
                witnesses=check_box(box,compiled,P,B,pi_upper,kind=="boundary")
                if witnesses is not None:
                    record=dict(kind=kind,panel=panel,path=path,box=list(map(str,box)),**witnesses)
                    stream.write(json.dumps(record,separators=(",",":"))+"\n")
                    counts[kind]+=1
                    require(counts["boundary"]+counts["annulus"]<=settings["max_leaves"],"Per-mode leaf bound exhausted")
                    for i,v in enumerate(witnesses["lower"]):
                        value=exact(v)
                        minima[i]=value if minima[i] is None else min(minima[i],value)
                    if kind=="boundary":
                        maximum=max(maximum,exact(witnesses["modulus_upper"]))
                    continue
                limit=16 if kind=="boundary" else 28
                require(len(path)<limit,f"Subdivision depth exhausted: {kind}, panel {panel}, box {box}")
                rlo,rhi,alo,ahi=box
                axis="r" if kind=="annulus" and rhi-rlo>=3*R*(ahi-alo) else "a"
                low,high=children(box,axis)
                stack.append((high,path+[[axis,1]]))
                stack.append((low,path+[[axis,0]]))
                counts["subdivisions"]+=1
    tail=maximum*rat(R)**(-2*settings["degree"]-3)/(1-rat(R)**(-2))
    return dict(pi_upper=pi_upper,maximum_boundary_modulus=endpoint(maximum,True),
                scalar_tail_upper=endpoint(tail,True),minimum_witnesses=[endpoint(v) for v in minima],
                counts=counts)


def one_mode(task):
    index,key,weight,data,formula,settings,destination,fingerprint=task
    ctx.prec=settings["bits"]
    ctx.cap=settings["degree"]+1
    destination=Path(destination)
    started=time.time()
    leaf_path=destination/f"mode_{index:04d}.leaves.jsonl.gz"
    try:
        primitives=specialize(formula,Q(data["theta"]),key)
        coefficients=formal_coefficients(data,primitives,settings["degree"])
        with gzip.open(leaf_path,"xt",encoding="utf8",compresslevel=3) as stream:
            cover=full_cover(data,primitives,settings,stream)
        record=dict(index=index,key=list(map(str,key)),weight=str(weight),fingerprint=fingerprint,
                    coefficients=coefficients,cover=cover,
                    leaf_file=leaf_path.name,leaf_sha256=sha256(leaf_path.read_bytes()).hexdigest(),
                    elapsed_seconds=time.time()-started)
        with (destination/f"mode_{index:04d}.json").open("x") as stream:
            json.dump(record,stream,indent=2);stream.write("\n")
        print(json.dumps(dict(mode=index,completed=True,elapsed_seconds=record["elapsed_seconds"])),flush=True)
        return index
    except BaseException:
        with (destination/f"mode_{index:04d}.failure.txt").open("x") as stream:
            traceback.print_exc(file=stream)
        raise


def run(args):
    require(args.cluster_resident_handoff,"Full replay is reserved for the existing cluster resident")
    require(not sys.platform.startswith("darwin"),"Refusing full computation on the local Mac")
    require(args.bits>=384 and args.degree>=420 and args.panels>=4096,"Replay precision/settings are insufficient")
    require(1<=args.workers<=4,"At most four resident CPUs")
    destination=Path(args.output)
    destination.mkdir(parents=True,exist_ok=False)
    data,formula=load_json(args.candidate),load_json(args.formulas)
    for key in ("P","B"):
        require(sum(abs(Q(v)) for v in data[key].values())<=1,"Inadmissible correlation map")
        require(all(int(n)>0 and int(n)%2 for n in data[key]),"Nonodd correlation map")
    require(Q(data["eta"])>=0 and data["phase_order"]==4,"Unsupported phase input")
    require(Q(data["target_gamma"])==Q(3559299,5000000),"Unexpected target")
    settings=dict(bits=args.bits,degree=args.degree,panels=args.panels,radius=args.radius,
                  max_leaves=args.max_leaves,omissions=[0,165],workers=args.workers)
    require(Q(args.radius)>1,"Cauchy radius must exceed one")
    files=[Path(__file__),Path(__file__).with_name("audit_exact.py"),Path(__file__).with_name("rational_core.py"),Path(args.formulas),Path(args.candidate)]
    hashes={p.name:sha256(p.read_bytes()).hexdigest() for p in files}
    environment=dict(python=sys.version,platform=platform.platform(),python_flint=importlib.metadata.version("python-flint"))
    payload=dict(data=data,formula=formula,settings=settings,source_hashes=hashes,environment=environment)
    fingerprint=sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()
    with (destination/"inputs.json").open("x") as f:
        json.dump(dict(**payload,fingerprint=fingerprint),f,indent=2);f.write("\n")
    modes=phase_modes(data,data["phase_order"])
    require(len(modes)==210,"Unexpected phase mode count")
    tasks=[(i,key,w,data,formula,settings,str(destination),fingerprint)
           for i,(key,w) in enumerate(modes) if i not in settings["omissions"]]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        indices=list(pool.map(one_mode,tasks,chunksize=1))
    with (destination/"completion.json").open("x") as f:
        json.dump(dict(completed_mode_indices=indices,fingerprint=fingerprint,
                       next_step="Run validate_replay.py independently; this file is not an acceptance certificate"),f,indent=2)
        f.write("\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    for key in ("candidate","formulas","output"):
        parser.add_argument("--"+key,required=True)
    parser.add_argument("--cluster-resident-handoff",action="store_true")
    parser.add_argument("--workers",type=int,default=4)
    parser.add_argument("--bits",type=int,default=512)
    parser.add_argument("--degree",type=int,default=480)
    parser.add_argument("--panels",type=int,default=8192)
    parser.add_argument("--radius",default="103/100")
    parser.add_argument("--max-leaves",type=int,default=200000)
    run(parser.parse_args())
