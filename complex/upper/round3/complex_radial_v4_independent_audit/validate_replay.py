"""Independent exact aggregation and geometric coverage validation.

This uses the Python standard library only. It checks every saved leaf,
including the entire rational binary subdivision tree, positive dyadic
witnesses, the squared Euler modulus inequality, and complete Cauchy tails.
It does not reimplement Arb's transcendental interval operations.
"""
import argparse
from fractions import Fraction as Q
import gzip
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
import sys

from audit_exact import abs_ball,ball,decimal,load_json,phase_modes,require,scaled
from rational_core import binary_ball,dyadic,pi_bounds


def check_partition(leaves,root,kind):
    trie={}
    for path,claimed_box in leaves:
        box=list(root)
        node=trie
        for step in path:
            require(isinstance(step,list) and len(step)==2,"Malformed subdivision step")
            axis,side=step
            require(axis in ("r","a") and type(side) is int and side in (0,1),"Invalid subdivision step")
            require(kind!="boundary" or axis=="a","Boundary subdivided radially")
            require("leaf" not in node,"A leaf overlaps its descendant")
            if "axis" in node:
                require(node["axis"]==axis,"Sibling partition axes disagree")
            else:
                node["axis"]=axis
            node=node.setdefault(side,{})
            left=0 if axis=="r" else 2
            mid=(box[left]+box[left+1])/2
            box[left+1 if side==0 else left]=mid
        require(not node,"Duplicate leaf or a leaf overlapping descendants")
        require(tuple(box)==tuple(claimed_box),"Saved coordinates differ from exact split path")
        node["leaf"]=True
    def complete(node):
        if "leaf" in node:
            require(set(node)=={"leaf"},"Leaf has children")
            return
        require(set(node)=={"axis",0,1},"Incomplete partition: an entire child is missing")
        complete(node[0]);complete(node[1])
    complete(trie)


def validate_leaves(path,record,settings,pi_rational_upper):
    cover=record["cover"]
    require(sha256(Path(path).read_bytes()).hexdigest()==record["leaf_sha256"],"Leaf hash mismatch")
    pi_upper=dyadic(cover["pi_upper"])
    require(pi_upper>=pi_rational_upper,"Uncertified upper bound for pi")
    M=dyadic(cover["maximum_boundary_modulus"])
    summary_minima=list(map(dyadic,cover["minimum_witnesses"]))
    require(len(summary_minima)==4 and all(x>0 for x in summary_minima),"Nonpositive summary witness")
    panels,R=settings["panels"],Q(settings["radius"])
    groups=iter((panel,kind) for panel in range(panels) for kind in ("boundary","annulus"))
    current=None
    leaves=[]
    counts={"boundary":0,"annulus":0}
    def finish(group,part):
        if group is None:
            return
        panel,kind=group
        left,right=Q(2*panel,panels),Q(2*(panel+1),panels)
        root=(R,R,left,right) if kind=="boundary" else (Q(1),R,left,right)
        check_partition(part,root,kind)
    with gzip.open(path,"rt",encoding="utf8") as stream:
        for line in stream:
            item=json.loads(line)
            group=item["panel"],item["kind"]
            if group!=current:
                finish(current,leaves)
                require(group==next(groups,None),"Missing, duplicated, or unordered initial panel")
                current=group;leaves=[]
            box=tuple(map(Q,item["box"]))
            require(len(box)==4,"Invalid rectangle")
            leaves.append((item["path"],box))
            require(len(item["path"])<=(16 if group[1]=="boundary" else 28),"Depth limit exceeded")
            ds=list(map(dyadic,item["lower"]))
            require(len(ds)==4 and all(x>0 for x in ds),"A leaf lacks strict positivity")
            require(all(s<=v for s,v in zip(summary_minima,ds)),"Invalid summary minimum")
            if group[1]=="boundary":
                p,m=dyadic(item["numerator_upper"]),dyadic(item["modulus_upper"])
                require(0<=p and 0<=m<=M,"Invalid boundary numerator/modulus")
                D,N1,N2,delta=ds
                require(16*m*m*D*D*N1*N2*delta>=pi_upper*pi_upper*p*p,
                        "Squared Euler boundary inequality failed")
            counts[group[1]]+=1
    finish(current,leaves)
    require(next(groups,None) is None,"Incomplete final groups")
    require(counts["boundary"]==cover["counts"]["boundary"] and counts["annulus"]==cover["counts"]["annulus"],"Leaf count mismatch")
    require(sum(counts.values())==2*panels+cover["counts"]["subdivisions"],"Subdivision count mismatch")
    tail=M*R**(-2*settings["degree"]-3)/(1-R**(-2))
    require(dyadic(cover["scalar_tail_upper"])>=tail,"Recorded Cauchy tail is too small")
    return counts,tail


def validate(args):
    require(args.cluster_resident_handoff and not sys.platform.startswith("darwin"),
            "Full leaf validation belongs on the existing cluster resident")
    directory=Path(args.records)
    payload=load_json(directory/"inputs.json")
    fingerprint=payload.pop("fingerprint")
    require(sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()==fingerprint,"Replay fingerprint mismatch")
    data,formula,settings=payload["data"],payload["formula"],payload["settings"]
    require(data==load_json(args.candidate),"Wrong replay candidate")
    require(formula==load_json(args.formulas),"Wrong determinant-generated formulas")
    audit=Path(__file__).parent
    files={p.name:p for p in [audit/"replay_independent.py",audit/"audit_exact.py",audit/"rational_core.py",Path(args.formulas),Path(args.candidate)]}
    require({name:sha256(path.read_bytes()).hexdigest() for name,path in files.items()}==payload["source_hashes"],"Replay source mismatch")
    require(settings["bits"]>=384 and settings["degree"]>=420 and settings["panels"]>=4096,"Insufficient replay")
    require(settings["omissions"]==[0,165] and Q(settings["radius"])==Q(103,100),"Wrong omission or radius")
    modes=phase_modes(data,data["phase_order"])
    require(len(modes)==210,"Wrong reconstructed mode count")
    expected=set(range(210))-{0,165}
    require({p.name for p in directory.glob("mode_*.json")}=={f"mode_{i:04d}.json" for i in expected},"Missing/additional mode records")
    require({p.name for p in directory.glob("mode_*.leaves.jsonl.gz")}=={f"mode_{i:04d}.leaves.jsonl.gz" for i in expected},"Missing/additional leaf files")
    n=settings["degree"]+1
    total=[[Q(0),Q(0)] for _ in range(n)]
    tail=Q(0)
    rows=[]
    comparisons=0
    pi_upper=pi_bounds()[1]
    for i in sorted(expected):
        key,weight=modes[i]
        record=load_json(directory/f"mode_{i:04d}.json")
        require(record["index"]==i and record["fingerprint"]==fingerprint,"Wrong mode identity")
        require(record["key"]==list(map(str,key)) and record["weight"]==str(weight),"Wrong exact mode")
        require(len(record["coefficients"])==n,"Wrong coefficient count")
        require(record["leaf_file"]==f"mode_{i:04d}.leaves.jsonl.gz","Wrong leaf file name")
        counts,mode_tail=validate_leaves(directory/record["leaf_file"],record,settings,pi_upper)
        for j,value in enumerate(record["coefficients"]):
            lo,hi=scaled(binary_ball(value),weight)
            total[j][0]+=lo;total[j][1]+=hi
        if args.compare_original:
            other=load_json(Path(args.compare_original)/f"mode_{i:04d}.json")
            require(other["key"]==record["key"] and other["weight"]==record["weight"],"Comparison mode differs")
            for j,text in enumerate(other["coefficients"]):
                require(j<n,"Comparison truncation is longer than independent replay")
                first,second=ball(text),binary_ball(record["coefficients"][j])
                require(max(first[0],second[0])<=min(first[1],second[1]),f"Disjoint coefficient intervals at mode {i}, coefficient {j}")
                comparisons+=1
        tail+=abs(weight)*mode_tail
        rows.append(dict(index=i,counts=counts,weighted_scalar_tail_upper=str(abs(weight)*mode_tail)))
        print(json.dumps(dict(validated_mode=i,counts=counts)),flush=True)
    norms={name:sum(abs(Q(v)) for v in data[name].values()) for name in ("P","B")}
    require(norms["P"]<1 and norms["B"]<=1,"Inadmissible maps")
    eta,lam=Q(data["eta"]),abs(Q(data["radial_lam"]))
    require(eta>=0 and data["phase_order"]==4,"Invalid phase")
    radial=Q(1) if eta==0 else min(Q(1),2*lam/(5*eta))
    A=abs(Q(data["epsilon"]))+abs(Q(data["radial_epsilon"]))*radial
    phase=(2*A)**5/factorial(5)
    omissions=sum(abs(modes[i][1])/(1+2*modes[i][0][1])/(1+2*modes[i][0][3]) for i in (0,165))
    nonlinear=sum(abs_ball(x)[1] for x in total[1:])
    gamma=total[0][0]-nonlinear-tail-phase-omissions
    target=Q(data["target_gamma"])
    result=dict(fingerprint=fingerprint,retained_modes=208,coefficients_per_mode=n,
                exact_P_norm=str(norms["P"]),exact_B_norm=str(norms["B"]),
                linear_lower=str(total[0][0]),finite_nonlinear_upper=str(nonlinear),
                complete_scalar_tail_upper=str(tail),complete_phase_tail=str(phase),omission_cost=str(omissions),
                gamma_lower=str(gamma),gamma_lower_decimal=decimal(gamma),target=str(target),
                strict_target_comparison=gamma>target,
                coefficient_intervals_compared=comparisons,
                exact_full_cover_validated=True,modes=rows,
                arithmetic_dependency="Arb supplied each local complex enclosure; the independent validator checks every exact partition and all subsequent inequalities with rational arithmetic")
    with open(args.output,"x") as stream:
        json.dump(result,stream,indent=2);stream.write("\n")
    require(gamma>target,"Full replay does not certify the target")
    print(json.dumps({k:v for k,v in result.items() if k not in ("modes","gamma_lower","complete_scalar_tail_upper","finite_nonlinear_upper","linear_lower")},indent=2),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("records","candidate","formulas","output"):
        p.add_argument("--"+name,required=True)
    p.add_argument("--cluster-resident-handoff",action="store_true")
    p.add_argument("--compare-original")
    validate(p.parse_args())
