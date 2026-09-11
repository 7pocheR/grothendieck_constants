"""Exact implications and coverage of stored interval rectangles; no Arb evaluation."""
from fractions import Fraction as Q
import gzip
import json
from exact_core import arc_witness,dyad,require


def split_box(panel,path,panels):
    lo,hi=Q(2*panel,panels),Q(2*(panel+1),panels)
    for side in path:
        require(type(side) is int and side in (0,1),"Invalid binary subdivision step")
        middle=(lo+hi)/2
        if side==0:
            hi=middle
        else:
            lo=middle
    return lo,hi


def fresh_geometry(path,record,settings):
    cover=record["geometry"]
    left=Q(0)
    count=0
    maximum=Q(0)
    minima=[None]*4
    sums=[[Q(0),Q(0)] for _ in range(3)]
    panels=settings["panels"]
    with gzip.open(path,"rt",encoding="utf8") as stream:
        for line in stream:
            item=json.loads(line)
            panel=item["panel"]
            require(type(panel) is int and 0<=panel<panels,"Invalid initial panel")
            require(len(item["path"])<=settings["max_depth"],"Fresh subdivision exceeds limit")
            expected=split_box(panel,item["path"],panels)
            actual=Q(item["left"]),Q(item["right"])
            require(expected==actual and actual[0]==left and left<actual[1]<=2,
                    "Fresh whole-circle partition has a gap, overlap, or incorrect path")
            require(len(item["values"])==5 and len(item["endpoint_ratios"])==3,
                    "Malformed primitive/endpoint enclosure")
            # Rebuild every rational implication from binary rectangles.
            witness=arc_witness(item["values"],item["Z"],item["endpoint_ratios"])
            require(witness==item["exact_witness"],"Fresh rational arc witnesses differ")
            M=dyad(witness["modulus_upper"])
            maximum=max(maximum,M)
            for j,pair in enumerate(witness["lower"]):
                value=dyad(pair)
                require(value>0,"A fresh leaf is not strictly positive")
                minima[j]=value if minima[j] is None else min(minima[j],value)
            for j,(lo,hi) in enumerate(witness["argument_increments"]):
                sums[j][0]+=dyad(lo);sums[j][1]+=dyad(hi)
            left=actual[1]
            count+=1
    require(left==2 and count==cover["leaves"]==panels+cover["refinements"],
            "Fresh circle coverage is incomplete")
    require(all(-3<lo<=hi<3 for lo,hi in sums),"Fresh zero winding not established")
    require([[str(lo),str(hi)] for lo,hi in sums]==cover["argument_sums"] and
            list(map(str,minima))==cover["minima"],"Fresh geometric summary differs")
    require(maximum==Q(cover["boundary_modulus_upper"]),"Fresh boundary maximum differs")
    R=Q(record["radius"])
    tail=maximum*R**(-2*record["degree"]-3)/(1-R**(-2))
    require(tail==Q(cover["complete_scalar_tail_upper"]),"Fresh complete scalar tail differs")
    return count,tail,minima


