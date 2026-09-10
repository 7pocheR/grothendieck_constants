"""Bounded adversarial tests of exact serialization and coverage validation."""
import ast
from fractions import Fraction as Q
import json
from pathlib import Path

from audit_exact import ball,require
from rational_core import binary_ball,pi_bounds,sqrt_bounds
from validate_replay import check_partition


def tests():
    folder=Path(__file__).parent
    for path in folder.glob("*.py"):
        ast.parse(path.read_text(),filename=str(path))
    require(ball("[+/- 1.16]")==(-Q(29,25),Q(29,25)),"Zero-centred ball parser")
    require(ball("[-1.23e-4 +/- 2e-7]")==(-Q(616,5000000),-Q(614,5000000)),"Signed scientific ball parser")
    require(binary_ball({"mid":[-3,-2],"rad":[1,-4]})==(-Q(13,16),-Q(11,16)),"Exact dyadic decoding")
    root=(Q(1),Q(103,100),Q(0),Q(1,4))
    valid=[
        ([["r",0]],(Q(1),Q(203,200),Q(0),Q(1,4))),
        ([["r",1],["a",0]],(Q(203,200),Q(103,100),Q(0),Q(1,8))),
        ([["r",1],["a",1]],(Q(203,200),Q(103,100),Q(1,8),Q(1,4))),
    ]
    check_partition(valid,root,"annulus")
    check_partition([([],root)],root,"annulus")
    failures={
        "missing_child":valid[:-1],
        "duplicate_leaf":valid+[valid[0]],
        "parent_overlaps_descendant":valid+[([],root)],
        "wrong_coordinate":[(valid[0][0],root)]+valid[1:],
        "mismatched_split_axes":[valid[0],([["a",1]],(Q(1),Q(103,100),Q(1,8),Q(1,4)))],
    }
    caught=[]
    for name,leaves in failures.items():
        try:
            check_partition(leaves,root,"annulus")
        except ValueError:
            caught.append(name)
        else:
            raise ValueError("Invalid cover accepted: "+name)
    pi=pi_bounds()
    require(Q("3.1415926535897932384626433832795028841971")<pi[0]<pi[1]<Q("3.1415926535897932384626433832795028841972"),"Machin interval sanity")
    lo,hi=sqrt_bounds(Q(17,23),digits=80)
    require(lo*lo<=Q(17,23)<=hi*hi,"Square root rational isolation")
    return dict(python_files_parse=True,exact_interval_parsing=True,valid_partitions_accepted=True,
                invalid_partitions_rejected=caught,pi_and_sqrt_rational_bounds=True,
                limitation="These tests do not execute the full Arb replay")


if __name__=="__main__":
    print(json.dumps(tests(),indent=2))
