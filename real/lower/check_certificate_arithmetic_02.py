"""Exact coverage and arithmetic of the extracted mathematical data.

This checks every index and sum, but does not evaluate any vertex integral,
scalar-cell integral, or quadrature panel.
"""
from fractions import Fraction as Q
from pathlib import Path
import itertools
import json
import sys

if not __debug__:
    raise RuntimeError("Assertions must be enabled")
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent


def box(path):
    result = [(Q(0),Q(2)), (Q(0),Q(2)), (Q(-2),Q(2)), (Q(-2),Q(2))]
    for bit in path:
        assert bit in "01"
        widths = [u-l for l,u in result]
        axis = widths.index(max(widths))
        lo,hi = result[axis]
        mid = (lo+hi)/2
        result[axis] = (lo,mid) if bit == "0" else (mid,hi)
    return result


def main():
    fiber = json.loads((HERE/"certificates_02/fibers.json").read_text())
    assert fiber["domain"] == [["0","2"],["0","2"],["-2","2"],["-2","2"]]
    assert fiber["odd_matrix"] == dict(p="139/118",q="63/59",r_squared="1323/6962")
    gaps, counts = {}, {}
    for name, data in fiber["targets"].items():
        p,q,r2,t,kappa = (Q(data[k]) for k in ("p","q","r_squared","t","kappa"))
        assert (p,t,kappa) == ((Q(21,20),Q(3,2),Q(26,125)) if name == "A" else (Q(13,10),Q(11,5),Q(1,6000)))
        assert q == (p+Q(13,8))/Q(21,8) and r2 == (p-1)*q
        paths = data["leaf_paths"]
        assert sorted(set(paths)) == paths
        assert all(not y.startswith(x) for x,y in zip(paths,paths[1:]))
        assert sum((Q(1,2**len(path)) for path in paths),Q(0)) == 1
        vertices = {tuple(map(Q,row[:4])):Q(row[4]) for row in data["vertices"]}
        assert len(vertices) == len(data["vertices"]) and all(len(row)==5 for row in data["vertices"])
        needed, minimum = set(), None
        for path in paths:
            intervals = box(path)
            corners = list(itertools.product(*intervals))
            needed.update(corners)
            correction = sum(((u-l)**2/4 for l,u in intervals),Q(0))
            gap = t+kappa-correction-max(vertices[v] for v in corners)
            assert gap > 0
            minimum = gap if minimum is None else min(minimum,gap)
        assert set(vertices) == needed
        gaps[name],counts[name] = str(minimum),dict(leaves=len(paths),vertices=len(vertices))
    assert gaps == {"A":"20720380967041001/144115188075855872000", "B":"311040589043071/27021597764222976000"}
    supports = []
    for i, simple_gap in enumerate((Q(1,5000),Q(7,10000),Q(3,10000)),1):
        data = json.loads((HERE/f"certificates_02/scalar_{i}.json").read_text())
        d = data["parameters"]
        cells = [(Q(l),Q(u),Q(value)) for l,u,value in data["cells"]]
        assert len(cells) == (2413 if i == 1 else 60001)
        assert cells[0][0] == 0 and cells[-1][1] == 12
        assert all(l<u and value>=0 for l,u,value in cells)
        assert all(a[1]==b[0] for a,b in zip(cells,cells[1:]))
        assert 2*sum((value for l,u,value in cells),Q(0))+Q(data["tail_upper"]) == Q(data["integral_upper"])
        panels = [tuple(map(Q,row)) for row in data["noise_panel_bounds"]]
        assert len(panels)==1024 and all(0<=l<=u for l,u in panels)
        error = Q(55*50**8,44100)*Q(3,1024)**8
        tangent = data["tangent_bounds"]
        assert Q(tangent["quadrature_error"]) == error
        assert sum((l for l,u in panels),Q(0))-error == Q(tangent["F_lower"])
        assert sum((u for l,u in panels),Q(0))+error == Q(tangent["F_upper"])
        assert 0 <= Q(tangent["intercept_lower"]) <= Q(tangent["intercept_upper"]) <= Q(d["d0"])
        assert 0 <= Q(tangent["slope_lower"]) <= Q(tangent["slope_upper"]) <= Q(d["d1"])
        assert Q(d["r"])**2 < Q(tangent["turning_radius_squared_lower"])
        b = Q(d["a"])-Q(d["alpha"])/Q(d["rho"])**2
        assert b == Q(data["quadratic"]) and b>0
        constant = Q(data["coarse_constant_upper"])
        assert constant >= Q(data["constant_upper"])
        endpoints = [b*Q(x)**2 + Q(d["ell"])*Q(x) + constant for x in data["x_interval"]]
        assert endpoints == list(map(Q,data["coarse_endpoint_upper"]))
        gap = Q(2917,2000)-max(endpoints)
        assert gap == Q(data["coarse_minimum_gap"]) and gap > simple_gap
        supports.append(dict(case=i,cells=len(cells),noise_panels=len(panels),gap=str(gap),simple_gap_lower=str(simple_gap)))
    returned = json.loads((HERE/"final_arithmetic_02/scalar.json").read_text())
    joint_coarse = Q(58306173,40000000)
    assert Q(returned["final"]["joint_upper"]) <= joint_coarse < Q(2917,2000)
    result = dict(fiber_counts=counts,fiber_minimum_gaps=gaps,scalar_supports=supports,
                  joint_upper=str(joint_coarse),joint_gap=str(Q(2917,2000)-joint_coarse),
                  integrals_evaluated=0)
    (HERE/"certificate_arithmetic_02.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
