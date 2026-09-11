"""Exact implications of returned outward enclosures; no floating or sampled tail estimates."""
from fractions import Fraction as Q
import gzip
from hashlib import sha256
import json
from pathlib import Path

from coupled_exact import (candidate, canonical, complete_modes, decode_primitives,
    input_facts, load, omission_cost, primitive, require, selected_modes)
from winding_exact_core import arc_witness, binary_ball, dyad


def validate_mode(record, directory, d, binding):
    require(record["status"] == "MODE_ARITHMETIC_PASS", "Mode arithmetic is incomplete")
    require(record["source_binding"] == binding, "Source binding differs")
    fresh = primitive(d, record["key"])
    require(fresh == decode_primitives(record["primitives"]), "Fresh exact primitive differs")
    degree = record["degree"]
    require(degree == record["settings"]["degree"], "Coefficient degree differs from settings")
    require(len(record["coefficients"]) == degree+1, "Missing formal coefficients")
    intervals = [binary_ball(v) for v in record["coefficients"]]
    width = abs(Q(record["weight"]))*sum(hi-lo for lo,hi in intervals)
    require(str(width) == record["weighted_finite_coefficient_interval_width"], "Coefficient width differs")
    proof = record["geometry"]
    R = Q(proof["radius"])
    require(R > 1 and proof["radius"] in record["settings"]["radii"], "Invalid tail radius")
    L2 = [sum(int(n)**2*abs(Q(c))*R**int(n) for n,c in d[name].items())
          for name in ("P","B","C")]
    require(list(map(str,L2)) == proof["map_second_angular_derivative_bounds"],
            "Map second derivative bound differs")
    with gzip.open(directory/record["geometry_file"], "rt", encoding="utf-8") as stream:
        arcs = [json.loads(line) for line in stream]
    arcs.sort(key=lambda a: Q(a["left"]))
    cursor, maximum = Q(0), Q(0)
    sums = [[Q(0), Q(0)] for _ in range(3)]
    minima = [None]*4
    for arc in arcs:
        left, right = Q(arc["left"]), Q(arc["right"])
        require(left == cursor and left < right <= 2, "Boundary gap or overlap")
        cursor = right
        witness = arc_witness(arc["values"], arc["Z"], arc["endpoint_ratios"])
        require(witness == arc["exact_witness"], "Exact arc implication differs")
        maximum = max(maximum, dyad(witness["modulus_upper"]))
        for j, pair in enumerate(witness["lower"]):
            value = dyad(pair)
            minima[j] = value if minima[j] is None else min(minima[j], value)
        for j, (lo, hi) in enumerate(witness["argument_increments"]):
            sums[j][0] += dyad(lo)
            sums[j][1] += dyad(hi)
    require(cursor == 2 and len(arcs) == proof["leaves"], "Incomplete boundary")
    require(len(arcs) == record["settings"]["panels"] + proof["refinements"],
            "Refinement accounting differs")
    require(all(-3 < lo <= hi < 3 for lo, hi in sums), "Nonzero winding not excluded")
    require(str(maximum) == proof["boundary_modulus_upper"], "Modulus maximum differs")
    require([[str(lo),str(hi)] for lo,hi in sums] == proof["argument_sums"], "Argument sums differ")
    require(list(map(str,minima)) == proof["minima"], "Boundary minimum differs")
    tail = maximum*R**(-2*degree-3)/(1-R**(-2))
    require(str(tail) == proof["scalar_tail_upper"], "Cauchy tail differs")
    return intervals, tail


def coefficient_implication(d, certified, target):
    """Account for EVERY finite mode, including each missing mode's full norm.

    certified maps a canonical key to (weight, coefficient intervals, tail).
    This is a sufficient upper-bound criterion, not an exact evaluation of K_G.
    """
    degree = max([len(v[1])-1 for v in certified.values()]+[0])
    lower, upper = [Q(0)]*(degree+1), [Q(0)]*(degree+1)
    tail = omissions = Q(0)
    visited = set()
    omitted_count = total_count = 0
    for key, weight in complete_modes(d):
        total_count += 1
        if key not in certified:
            omissions += omission_cost(d, key, weight)
            omitted_count += 1
            continue
        expected, intervals, remainder = certified[key]
        require(expected == weight, "Mode coefficient binding differs")
        require(key not in visited and remainder >= 0 and intervals and
                all(lo <= hi for lo,hi in intervals), "Malformed coefficient record")
        visited.add(key)
        for j, (lo, hi) in enumerate(intervals):
            a, b = sorted((weight*lo, weight*hi))
            lower[j] += a
            upper[j] += b
        tail += abs(weight)*remainder
    require(visited == set(certified), "Unrecognized or duplicate finite mode")
    nonlinear = sum(max(abs(a), abs(b)) for a,b in zip(lower[1:],upper[1:]))
    phase = Q(input_facts(d)["phase_error"])
    gamma = lower[0]-nonlinear-tail-omissions-phase
    require(target > 0, "A positive target is required for an upper-bound implication")
    return dict(total_modes=total_count, certified_modes=len(visited), omitted_modes=omitted_count,
        linear_lower=str(lower[0]), nonlinear_upper=str(nonlinear), complete_scalar_tail=str(tail),
        complete_omission_cost=str(omissions), complete_phase_error=str(phase), gamma_lower=str(gamma),
        target=str(target), sufficient_criterion_passed=gamma>target,
        qualification="Conditional on the specified outward-arithmetic execution; no sharpness or exact-constant assertion.")


def validate_selected(output, task, binding):
    output = Path(output)
    d = candidate()
    rows = selected_modes(d, task["labels"])
    validated = []
    for row in rows:
        directory = output/row["label"]
        record = load(directory/"mode.json")
        for name in ("label", "key", "weight", "omission_cost"):
            require(record[name] == row[name], "Selected mode binding differs")
        require(record["settings"] == task["settings"], "Selected mode settings differ")
        intervals, tail = validate_mode(record, directory, d, binding)
        validated.append(dict(label=row["label"], coefficient_count=len(intervals),
            radius=record["geometry"]["radius"], complete_unweighted_tail=str(tail),
            weighted_finite_coefficient_interval_width=record["weighted_finite_coefficient_interval_width"],
            suggested_degree_for_per_mode_tail=record["suggested_degree_for_per_mode_tail"],
            mode_sha256=sha256((directory/"mode.json").read_bytes()).hexdigest(),
            geometry_sha256=sha256((directory/record["geometry_file"]).read_bytes()).hexdigest()))
    result = dict(status=task["stage"].upper()+"_VALIDATED", source_binding=binding,
        selected_modes=validated, full_candidate_certified=False,
        trusted_dependency="Correct outward enclosures from the recorded Python, python-flint, FLINT and Arb execution.")
    (output/"validation.json").write_text(json.dumps(result, indent=2)+"\n")
    return result
