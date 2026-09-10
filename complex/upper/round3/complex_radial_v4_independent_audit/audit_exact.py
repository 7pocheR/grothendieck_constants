"""Independent rational validation of the saved radial certificate records.

This program does not import construction code or trust status fields. It checks
all records, reconstructs the phase weights by exponential multiplication, and
aggregates decimal interval endpoints with fractions.Fraction. Counts do not
replace the missing geometric leaf records or a fresh analytic replay.
"""
import argparse
from collections import defaultdict
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
import re


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_json(path):
    def pairs(items):
        out = {}
        for k, v in items:
            require(k not in out, f"Duplicate JSON key: {k}")
            out[k] = v
        return out
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs)


def phase_modes(data, order):
    """Expand exp(t L) through total order N, using eight elementary terms.

    A key gives (row frequency shift, row damping, column frequency shift,
    column damping). Repeated multiplication and division by n produces L^n/n!.
    The two sine profiles on each side each supply two elementary terms.
    """
    require(type(order) is int and order >= 0, "Invalid Taylor order")
    k, eta = Q(data["k"]), Q(data["eta"])
    require(eta >= 0, "Negative damping")
    elementary = []
    for side in (0, 2):
        for epsilon, frequency, damping in (
            (Q(data["epsilon"]), Q(data["lam"]), Q(0)),
            (Q(data["radial_epsilon"]), Q(data["radial_lam"]), eta),
        ):
            for sign in (-1, 1):
                key = [Q(0)] * 4
                key[side] = sign * frequency
                key[side + 1] = damping
                elementary.append((tuple(key), -sign * epsilon / 2))
    term = {(Q(0),) * 4: Q(1)}
    total = defaultdict(Q, term)
    for n in range(1, order + 1):
        nxt = defaultdict(Q)
        for x, c in term.items():
            for y, d in elementary:
                nxt[tuple(a + b for a, b in zip(x, y))] += c * d / n
        term = {key: c for key, c in nxt.items() if c}
        for key, c in term.items():
            total[key] += c
    symmetric = defaultdict(Q)
    for (dr, alpha, ds, beta), c in total.items():
        left, right = sorted(((k + dr, alpha), (k + ds, beta)))
        symmetric[(*left, *right)] += c
    return sorted((key, c) for key, c in symmetric.items() if c)


NUMBER = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\Z")


def ball(text):
    require(isinstance(text, str), "Interval must be a string")
    s = text.strip()
    if s.startswith("["):
        require(s.endswith("]"), "Malformed interval")
        s = s[1:-1].strip()
    if "+/-" in s:
        mid, rad = (x.strip() for x in s.split("+/-"))
        mid = mid or "0"
        require(bool(NUMBER.fullmatch(mid)) and bool(NUMBER.fullmatch(rad)), "Nonfinite or malformed ball")
        m, r = Q(mid), Q(rad)
        require(r >= 0, "Negative radius")
        return m - r, m + r
    require(bool(NUMBER.fullmatch(s)), f"Malformed exact scalar: {s}")
    return Q(s), Q(s)


def scaled(interval, q):
    x, y = (v * q for v in interval)
    return min(x, y), max(x, y)


def abs_ball(interval):
    lo, hi = interval
    return max(Q(0), lo, -hi), max(abs(lo), abs(hi))


def decimal(q, upward=False, digits=40):
    with localcontext() as ctx:
        ctx.prec = digits
        ctx.rounding = ROUND_CEILING if upward else ROUND_FLOOR
        return str(Decimal(q.numerator) / Decimal(q.denominator))


def interval_json(x):
    return dict(lower=str(x[0]), upper=str(x[1]),
                lower_decimal=decimal(x[0]), upper_decimal=decimal(x[1], True))


def validate(source, record_dir, allow_higher_precision=False):
    source, record_dir = Path(source), Path(record_dir)
    inputs = load_json(record_dir / "inputs.json")
    cert = load_json(record_dir / "certificate.json")
    data = load_json(source / "radial_low_candidate.json")
    require(data == inputs["data"], "Input candidate differs from source copy")
    hashes = {name: sha256((source / name).read_bytes()).hexdigest()
              for name in ("certify_radial.py", "certify_radial_v4.py")}
    require(hashes == inputs["source_hashes"] == cert["source_hashes"], "Source hashes differ")
    settings = inputs["settings"]
    require(settings == cert["settings"], "Certificate settings differ")
    if allow_higher_precision:
        require(settings["bits"] >= 384 and settings["degree"] >= 420
                and settings["panels"] >= 4096 and Q(settings["radius"]) == Q(103,100),
                "Insufficient fresh replay settings")
    else:
        require(settings["bits"] == 256 and settings["degree"] == 420
                and settings["panels"] == 4096 and Q(settings["radius"]) == Q(103, 100),
                "Unexpected original settings")
    require(settings["phase_order"] == data["phase_order"] == 4, "Unexpected order")
    fingerprint = sha256(json.dumps(dict(data=data, settings=settings, source_hashes=hashes),
                                     sort_keys=True).encode()).hexdigest()
    require(fingerprint == inputs["fingerprint"] == cert["fingerprint"], "Fingerprint mismatch")
    norms = {}
    for name in ("P", "B"):
        powers = [int(k) for k in data[name]]
        require(all(str(n) in data[name] and n > 0 and n % 2 == 1 for n in powers), "Not an odd polynomial")
        require(len(set(powers)) == len(powers), "Duplicate polynomial exponent")
        require(all(isinstance(c, str) for c in data[name].values()), "Nonexact polynomial input")
        norms[name] = sum(abs(Q(c)) for c in data[name].values())
        require(norms[name] <= 1, f"Inadmissible {name}")
        require(norms[name] == Q(data[f"exact_{name}_norm"]) == Q(cert[f"exact_{name}_norm"]), "Norm mismatch")
    require(norms["P"] < 1 and norms["B"] == 1, "Unexpected input norms")
    for name in ("theta", "k", "epsilon", "lam", "radial_epsilon", "radial_lam", "eta", "target_gamma"):
        require(isinstance(data[name], str), f"Nonexact parameter: {name}")
        Q(data[name])
    require(Q(data["target_gamma"]) == Q(3559299, 5000000), "Wrong target convention")
    modes = phase_modes(data, data["phase_order"])
    require(len(modes) == 210, "Wrong independent mode count")
    omissions = settings["omitted_mode_indices"]
    require(omissions == [0, 165], "Unexpected omissions")
    costs, omitted_records = Q(0), []
    for i in omissions:
        key, weight = modes[i]
        cost = abs(weight) / (1 + 2 * key[1]) / (1 + 2 * key[3])
        costs += cost
        omitted_records.append(dict(index=i, key=list(map(str, key)), weight=str(weight), exact_norm_cost=str(cost)))
    require(omitted_records == inputs["omitted_modes"] == cert["omitted_modes"], "Omitted modes differ")
    require(costs == Q(data["epsilon"])**4 / 96 == Q(cert["exact_omitted_primitive_norm_cost"]), "Omission cost mismatch")
    expected = set(range(210)) - set(omissions)
    paths = list(record_dir.glob("mode_*.json"))
    require({p.name for p in paths} == {f"mode_{i:04d}.json" for i in expected}, "Missing or additional mode files")
    require(len(paths) == 208 == cert["retained_primitive_count"], "Wrong retained mode count")
    d, panels = settings["degree"], settings["panels"]
    R = Q(settings["radius"])
    factor = R**(-2*d-3) / (1 - R**(-2))
    total = [[Q(0), Q(0)] for _ in range(d+1)]
    saved_tails = [Q(0), Q(0)]
    reconstructed_tail = Q(0)
    rows = []
    wide_minima = defaultdict(list)
    raw_hashes = {}
    for i in sorted(expected):
        path = record_dir / f"mode_{i:04d}.json"
        mode = load_json(path)
        raw_hashes[path.name] = sha256(path.read_bytes()).hexdigest()
        key, weight = modes[i]
        require(mode["index"] == i and mode["fingerprint"] == fingerprint, f"Wrong identity {i}")
        require(mode["key"] == list(map(str, key)) and mode["weight"] == str(weight), f"Wrong exact mode {i}")
        require(key[1] >= 0 and key[3] >= 0, f"Negative mode damping {i}")
        require(len(mode["coefficients"]) == d + 1, f"Coefficient count {i}")
        for j, text in enumerate(mode["coefficients"]):
            lo, hi = scaled(ball(text), weight)
            total[j][0] += lo
            total[j][1] += hi
        tail = ball(mode["unweighted_scalar_tail"])
        require(tail[1] > 0, f"Invalid scalar tail {i}")
        for side in (0, 1):
            saved_tails[side] += abs(weight) * tail[side]
        proof = mode["proof"]
        M = ball(proof["boundary_modulus_bound"])
        require(M[1] > 0, f"Invalid boundary upper bound {i}")
        reconstructed_tail += abs(weight) * M[1] * factor
        # Necessary bookkeeping condition only, not an independent leaf cover.
        b, a, n = (proof[k] for k in ("final_boundary_arcs", "final_annulus_rectangles", "total_subdivisions"))
        require(all(type(x) is int for x in (a, b, n)), f"Nonintegral counts {i}")
        require(proof["initial_boundary_arcs"] == panels and a >= panels and b >= panels
                and a + b == 2*panels + n, f"Inconsistent subdivision count {i}")
        require(b <= panels*2**12 and a <= panels*2**24, f"Exceeds depth bounds {i}")
        for field in ("annulus_minimum_D_modulus", "annulus_minimum_N1_modulus", "annulus_minimum_N2_modulus",
                      "minimum_Euler_distance_boundary", "minimum_Euler_distance_annulus"):
            interval = ball(proof[field])
            require(interval[1] > 0, f"Nonpositive upper endpoint in {field}, mode {i}")
            if interval[0] <= 0:
                wide_minima[field].append(i)
        rows.append(dict(index=i, key=list(map(str, key)), weight=str(weight), coefficients_checked=d+1,
                         boundary_arcs=b, annulus_rectangles=a, subdivisions=n,
                         boundary_modulus_upper=str(M[1]), saved_tail_upper=str(tail[1])))
    nonlinear = [Q(0), Q(0)]
    for x in total[1:]:
        lo, hi = abs_ball(x)
        nonlinear[0] += lo
        nonlinear[1] += hi
    eta, lam = Q(data["eta"]), abs(Q(data["radial_lam"]))
    radial = Q(1) if eta == 0 else min(Q(1), 2*lam/(5*eta))
    phase_bound = abs(Q(data["epsilon"])) + abs(Q(data["radial_epsilon"]))*radial
    phase_tail = (2*phase_bound)**5/factorial(5)
    require(phase_bound == Q(cert["exact_phase_multiplier_bound"]), "Phase bound mismatch")
    require(phase_tail == Q(cert["complete_phase_Taylor_tail"]), "Phase remainder mismatch")
    gamma_saved = total[0][0] - nonlinear[1] - saved_tails[1] - phase_tail - costs
    gamma_reconstructed = total[0][0] - nonlinear[1] - reconstructed_tail - phase_tail - costs
    target = Q(data["target_gamma"])
    require(gamma_saved > target and gamma_reconstructed > target, "Independent rational aggregate misses target")
    comparisons = {}
    for name, ours in (("linear_coefficient", total[0]), ("finite_nonlinear_norm", nonlinear),
                       ("complete_scalar_degree_tail", saved_tails)):
        theirs = ball(cert[name])
        require(theirs[0] <= ours[0] and ours[1] <= theirs[1], f"Saved aggregate does not enclose exact decimal aggregation: {name}")
        comparisons[name] = interval_json(ours)
    return dict(
        conclusion="All 208 saved mode records pass independent rational input, weight, schema, and aggregation checks. Geometric leaves and fresh replay remain unverified.",
        input_fingerprint=fingerprint, source_hashes=hashes, mode_file_hashes=raw_hashes,
        exact_P_norm=str(norms["P"]), exact_B_norm=str(norms["B"]), modes=rows,
        aggregate=comparisons, omission_cost=str(costs), phase_bound=str(phase_bound), phase_tail=str(phase_tail),
        reconstructed_scalar_tail_upper=str(reconstructed_tail), reconstructed_scalar_tail_upper_decimal=decimal(reconstructed_tail, True),
        gamma_from_saved_tails_lower=str(gamma_saved), gamma_from_saved_tails_lower_decimal=decimal(gamma_saved),
        gamma_from_boundary_bounds_lower=str(gamma_reconstructed), gamma_from_boundary_bounds_lower_decimal=decimal(gamma_reconstructed),
        target=str(target), margin_lower_decimal=decimal(gamma_reconstructed-target),
        conservative_gamma_rational="71185999/100000000",
        conservative_gamma_verified=gamma_reconstructed > Q(71185999, 100000000),
        serialized_minima_containing_zero=dict(wide_minima),
        coverage_note="Only source traversal and the necessary leaf/subdivision count identity can be checked from these records. Exact leaf rectangles and leaf arithmetic were not saved.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--records", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--allow-higher-precision", action="store_true")
    args = parser.parse_args()
    result = validate(args.source, args.records, args.allow_higher_precision)
    with open(args.output, "x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("conclusion", "gamma_from_boundary_bounds_lower_decimal", "margin_lower_decimal", "conservative_gamma_verified")}, indent=2))
