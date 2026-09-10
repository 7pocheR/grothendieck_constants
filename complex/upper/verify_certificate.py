#!/usr/bin/env python3
"""Check the supplied complex Grothendieck certificate with exact arithmetic.

This entry point checks every selected file hash, all 100048 coefficient
records, and the complete rational aggregation. It does not reevaluate
the local Arb enclosures or traverse the geometric subdivision trees.
See README.md for the separate full validation and recomputation commands.
"""
import argparse
from collections import defaultdict
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from hashlib import sha256
import json
from math import factorial, isqrt
from pathlib import Path, PurePosixPath
import re
import sys
import time


AUDIT = Path("round3/complex_radial_v4_independent_audit")
RECORDS = Path("cluster/resident_20260907/results/complex_upper_radial_independent/cluster_attempt/independent_512_480_8192")
PINNED = {
    "proof_final.md": "aa2c85749501e50155dcbd54c1eda7401bcaf697afa4314e396ea5b611271bc3",
    "certificate_final.json": "00966c5de73e6d1152ff9b866d00ac48c6c675e128f122e3ab3c987b2a1d3b57",
    "primitive_polynomials.json": "df0c52e37de5eaa066b088ec5c7ca8f12755cf988ff04ce78553a637b3c5804c",
    "source_snapshot/radial_low_candidate.json": "6f944ffcbeebca2c850fdf3bbe8e6bcb9ca47d9962b4c8da1df18f0486698ce3",
}
REQUIRED_SOURCE_FILES = (
    "proof_final.md", "certificate_final.json", "primitive_polynomials.json",
    "audit_exact.py", "rational_core.py", "replay_independent.py",
    "validate_replay.py", "check_identities.py", "check_validator.py",
    "check_rational_tails.py", "check_rational_tails_v2.py",
    "source_snapshot/certify_radial.py", "source_snapshot/radial_low_candidate.json",
)


def need(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h = sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def unique_pairs(pairs):
    value = {}
    for key, item in pairs:
        need(key not in value, "Duplicate JSON object key: " + key)
        value[key] = item
    return value


def load(path):
    def bad_constant(value):
        raise ValueError("Nonfinite JSON number: " + value)
    return json.loads(Path(path).read_text(encoding="utf8"),
                      object_pairs_hook=unique_pairs, parse_constant=bad_constant)


def no_floats(value):
    need(type(value) is not float, "Floating-point value in mathematical certificate")
    if isinstance(value, dict):
        for item in value.values():
            no_floats(item)
    elif isinstance(value, list):
        for item in value:
            no_floats(item)


def dyad(value):
    need(type(value) is list and len(value) == 2 and
         all(type(x) is int for x in value), "Malformed exact dyadic")
    m, e = value
    return F(m * 2**e) if e >= 0 else F(m, 2**(-e))


def interval(value):
    need(type(value) is dict and set(value) == {"mid", "rad"}, "Malformed coefficient ball")
    m, r = dyad(value["mid"]), dyad(value["rad"])
    need(r >= 0, "Negative coefficient radius")
    return m-r, m+r


def decimal_bound(value, upper=False):
    with localcontext() as context:
        context.prec = 40
        context.rounding = ROUND_CEILING if upper else ROUND_FLOOR
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def verify(root, manifest_path):
    started = time.monotonic()
    root = Path(root).resolve()
    manifest = load(manifest_path)
    need(manifest["schema_version"] == 1, "Unsupported integrity manifest")
    files = manifest["files"]
    need(type(files) is list and files, "Empty integrity manifest")
    checked = {}
    total_bytes = 0
    for entry in files:
        name = entry["path"]
        need(type(name) is str and "\\" not in name, "Invalid relative file path")
        relative = PurePosixPath(name)
        need(not relative.is_absolute() and relative.parts and
             all(p not in (".", "..") for p in relative.parts), "Unsafe relative file path")
        need(name == str(relative) and name not in checked, "Duplicate or noncanonical file path")
        path = root / name
        need(path.is_file() and not path.is_symlink(), "Missing or symbolic-link file: " + name)
        need(path.resolve().is_relative_to(root), "File resolves outside package")
        need(type(entry["bytes"]) is int and entry["bytes"] >= 0, "Invalid recorded file size")
        need(re.fullmatch(r"[0-9a-f]{64}", entry["sha256"]) is not None, "Invalid SHA256")
        need(path.stat().st_size == entry["bytes"], "File size changed: " + name)
        actual = digest(path)
        need(actual == entry["sha256"], "File hash changed: " + name)
        checked[name] = actual
        total_bytes += entry["bytes"]

    audit, records = root / AUDIT, root / RECORDS
    expected_indices = set(range(210)) - {0,165}
    required_paths = {str(AUDIT / p) for p in REQUIRED_SOURCE_FILES}
    required_paths.add(str(RECORDS / "inputs.json"))
    for i in expected_indices:
        required_paths.add(str(RECORDS / f"mode_{i:04d}.json"))
        required_paths.add(str(RECORDS / f"mode_{i:04d}.leaves.jsonl.gz"))
    need(required_paths <= set(checked), "Integrity manifest omits a required source or data file")
    for relative, expected in PINNED.items():
        need(checked.get(str(AUDIT / relative)) == expected, "A fixed mathematical input changed: " + relative)
    certificate = load(audit / "certificate_final.json")
    candidate = load(audit / "source_snapshot/radial_low_candidate.json")
    primitives = load(audit / "primitive_polynomials.json")
    payload = load(records / "inputs.json")
    fingerprint = payload.pop("fingerprint")
    need(sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest() == fingerprint,
         "Computation input fingerprint differs")
    need(payload["data"] == candidate and payload["formula"] == primitives, "Computation inputs differ")
    settings = payload["settings"]
    need(settings == dict(bits=512, degree=480, panels=8192, radius="103/100",
                          max_leaves=200000, omissions=[0,165], workers=4), "Computation settings differ")
    source_paths = [audit / x for x in ("replay_independent.py", "audit_exact.py", "rational_core.py",
                                       "primitive_polynomials.json", "source_snapshot/radial_low_candidate.json")]
    need({p.name: checked.get(str(p.relative_to(root))) for p in source_paths} == payload["source_hashes"],
         "Frozen computation source hashes differ")
    # Import only after its bytes have been checked against the fixed input.
    sys.path.insert(0, str(audit))
    from audit_exact import phase_modes

    no_floats(certificate)
    keys = {"theta", "k", "epsilon", "lam", "radial_epsilon", "radial_lam", "eta", "P", "B"}
    need(set(certificate["parameters"]) == keys and
         certificate["parameters"] == {k: candidate[k] for k in keys}, "Certificate parameters differ")
    need(certificate["primitive_definitions"] == "primitive_polynomials.json" and
         primitives["variables"] == ["a", "b", "h", "r", "s", "alpha", "beta"], "Primitive definitions differ")
    need(certificate["phase_order"] == candidate["phase_order"] == 4 and
         certificate["scalar_degree_parameter"] == 480 and
         certificate["initial_angular_panels"] == 8192 and
         F(certificate["cauchy_radius"]) == F(103,100), "Mathematical truncation differs")
    norms = {}
    for name in ("P", "B"):
        need(candidate[name] and all(str(int(n)) == n and int(n) > 0 and int(n) % 2
                                     for n in candidate[name]), "Invalid correlation exponent")
        norms[name] = sum(abs(F(c)) for c in candidate[name].values())
        need(norms[name] <= 1, "Inadmissible tensor polynomial")
    need(norms["P"] == F("0.999999680074652152") and norms["B"] == 1, "Coefficient norms differ")

    modes = phase_modes(candidate,4)
    expected = set(range(210)) - {0,165}
    need(len(modes) == 210 and len(certificate["modes"]) == 208 and
         {m["index"] for m in certificate["modes"]} == expected, "Incomplete or duplicated mode set")
    need({p.name for p in records.glob("mode_*.json")} == {f"mode_{i:04d}.json" for i in expected},
         "Missing or additional mode record")
    need({p.name for p in records.glob("mode_*.leaves.jsonl.gz")} ==
         {f"mode_{i:04d}.leaves.jsonl.gz" for i in expected}, "Missing or additional full cover")

    R = F(103,100)
    Ar, Br = (sum(abs(F(c))*R**int(n) for n,c in candidate[name].items()) for name in ("P", "B"))
    total = [[F(0),F(0)] for _ in range(481)]
    tail, coefficient_count = F(0), 0
    minimum = [None]*4
    factor = R**-963/(1-R**-2)
    counts = dict(boundary=0, annulus=0, subdivisions=0)
    for mode in certificate["modes"]:
        i = mode["index"]
        key, weight = modes[i]
        raw_name = str(RECORDS / f"mode_{i:04d}.json")
        need(raw_name in checked, "Mode absent from integrity manifest")
        raw = load(root / raw_name)
        need(raw["index"] == i and raw["fingerprint"] == fingerprint, "Raw mode identity differs")
        need(mode["parameters_r_alpha_s_beta"] == raw["key"] == list(map(str,key)) and
             mode["weight"] == raw["weight"] == str(weight), "Mode parameters or weight differ")
        need(mode["coefficients_of_z_to_2j_plus_1"] == raw["coefficients"] and
             len(raw["coefficients"]) == 481, "Coefficient records differ")
        for j,value in enumerate(raw["coefficients"]):
            lo,hi = interval(value)
            lo,hi = (lo*weight,hi*weight) if weight >= 0 else (hi*weight,lo*weight)
            total[j][0] += lo
            total[j][1] += hi
            coefficient_count += 1
        need(raw["leaf_file"] == f"mode_{i:04d}.leaves.jsonl.gz", "Incorrect leaf filename")
        leaf = records / raw["leaf_file"]
        need((audit / mode["leaf_file"]).resolve() == leaf.resolve(), "Certificate leaf link differs")
        need(checked.get(str(leaf.relative_to(root))) == mode["leaf_sha256"] == raw["leaf_sha256"],
             "Complete compressed cover hash differs")
        lower = mode["positive_lower_bounds_D_N1_N2_Euler_distance"]
        need(lower == raw["cover"]["minimum_witnesses"], "Positive mode minima differ")
        ds = list(map(dyad,lower))
        need(len(ds) == 4 and all(x > y for x,y in zip(ds,[F(3,100),F(1,10),F(7,50),F(1,200000)])),
             "A stated conservative positivity bound fails")
        for j,v in enumerate(ds):
            minimum[j] = v if minimum[j] is None else min(minimum[j],v)
        for name in counts:
            n = raw["cover"]["counts"][name]
            need(type(n) is int and n >= 0, "Invalid reported leaf count")
            counts[name] += n
        polynomial = defaultdict(F)
        values = [F(candidate["theta"]),key[0],key[2],key[1],key[3]]
        for term in primitives["polynomials"]["P0"]:
            powers, c = term["powers"],F(term["coefficient"])
            need(len(powers) == 7 and all(type(n) is int and n >= 0 for n in powers), "Invalid primitive monomial")
            for value,exponent in zip(values,powers[2:]):
                c *= value**exponent
            polynomial[tuple(powers[:2])] += c
        U = sum(abs(c)*Ar**j*Br**k for (j,k),c in polynomial.items())
        need(U == F(mode["polynomial_numerator_upper"]), "Rational numerator upper bound differs")
        D,N1,N2,delta = ds
        radicand,scale = N1*N2*delta,10**35
        integer_root = isqrt(radicand.numerator*scale*scale//radicand.denominator)
        lower_root = F(integer_root,scale)
        need(lower_root > 0 and lower_root*lower_root <= radicand < F(integer_root+1,scale)**2,
             "Square-root isolation fails")
        need(lower_root == F(mode["rational_square_root_lower"]), "Recorded square root differs")
        M = F(22,7)*U/(4*D*lower_root)
        need(M == F(mode["rational_boundary_modulus_upper"]), "Rational boundary bound differs")
        tail += abs(weight)*M*factor

    # Polynomial division in the elementary positive integral proving pi < 22/7.
    quotient = [F(4),F(0),F(-4),F(0),F(5),F(-4),F(1)]
    product = [F(0)]*9
    for j,c in enumerate(quotient):
        product[j] += c
        product[j+2] += c
    product[0] -= 4
    need(product == [0,0,0,0,1,-4,6,-4,1] and
         sum(c/F(j+1) for j,c in enumerate(quotient)) == F(22,7), "Positive integral identity fails")
    L = total[0][0]
    Q = sum(max(abs(lo),abs(hi)) for lo,hi in total[1:])
    eta = F(candidate["eta"])
    need(eta >= 0, "Negative radial damping")
    damped = F(1) if eta == 0 else min(F(1),2*abs(F(candidate["radial_lam"]))/(5*eta))
    amplitude = abs(F(candidate["epsilon"]))+abs(F(candidate["radial_epsilon"]))*damped
    phase = (2*amplitude)**5/factorial(5)
    omission = sum(abs(modes[i][1])/((1+2*modes[i][0][1])*(1+2*modes[i][0][3])) for i in (0,165))
    need(omission == F(candidate["epsilon"])**4/96, "Omitted primitive cost differs")
    actual = [L,Q,tail,phase,omission]
    names = ("linear_coefficient_strict_lower", "finite_nonlinear_norm_strict_upper",
             "complete_scalar_tail_strict_upper", "complete_phase_tail_strict_upper",
             "omitted_primitive_cost_strict_upper")
    bounds = [F(certificate[k]) for k in names]
    need(actual[0] > bounds[0] and all(x < y for x,y in zip(actual[1:],bounds[1:])),
         "A simplified proof estimate fails")
    gamma = L-Q-tail-phase-omission
    closing = bounds[0]-sum(bounds[1:])
    theorem = F(71185999,100000000)
    need(gamma > closing == F("0.711859994603") > theorem == F(certificate["certified_gamma"]) >
         F(certificate["target_gamma"]) == F(3559299,5000000), "Exact final comparison fails")
    need(F("1.404770620695")*theorem > 1, "Strict decimal reciprocal comparison fails")
    need(coefficient_count == 100048 and counts == dict(boundary=1703936,annulus=1832072,subdivisions=128136),
         "Coefficient or reported covering counts differ")
    return dict(result="PASS", theorem="K_G^C <= 100000000/71185999 < 1.404770620695",
                integrity_files=len(checked), integrity_bytes=total_bytes,
                coefficient_records_compared=coefficient_count, full_cover_files_hashed=208,
                reported_cover_counts=counts, linear_lower=decimal_bound(L),
                finite_nonlinear_upper=decimal_bound(Q,True), complete_scalar_tail_upper=decimal_bound(tail,True),
                complete_phase_tail_upper=decimal_bound(phase,True), omitted_primitive_cost_upper=decimal_bound(omission,True),
                gamma_lower=decimal_bound(gamma), exact_proof_closing_bound=str(closing),
                exact_theorem_gamma=str(theorem), global_positive_minima=list(map(str,minimum)),
                full_geometric_tree_traversal_performed=False, local_Arb_enclosures_recomputed=False,
                scope="Integrity and exact coefficient aggregation; the local enclosures require the separate full recomputation described in README.md.",
                elapsed_seconds=round(time.monotonic()-started,3))


if __name__ == "__main__":
    if hasattr(sys,"set_int_max_str_digits"):
        sys.set_int_max_str_digits(40000)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify(args.root,args.manifest or args.root/"SHA256SUMS.json")
    rendered = json.dumps(result,indent=2)+"\n"
    if args.output is not None:
        with args.output.open("x",encoding="utf8") as stream:
            stream.write(rendered)
    print(rendered,end="")
