"""Independent replay of the three scalar supports used on [1/10,11/50].

No producer code or earlier checker result is imported. Tangent integrals
use a four-node Gaussian quadrature with an analytic eighth-derivative
remainder. Positive-part cells use Taylor ranges and monomial moments.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import time

sys.dont_write_bytecode = True


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_cases(campaign):
    new = campaign/"cluster/resident_20260907/results/real_lower_two_correlations_certificate"
    old = campaign/"round3/real_lower_frontier_continuation/agreement_degree9_certificate_01.json"
    cases = []
    for label, path, interval in (("degree9", old, ("1/10", "4/25")),
                                  ("new0", new/"scalar_0.json", ("4/25", "1/5")),
                                  ("new1", new/"scalar_1.json", ("1/5", "11/50"))):
        raw = json.loads(path.read_text())
        data = raw["input"].copy()
        if label != "degree9":
            candidate = new/("candidate_"+label[-1]+".json")
            assert json.loads(candidate.read_text()) == data
            assert raw["input_sha256"] == digest(candidate)
        else:
            data["tau"] = data["rho"]
            data["r"] = "1/10"
        assert Q(data["a"]) == Q(13, 8)
        assert 0 < Q(data["rho"]) < 1 and 0 < Q(data["tau"]) < 1
        assert 0 < Q(data["alpha"]) <= 1
        assert Q(data["d0"]) >= 0 and Q(data["d1"]) >= 0
        degrees = set(map(int, data["z"]))
        assert degrees == set(range(3, 10 if label == "degree9" else 18, 2))
        cells = [(Q(r["lo"]), Q(r["hi"])) for r in raw["leaves"]]
        assert cells and cells[0][0] == 0 and cells[-1][1] == 12
        assert all(l < u for l, u in cells)
        assert all(a[1] == b[0] for a, b in zip(cells, cells[1:]))
        cases.append(dict(label=label, data=data, cells=cells, interval=interval,
                          sha256=digest(path), file=str(path)))
    return cases


def arithmetic():
    from flint import arb, ctx
    ctx.prec = 224

    def aq(x):
        x = Q(x)
        return arb(x.numerator)/x.denominator

    def up(x, bits=110):
        den = 1 << bits
        n = (x.upper()*den).ceil().unique_fmpz()
        assert n is not None
        q = Q(int(n), den)
        assert aq(q) >= x.upper()
        return q

    def down(x, bits=110):
        return -up(-x, bits)

    return arb, aq, up, down


def noise_certificate(case):
    arb, aq, up, down = arithmetic()
    d = case["data"]
    tau, r = aq(d["tau"]), aq(d["r"])
    assert 0 < tau < aq("99/100") and 0 <= r <= aq("1/8")
    angle = tau.asin()
    assert angle < aq("143/100") and arb.pi() > aq("157/50")
    # On the complex disks of radius 1/50 about this real interval,
    # |1+sin(z)| = cosh(Im z)+sin(Re z) > 3/625.
    # Thus |exp(-r^2/(1+sin(z)))| < exp(4) < 55.
    assert 1-aq("3/25").cos() > aq("3/625")
    assert r*r/aq("3/625") < 4 and arb(4).exp() < 55
    sqrt30 = arb(30).sqrt()
    nodes = []
    for sign in (-1, 1):
        square = (3+sign*2*sqrt30/5)/7
        weight = (18-sign*sqrt30)/36
        root = square.sqrt()
        nodes.extend([(-root, weight), (root, weight)])
    N = 1024
    step = 2*angle/N
    total = arb(0)
    panels = []
    for i in range(N):
        center = -angle+(arb(i)+aq("1/2"))*step
        value = arb(0)
        for node, weight in nodes:
            theta = center+step*node/2
            value += weight*(-r*r/(1+theta.sin())).exp()
        contribution = step*value/(2*arb.pi())
        total += contribution
        panels.append([str(down(contribution)), str(up(contribution))])
    # For a cell of length h, integral of the squared monic node
    # polynomial is h^9/44100. Hermite interpolation and Cauchy give
    # total error <= (L/pi)*h^8*55*50^8/44100.
    error = Q(55*50**8, 44100)*Q(3, N)**8  # L<3, pi>3
    value = total+arb(0, aq(error).upper())
    c = ((1-tau)/(1+tau)).sqrt()
    rstar2 = -4*c.log()/(1/(c*c)-c*c)
    mass = (r/arb(2).sqrt()).erfc()
    slope = (r/(c*arb(2).sqrt())).erf()-(c*r/arb(2).sqrt()).erf()
    intercept = value-slope*mass
    assert r*r <= rstar2
    assert intercept >= 0 and slope >= 0
    assert aq(d["d0"]) >= intercept and aq(d["d1"]) >= slope
    return dict(label=case["label"], tau=d["tau"], r=d["r"], panels=panels,
                F_lower=str(down(value)), F_upper=str(up(value)),
                quadrature_error=str(error),
                turning_radius_squared_lower=str(down(rstar2)),
                turning_radius_squared_upper=str(up(rstar2)),
                slope_lower=str(down(slope)), slope_upper=str(up(slope)),
                intercept_lower=str(down(intercept)), intercept_upper=str(up(intercept)))


def scalar_calculator(case):
    arb, aq, up, down = arithmetic()
    d = case["data"]
    nu = (2/arb.pi()).sqrt()
    threshold = (1-aq(d["alpha"]))*aq(d["d1"])/(nu*aq(d["tau"])**3)
    degree = max(map(int, d["z"]))
    hermite_coefficients = {1: aq(d["ell"]),
                            **{int(j):2*aq(v) for j, v in d["z"].items()}}
    p = [arb(0) for _ in range(degree+1)]
    for j, coef in hermite_coefficients.items():
        for k in range(j//2+1):
            integer = Q((-1)**k*math.factorial(j), 2**k*math.factorial(k)*math.factorial(j-2*k))
            p[j-2*k] += coef*aq(integer)/arb(math.factorial(j)).sqrt()
    sqrt2 = arb(2).sqrt()

    def phi(x):
        return nu*(-x*x/2).exp()/2

    def moments(lo, hi):
        l, u = aq(lo), aq(hi)
        fl, fu = phi(l), phi(u)
        out = [((l/sqrt2).erfc()-(u/sqrt2).erfc())/2, fl-fu]
        for j in range(2, degree+1):
            out.append(l**(j-1)*fl-u**(j-1)*fu+(j-1)*out[j-2])
        return out

    def cell(index):
        lo, hi = case["cells"][index]
        m, radius = aq((lo+hi)/2), aq((hi-lo)/2)
        # Closed polynomial formula, with Taylor coefficients evaluated
        # by Horner's rule separately for each derivative order.
        taylor = []
        for k in range(degree+1):
            v = arb(0)
            for j in range(degree, k-1, -1):
                v = v*m+p[j]*math.comb(j, k)
            taylor.append(v)
        remainder = arb(0)
        for k in range(degree, 0, -1):
            remainder = (remainder+taylor[k].abs_upper())*radius
        lower = taylor[0].lower()-remainder.upper()
        upper = taylor[0].upper()+remainder.upper()
        supremum = max(lower.abs_upper(), upper.abs_upper())
        if supremum <= threshold.lower():
            bound, kind = Q(0), "zero"
        elif lower > threshold.upper() or upper < -threshold.upper():
            ms = moments(lo, hi)
            integral = sum((coef*moment for coef, moment in zip(p, ms)), arb(0))
            sign = 1 if lower > threshold.upper() else -1
            bound = up(sign*integral-threshold*ms[0])
            kind = "positive" if sign == 1 else "negative"
        else:
            mass = ((aq(lo)/sqrt2).erfc()-(aq(hi)/sqrt2).erfc())/2
            bound = up((supremum-threshold.lower()).max(0)*mass.upper())
            kind = "unresolved"
        assert bound >= 0
        return dict(case=case["label"], index=index, upper=str(bound), kind=kind)

    R = arb(12)
    e = (-R*R/2).exp()
    tails = [(R/sqrt2).erfc(), nu*e]
    for j in range(2, degree+1):
        tails.append(nu*R**(j-1)*e+(j-1)*tails[j-2])
    tail = up(sum((coef.abs_upper()*moment for coef, moment in zip(p, tails)), arb(0)))

    def complete(finite_upper):
        integral = 2*finite_upper+tail
        rho, alpha, tau = map(aq, (d["rho"], d["alpha"], d["tau"]))
        constant = alpha*rho.asin()/rho**3+(1-alpha)*(aq(d["d0"])+aq(d["d1"]))/((nu*nu)*tau**3)
        for j, value in d["z"].items():
            j = int(j)
            constant += aq(value)**2/(1+alpha if j == 3 else alpha*rho**(j-3))
        constant_upper = up(constant+aq(integral)/nu)
        quadratic = Q(d["a"])-Q(d["alpha"])/Q(d["rho"])**2
        assert quadratic > 0
        endpoints = [quadratic*Q(x)**2+Q(d["ell"])*Q(x)+constant_upper for x in case["interval"]]
        assert max(endpoints) < Q(2917, 2000)
        return dict(label=case["label"], integral_upper=str(integral), tail_upper=str(tail),
                    constant_upper=str(constant_upper), quadratic=str(quadratic),
                    interval=case["interval"], endpoint_upper=list(map(str, endpoints)))

    return cell, complete


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("noise", "part", "finish"))
    ap.add_argument("--campaign", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--part", type=int, default=0)
    ap.add_argument("--parts", type=int, default=16)
    ap.add_argument("--max-cells", type=int, default=1000000)
    ap.add_argument("--seconds", type=float, default=1200)
    opt = ap.parse_args()
    assert 0 <= opt.part < opt.parts and opt.parts > 0
    assert opt.seconds > 0 and opt.max_cells > 0
    cases = load_cases(opt.campaign)
    binding = {c["label"]:c["sha256"] for c in cases}
    checker_hash = digest(__file__)
    opt.output.mkdir(parents=True, exist_ok=True)
    if opt.mode == "noise":
        start = time.monotonic()
        result = dict(binding=binding, checker_sha256=checker_hash,
                       certificates=[noise_certificate(c) for c in cases])
        result["elapsed_seconds"] = time.monotonic()-start
        (opt.output/"noise.json").write_text(json.dumps(result, indent=2)+"\n")
        print(json.dumps(dict(noise_supports=3, elapsed_seconds=result["elapsed_seconds"])))
        return
    tasks = [(ci, j) for ci, case in enumerate(cases) for j in range(len(case["cells"]))]
    if opt.mode == "part":
        assigned = tasks[opt.part::opt.parts]
        calculators = [scalar_calculator(c)[0] for c in cases]
        start = time.monotonic()
        result = dict(binding=binding, checker_sha256=checker_hash,
                       part=opt.part, parts=opt.parts, assigned=len(assigned), records=[])
        for ci, j in assigned:
            if len(result["records"]) >= opt.max_cells or time.monotonic()-start >= opt.seconds:
                break
            result["records"].append(calculators[ci](j))
        result["elapsed_seconds"] = time.monotonic()-start
        result["complete"] = len(result["records"]) == len(assigned)
        (opt.output/f"scalar_part_{opt.part:02}.json").write_text(json.dumps(result, indent=2)+"\n")
        print(json.dumps({k:v for k,v in result.items() if k != "records"}))
        return
    noise = json.loads((opt.output/"noise.json").read_text())
    assert noise["binding"] == binding and noise["checker_sha256"] == checker_hash
    assert [x["label"] for x in noise["certificates"]] == [c["label"] for c in cases]
    for case, row in zip(cases, noise["certificates"]):
        d = case["data"]
        assert row["tau"] == d["tau"] and row["r"] == d["r"]
        assert Q(d["r"])**2 <= Q(row["turning_radius_squared_lower"])
        assert Q(row["intercept_lower"]) >= 0 and Q(row["slope_lower"]) >= 0
        assert Q(d["d0"]) >= Q(row["intercept_upper"])
        assert Q(d["d1"]) >= Q(row["slope_upper"])
    values = {}
    for i in range(opt.parts):
        part = json.loads((opt.output/f"scalar_part_{i:02}.json").read_text())
        assert part["binding"] == binding and part["checker_sha256"] == checker_hash
        assert part["part"] == i and part["parts"] == opt.parts
        expected = {(cases[ci]["label"], j) for ci,j in tasks[i::opt.parts]}
        actual = {(r["case"], r["index"]) for r in part["records"]}
        assert actual == expected and len(part["records"]) == len(expected)
        for r in part["records"]:
            key = r["case"], r["index"]
            assert key not in values and Q(r["upper"]) >= 0
            values[key] = Q(r["upper"])
    results = []
    for case in cases:
        total = sum((values[(case["label"], j)] for j in range(len(case["cells"]))), Q(0))
        results.append(scalar_calculator(case)[1](total))
    arb, aq, up, down = arithmetic()
    rho = aq("97/100")
    joint = up(rho.asin()/rho**3+(aq("13/8")-1/rho**2)*aq("1/10")**2)
    assert Q(13,8)-1/Q(97,100)**2 > 0 and joint < Q(2917,2000)
    result = dict(binding=binding, checker_sha256=checker_hash, checked_cells=len(values),
                   scalar_supports=results, joint_interval=["0","1/10"], joint_upper=str(joint))
    (opt.output/"scalar_verified.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
