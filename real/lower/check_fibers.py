"""Independent, bounded replay of every returned affine-fiber leaf.

The input is treated as a proposed subdivision and proposed rational bounds.
No producer module, checkpoint status, or optimization result is imported.
The integral verifier uses midpoint Taylor ranges and Gaussian moment
recurrences; it does not use the producer's Bernstein or primitive routines.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import heapq
import itertools
import json
import math
import sys
import time

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
TARGETS = {
    "A": (Q(21, 20), Q(107, 105), Q(107, 2100), Q(3, 2), Q(26, 125)),
    "B": (Q(13, 10), Q(39, 35), Q(117, 350), Q(11, 5), Q(1, 6000)),
}
DOMAIN = ((Q(0), Q(2)), (Q(0), Q(2)), (Q(-2), Q(2)), (Q(-2), Q(2)))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def box_for(path):
    box = list(DOMAIN)
    for bit in path:
        if bit not in "01":
            raise ValueError("Nonbinary subdivision path")
        widths = [u-l for l, u in box]
        axis = widths.index(max(widths))
        l, u = box[axis]
        mid = (l+u)/2
        box[axis] = (l, mid) if bit == "0" else (mid, u)
    return box


def matrix_checks():
    a = Q(13, 8)
    evidence = []
    for name, (p, q, r2, t, kappa) in TARGETS.items():
        assert q == (p+a)/(a+1) and r2 == (p-1)*q
        for label, P, R2, S in (("even", p, r2, q),
                               ("odd", Q(139, 118), 6*Q(21, 118)**2, Q(63, 59))):
            for d0, d1 in ((0, 0), (1, 0), (-a, 1)):
                determinant = (P-d0)*(S-d1)-R2
                assert P-d0 >= 0 and S-d1 >= 0 and determinant >= 0
                evidence.append([name, label, str(d0), str(d1), str(determinant)])
            assert P > 0 and P*S-R2 > 0
            assert 2-P > 0 and 2-S > 0 and (2-P)*(2-S)-R2 > 0
    return evidence


def inventory(evidence):
    vertices = {}
    leaves = []
    digests = {}
    counts = {}
    for name in TARGETS:
        all_paths = []
        raw_lines = 0
        for bits in itertools.product("01", repeat=3):
            shard = "".join(bits)
            path = evidence / (name+"_"+shard) / "leaves.jsonl"
            digests[str(path.relative_to(evidence))] = sha(path)
            seen = {}
            for line in path.read_text().splitlines():
                row = json.loads(line)
                raw_lines += 1
                p = row["path"]
                assert p.startswith(shard)
                if p in seen:
                    assert row == seen[p], "Conflicting duplicate leaf"
                    continue
                seen[p] = row
                box = box_for(p)
                correction = sum(((u-l)**2/4 for l, u in box), Q(0))
                assert correction == Q(row["correction"])
                proposed = list(map(Q, row["vertex_upper"]))
                assert len(proposed) == 16
                assert Q(row["upper"]) == max(proposed)+correction
                assert Q(row["upper"]) <= TARGETS[name][3]+TARGETS[name][4]
                keys = []
                for v, bound in zip(itertools.product(*box), proposed):
                    key = name+":"+",".join(map(str, v))
                    keys.append(key)
                    if key not in vertices:
                        vertices[key] = dict(key=key, target=name, v=list(map(str, v)),
                                             proposed=str(bound))
                    else:
                        vertices[key]["proposed"] = str(min(bound, Q(vertices[key]["proposed"])))
                leaves.append(dict(target=name, path=p, keys=keys, correction=str(correction)))
            paths = sorted(seen)
            assert paths
            assert all(not q.startswith(p) for p, q in zip(paths, paths[1:])), "Overlapping leaves"
            assert sum((Q(1, 2**(len(p)-3)) for p in paths), Q(0)) == 1, "Incomplete shard"
            all_paths.extend(paths)
        assert sum((Q(1, 2**len(p)) for p in all_paths), Q(0)) == 1
        counts[name] = dict(raw_lines=raw_lines, unique_leaves=len(all_paths),
                            max_depth=max(map(len, all_paths)),
                            unique_vertices=sum(v["target"] == name for v in vertices.values()))
    payload = dict(input_sha256=digests, counts=counts,
                   matrices=matrix_checks(), vertices=[vertices[k] for k in sorted(vertices)], leaves=leaves)
    return payload


def make_verifier(precision, tolerance):
    from flint import arb, ctx
    ctx.prec = precision

    def aq(x):
        x = Q(x)
        return arb(x.numerator)/x.denominator

    nu = (arb(2)/arb.pi()).sqrt()
    sqrt2 = arb(2).sqrt()
    cutoff = Q(12)

    def phi(x):
        return (-x*x/2).exp()/((2*arb.pi()).sqrt())

    def mass(l, u):
        # erfc avoids subtracting two values almost equal to one in the tails.
        if l >= 0:
            return ((l/sqrt2).erfc()-(u/sqrt2).erfc())/2
        if u <= 0:
            return ((-u/sqrt2).erfc()-(-l/sqrt2).erfc())/2
        return ((u/sqrt2).erf()-(l/sqrt2).erf())/2

    def integral_signed(c, lo, hi):
        l, u = aq(lo), aq(hi)
        fl, fu = phi(l), phi(u)
        moments = [mass(l, u), fl-fu]
        for j in range(2, len(c)):
            moments.append(l**(j-1)*fl-u**(j-1)*fu+(j-1)*moments[j-2])
        return sum((coef*moment for coef, moment in zip(c, moments)), arb(0))

    def cell(c, lo, hi):
        mid, rad = aq((lo+hi)/2), aq((hi-lo)/2)
        # Exact Taylor identity at the midpoint, with triangle bounds for all
        # nonconstant terms; all arithmetic remains outward.
        taylor = [sum((c[j]*math.comb(j, k)*mid**(j-k)
                        for j in range(k, len(c))), arb(0)) for k in range(len(c))]
        radius = sum((taylor[k].abs_upper()*rad**k for k in range(1, len(c))), arb(0))
        lower = taylor[0].lower()-radius.upper()
        upper = taylor[0].upper()+radius.upper()
        if lower >= 0:
            return True, integral_signed(c, lo, hi)
        if upper <= 0:
            return True, -integral_signed(c, lo, hi)
        supremum = max(lower.abs_upper(), upper.abs_upper())
        return False, (supremum*mass(aq(lo), aq(hi)).upper()).upper()

    def absolute_integral(c):
        known = arb(0)
        unknown = arb(0)
        heap = []
        serial = 0

        def insert(lo, hi):
            nonlocal known, unknown, serial
            certain, value = cell(c, lo, hi)
            if certain:
                known += value
            else:
                unknown += value
                heapq.heappush(heap, (-float(value), serial, lo, hi, value))
                serial += 1

        # This fixed partition does not depend on any computed polynomial roots.
        for i in range(-12, 12):
            insert(Q(i), Q(i+1))
        splits = 0
        while heap and not unknown.upper() <= aq(tolerance):
            if splits == 20000:
                break  # The remainder is still included; early stopping is safe.
            _, _, lo, hi, old = heapq.heappop(heap)
            unknown -= old
            mid = (lo+hi)/2
            insert(lo, mid)
            insert(mid, hi)
            splits += 1
        R = aq(cutoff)
        e = (-R*R/2).exp()
        tails = [(R/sqrt2).erfc(), nu*e]
        for j in range(2, len(c)):
            tails.append(nu*R**(j-1)*e+(j-1)*tails[j-2])
        tail = sum((coef.abs_upper()*moment for coef, moment in zip(c, tails)), arb(0))
        upper = (known.upper()+unknown.upper()+tail.upper()).upper()
        return upper, dict(splits=splits, finite_remainder=str(unknown.upper()), tail=str(tail.upper()))

    def ceil_rational(x, bits=60):
        den = 1 << bits
        integer = (x.upper()*den).ceil().unique_fmpz()
        if integer is None:
            raise ArithmeticError("Unable to extract an exact upward integer")
        result = Q(int(integer), den)
        assert aq(result) >= x.upper()
        return result

    def vertex(task):
        p, q, r2, t, kappa = TARGETS[task["target"]]
        po, qo, ro2 = Q(139, 118), Q(63, 59), 6*Q(21, 118)**2
        v = list(map(Q, task["v"]))
        # Multiplication by a Cholesky factor reconstructed from its defining
        # matrices. Ordering is Hermite degrees 0,1,2,3.
        d0 = 2*aq(p).sqrt()*aq(v[0])
        d1 = 2*aq(po).sqrt()*aq(v[1])
        d2 = 2*(aq(r2/p).sqrt()*aq(v[0])+aq(q-r2/p).sqrt()*aq(v[2]))
        d3 = 2*(aq(ro2/po).sqrt()*aq(v[1])+aq(qo-ro2/po).sqrt()*aq(v[3]))
        c = [d0-d2/sqrt2, d1-3*d3/arb(6).sqrt(), d2/sqrt2, d3/arb(6).sqrt()]
        plus, minus = c.copy(), c.copy()
        plus[1] += aq(t)
        minus[1] -= aq(t)
        up1, info1 = absolute_integral(plus)
        up2, info2 = absolute_integral(minus)
        upper = ceil_rational((up1+up2)/(2*nu)-aq(sum((s*s for s in v), Q(0))))
        return dict(key=task["key"], upper=str(upper),
                    proposed_bound_verified=upper <= Q(task["proposed"]), integrals=[info1, info2])

    return vertex


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("inventory", "part", "finish"))
    ap.add_argument("--evidence", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--parts", type=int, default=16)
    ap.add_argument("--part", type=int, default=0)
    ap.add_argument("--max-vertices", type=int, default=100000)
    ap.add_argument("--seconds", type=float, default=1200)
    ap.add_argument("--precision", type=int, default=192)
    ap.add_argument("--tolerance", default="1/4503599627370496")
    opt = ap.parse_args()
    assert opt.parts > 0 and 0 <= opt.part < opt.parts
    assert opt.precision >= 160 and Q(opt.tolerance) > 0
    assert opt.seconds > 0 and opt.max_vertices > 0
    opt.output.mkdir(parents=True, exist_ok=True)
    payload = inventory(opt.evidence)
    input_digest = hashlib.sha256(json.dumps(payload["input_sha256"], sort_keys=True).encode()).hexdigest()
    if opt.mode == "inventory":
        result = dict(counts=payload["counts"], input_sha256=payload["input_sha256"],
                       input_digest=input_digest, matrices=payload["matrices"],
                       vertex_integrals_recomputed=False)
        (opt.output/"inventory.json").write_text(json.dumps(result, indent=2)+"\n")
        print(json.dumps(result["counts"]))
        return
    if opt.mode == "part":
        assigned = payload["vertices"][opt.part::opt.parts]
        vertex = make_verifier(opt.precision, Q(opt.tolerance))
        started = time.monotonic()
        result = dict(input_digest=input_digest, checker_sha256=sha(__file__),
                       part=opt.part, parts=opt.parts, assigned=len(assigned),
                       precision=opt.precision, tolerance=opt.tolerance, records=[])
        for task in assigned:
            if len(result["records"]) >= opt.max_vertices or time.monotonic()-started >= opt.seconds:
                break
            result["records"].append(vertex(task))
        result["complete"] = len(result["records"]) == len(assigned)
        result["elapsed_seconds"] = time.monotonic()-started
        result["all_proposed_bounds_verified"] = all(r["proposed_bound_verified"] for r in result["records"])
        (opt.output/f"part_{opt.part:02}.json").write_text(json.dumps(result, indent=2)+"\n")
        print(json.dumps({k:v for k,v in result.items() if k != "records"}))
        return
    records = {}
    for i in range(opt.parts):
        result = json.loads((opt.output/f"part_{i:02}.json").read_text())
        assert result["input_digest"] == input_digest
        assert result["checker_sha256"] == sha(__file__)
        assert result["part"] == i and result["parts"] == opt.parts
        # Do not trust the complete flag: compare the exact key sets.
        assigned = {v["key"] for v in payload["vertices"][i::opt.parts]}
        assert len(result["records"]) == len(assigned)
        assert {r["key"] for r in result["records"]} == assigned
        for row in result["records"]:
            assert row["key"] not in records
            records[row["key"]] = Q(row["upper"])
    minimum_gap = {t: None for t in TARGETS}
    for leaf in payload["leaves"]:
        upper = max(records[k] for k in leaf["keys"])+Q(leaf["correction"])
        target = leaf["target"]
        gap = TARGETS[target][3]+TARGETS[target][4]-upper
        assert gap >= 0, (target, leaf["path"], str(gap))
        minimum_gap[target] = gap if minimum_gap[target] is None else min(gap, minimum_gap[target])
    proposed_verified = sum(records[v["key"]] <= Q(v["proposed"]) for v in payload["vertices"])
    result = dict(input_digest=input_digest, checker_sha256=sha(__file__),
                   counts=payload["counts"], checked_vertices=len(records),
                   checked_leaves=len(payload["leaves"]),
                   proposed_vertex_bounds_verified=proposed_verified,
                   exact_minimum_leaf_gap={k:str(v) for k,v in minimum_gap.items()})
    (opt.output/"verified.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
