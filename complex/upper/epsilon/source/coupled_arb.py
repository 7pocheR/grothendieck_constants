"""Outward coefficients, complete boundary arcs, and Cauchy tails for six-coordinate modes."""
from fractions import Fraction as Q
import gzip
import json
from pathlib import Path
import signal
import time

from flint import arb, acb, arb_series, acb_poly, ctx
from coupled_exact import (NAMES, ZERO, candidate, encode_primitives, input_facts,
                            primitive, require)
from winding_exact_core import (arc_witness, binary_ball, ceil_dyadic, dyad, floor_dyadic)


def rat(value):
    value = Q(value)
    return arb(value.numerator)/value.denominator


def pack(value):
    require(value.is_finite(), "Nonfinite arithmetic enclosure")
    return dict(mid=list(map(int, value.mid().man_exp())), rad=list(map(int, value.rad().man_exp())))


def pack_geometry_real(value):
    require(value.is_finite(), "Nonfinite geometric enclosure")
    lo = dyad(list(map(int, value.lower().man_exp())))
    hi = dyad(list(map(int, value.upper().man_exp())))
    left, right = floor_dyadic(lo, 80)[0], ceil_dyadic(hi, 80)[0]
    return dict(mid=[left+right, -81], rad=[right-left, -81])


def rectangle(value):
    return dict(real=pack_geometry_real(value.real), imag=pack_geometry_real(value.imag))


def formal_coefficients(d, polynomials, degree):
    ctx.cap = degree+1
    maps = []
    for name in ("P", "B", "C"):
        values = [arb(0)]*(degree+1)
        for exponent, coefficient in d[name].items():
            j = (int(exponent)-1)//2
            if j <= degree:
                values[j] = rat(coefficient)
        maps.append(arb_series(values))
    y = arb_series([0, 1])
    powers = []
    for i, f in enumerate(maps):
        highest = max(k[i] for p in polynomials.values() for k in p)
        powers.append([f**j for j in range(highest+1)])
    result = {}
    for name, polynomial in polynomials.items():
        parity = int(name in ("P0", "Q0"))
        value = arb_series([0])
        for k, coefficient in polynomial.items():
            require(sum(k) % 2 == parity, "Wrong formal parity")
            term = rat(coefficient)*y**((sum(k)-parity)//2)
            for i in range(3):
                term *= powers[i][k[i]]
            value += term
        result[name] = value
    D, N1, N2, P0, Q0 = (result[name] for name in NAMES)
    T0 = polynomials["N1"][ZERO]*polynomials["N2"][ZERO]
    require(T0 > 0 and D[0] > 0, "Origin branch not positive")
    T = N1*N2
    Z = y*P0*Q0/T
    require(Z[0] == 0, "Nonzero hypergeometric constant")
    h = [Q(1)]
    for j in range(1, degree+1):
        h.append(h[-1]*Q((2*j-1)**2, 4*j*(j+1)))
    H = arb_series([rat(h[-1])])
    for coefficient in reversed(h[:-1]):
        H = H*Z+rat(coefficient)
    F = arb.pi()/4 * P0 / D / rat(T0).sqrt() / (T/rat(T0)).sqrt() * H
    return [pack(F[j]) for j in range(degree+1)]


def compiled(polynomials):
    return [[(k, rat(c)) for k, c in polynomials[name].items()] for name in NAMES]


def evaluate(polynomials, *values):
    bounds = [max(k[i] for p in polynomials for k, c in p) for i in range(3)]
    powers = [[value**j for j in range(bound+1)] for value, bound in zip(values, bounds)]
    result = []
    for polynomial in polynomials:
        value = acb(0)
        for (i, j, k), c in polynomial:
            value += c*powers[0][i]*powers[1][j]*powers[2][k]
        result.append(value)
    return result


def map_polynomials(d, R):
    result = []
    for name in ("P", "B", "C"):
        length = (max(map(int, d[name]))+1)//2
        values, derivative = [arb(0)]*length, [arb(0)]*length
        L2 = Q(0)
        for exponent, coefficient in d[name].items():
            exponent, coefficient = int(exponent), Q(coefficient)
            values[(exponent-1)//2] = rat(coefficient)
            derivative[(exponent-1)//2] = rat(exponent*coefficient)
            L2 += exponent**2*abs(coefficient)*R**exponent
        result.append((acb_poly(values), acb_poly(derivative), L2))
    return result


def point(angle, R):
    return acb(rat(R)) if angle in (Q(0), Q(2)) else rat(R)*acb(0, arb.pi()*rat(angle)).exp()


def map_arc(polynomial, left, right, R):
    P, derivative, L2 = polynomial
    z = point((left+right)/2, R)
    center = z*P(z*z)
    first = acb(0, 1)*z*derivative(z*z)
    width = Q(22, 7)*(right-left)/2
    displacement = arb(0, rat(width).upper())
    error = arb(0, rat(width**2*L2/2).upper())
    return center+first*displacement+acb(error, error)


def geometry(d, polynomials, R, degree, settings, stream):
    require(R > 1, "Cauchy radius must exceed one")
    deadline = time.monotonic()+settings["radius_seconds"]
    maps = map_polynomials(d, R)
    polys = compiled(polynomials)
    endpoints = {}
    def endpoint(angle):
        if angle == 2:
            angle = Q(0)
        if angle not in endpoints:
            z = point(angle, R)
            endpoints[angle] = evaluate(polys, *(z*P(z*z) for P, D, L in maps))
        return endpoints[angle]
    maximum = Q(0)
    sums = [[Q(0), Q(0)] for _ in range(3)]
    minima = [None]*4
    leaves = refinements = evaluations = 0
    panels = settings["panels"]
    for panel in range(panels):
        stack = [(Q(2*panel, panels), Q(2*(panel+1), panels), 0)]
        while stack:
            require(time.monotonic() < deadline, "Radius work limit reached")
            require(evaluations < settings["max_evaluations"], "Radius evaluation limit reached")
            left, right, depth = stack.pop()
            values = evaluate(polys, *(map_arc(P, left, right, R) for P in maps))
            evaluations += 1
            try:
                Z = values[3]*values[4]/(values[1]*values[2])
                first, last = endpoint(left), endpoint(right)
                ratios = [b/a for a, b in zip(first[:3], last[:3])]
                encoded, encoded_Z, encoded_ratios = list(map(rectangle, values)), rectangle(Z), list(map(rectangle, ratios))
                witness = arc_witness(encoded, encoded_Z, encoded_ratios)
            except (ValueError, ZeroDivisionError):
                require(depth < settings["max_depth"], "Complete arc unresolved at subdivision limit")
                middle = (left+right)/2
                stack.extend([(middle, right, depth+1), (left, middle, depth+1)])
                refinements += 1
                continue
            leaves += 1
            require(leaves <= settings["max_leaves"], "Complete arc count exceeds budget")
            stream.write(json.dumps(dict(left=str(left), right=str(right), depth=depth,
                values=encoded, Z=encoded_Z, endpoint_ratios=encoded_ratios,
                exact_witness=witness), separators=(",", ":"))+"\n")
            maximum = max(maximum, dyad(witness["modulus_upper"]))
            for j, pair in enumerate(witness["lower"]):
                value = dyad(pair)
                minima[j] = value if minima[j] is None else min(minima[j], value)
            for j, (lo, hi) in enumerate(witness["argument_increments"]):
                sums[j][0] += dyad(lo)
                sums[j][1] += dyad(hi)
    require(all(-3 < lo <= hi < 3 for lo, hi in sums), "Zero winding not proved")
    require(leaves == panels+refinements, "Subdivision accounting differs")
    tail = maximum*R**(-2*degree-3)/(1-R**(-2))
    return dict(radius=str(R), boundary_modulus_upper=str(maximum), scalar_tail_upper=str(tail),
        minima=list(map(str, minima)), argument_sums=[[str(lo), str(hi)] for lo, hi in sums],
        leaves=leaves, refinements=refinements, evaluations=evaluations,
        map_second_angular_derivative_bounds=[str(p[2]) for p in maps])


def suggested_degree(weight, maximum, radius, budget, cap=16384):
    def cost(n):
        return abs(weight)*maximum*radius**(-2*n-3)/(1-radius**(-2))
    left, right = 0, 1
    while right < cap and cost(right) > budget:
        right = min(cap, 2*right)
    if cost(right) > budget:
        return None
    while left < right:
        middle = (left+right)//2
        if cost(middle) <= budget:
            right = middle
        else:
            left = middle+1
    return left


def one_mode(task):
    row, d, settings, output, binding = task
    import resource
    resource.setrlimit(resource.RLIMIT_AS, (settings["worker_memory_bytes"], settings["worker_memory_bytes"]))
    ctx.prec = settings["bits"]
    started = time.monotonic()
    directory = Path(output)/row["label"]
    directory.mkdir(exist_ok=False)
    record = dict(row, source_binding=binding, settings=settings, degree=settings["degree"],
                  status="UNRESOLVED", radius_attempts=[])
    def alarm(signum, frame):
        raise TimeoutError("Mode wall-clock limit reached")
    signal.signal(signal.SIGALRM, alarm)
    signal.alarm(settings["mode_seconds"])
    try:
        polynomials = primitive(d, row["key"])
        record["primitives"] = encode_primitives(polynomials)
        for radius in settings["radii"]:
            name = "geometry_"+radius.replace("/", "_")+".jsonl.gz"
            try:
                with gzip.open(directory/name, "wt", encoding="utf-8") as stream:
                    proof = geometry(d, polynomials, Q(radius), settings["degree"], settings, stream)
                record.update(geometry=proof, geometry_file=name)
                record["coefficients"] = formal_coefficients(d, polynomials, settings["degree"])
                width = sum(hi-lo for lo,hi in map(binary_ball, record["coefficients"]))
                record["weighted_finite_coefficient_interval_width"] = str(abs(Q(row["weight"]))*width)
                record["suggested_degree_for_per_mode_tail"] = suggested_degree(Q(row["weight"]),
                    Q(proof["boundary_modulus_upper"]), Q(radius), Q(settings["projection_tail_budget"]))
                record["status"] = "MODE_ARITHMETIC_PASS"
                break
            except (ValueError, ZeroDivisionError) as error:
                record["radius_attempts"].append(dict(radius=radius, reason=str(error), file=name,
                    qualification="Failure to establish a sufficient test; not evidence of a singularity."))
    except (Exception, MemoryError) as error:
        record["error"] = type(error).__name__+": "+str(error)
    finally:
        signal.alarm(0)
    record["elapsed_seconds"] = time.monotonic()-started
    (directory/"mode.json").write_text(json.dumps(record, indent=2)+"\n")
    return {k: record.get(k) for k in ("label", "key", "weight", "status", "error", "degree",
        "suggested_degree_for_per_mode_tail", "weighted_finite_coefficient_interval_width", "elapsed_seconds")}
