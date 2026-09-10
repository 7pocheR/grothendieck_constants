"""Complete outward-interval certificate for a rational radial-phase candidate.

The analytic proof is in certificate_argument.md and gaussian_expansion.md.
Every retained mode and boundary arc is covered; no sampled maximum is used.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
import json
from math import factorial
import os
from pathlib import Path
import time

from flint import arb, acb, arb_series, acb_poly, ctx


def rat(value):
    q = F(value)
    return arb(q.numerator)/q.denominator


def exact_modes(data, order):
    e0, e1 = F(data["epsilon"]), F(data["radial_epsilon"])
    l0, l1, eta, k = (F(data[x]) for x in ["lam", "radial_lam", "eta", "k"])
    assert eta >= 0
    atoms = []
    for n in range(order+1):
        for n1 in range(n+1):
            n0 = n-n1
            for j0 in range(n0+1):
                for j1 in range(n1+1):
                    c = ((e0/2)**n0*(e1/2)**n1*(-1)**(j0+j1)
                         /factorial(j0)/factorial(n0-j0)
                         /factorial(j1)/factorial(n1-j1))
                    if c:
                        atoms.append((n, c, k-(n0-2*j0)*l0-(n1-2*j1)*l1, n1*eta))
    modes = {}
    for i, (ni, ci, ri, ai) in enumerate(atoms):
        for j in range(i, len(atoms)):
            nj, cj, rj, aj = atoms[j]
            if ni+nj > order:
                break
            left, right = sorted([(ri, ai), (rj, aj)])
            key = (*left, *right)
            modes[key] = modes.get(key, F(0))+ci*cj*(1 if i == j else 2)
    return [(key, weight) for key, weight in sorted(modes.items()) if weight]


def polys(a, b, h, r, s, ar, ac):
    rr, cc = ar*ar+r*r, ac*ac+s*s
    fr, fc = rr+ar, cc+ac
    a2, b2, h2 = a*a, b*b, h*h
    D = ((1+ar)**2+r*r)*((1+ac)**2+s*s)-fr*fc*(a2+b2)+rr*cc*a2*b2+2*a*b*r*s
    N1 = ((1+h2)*(1+ar)*((1+ac)**2+s*s)
          -fc*((1+ar*(1+h2))*a2+(h2+ar*(1+h2))*b2)
          +ar*(1+h2)*cc*a2*b2)
    N2 = ((1+h2)*(1+ac)*((1+ar)**2+r*r)
          -fr*((1+ac*(1+h2))*a2+(h2+ac*(1+h2))*b2)
          +ac*(1+h2)*rr*a2*b2)
    V = (a*((1+ar)*(1+ac)+h2*r*s)-b*(h2*(1+ar)*(1+ac)+r*s)
         +a2*b*(h2*ar*ac+r*s)-a*b2*(ar*ac+h2*r*s))
    W = h*(a+b)*(r+s+(1-a*b)*(ac*r+ar*s))
    return D, N1, N2, V+W, V-W


def series(A, B, h, r, s, ar, ac, degree):
    y = arb_series([0, 1])
    rr, cc = ar*ar+r*r, ac*ac+s*s
    fr, fc = rr+ar, cc+ac
    h2, A2, B2 = h*h, A*A, B*B
    D = ((1+ar)**2+r*r)*((1+ac)**2+s*s)-fr*fc*y*(A2+B2)+rr*cc*y*y*A2*B2+2*y*A*B*r*s
    N1 = ((1+h2)*(1+ar)*((1+ac)**2+s*s)
          -fc*y*((1+ar*(1+h2))*A2+(h2+ar*(1+h2))*B2)
          +ar*(1+h2)*cc*y*y*A2*B2)
    N2 = ((1+h2)*(1+ac)*((1+ar)**2+r*r)
          -fr*y*((1+ac*(1+h2))*A2+(h2+ac*(1+h2))*B2)
          +ac*(1+h2)*rr*y*y*A2*B2)
    V = (A*((1+ar)*(1+ac)+h2*r*s)-B*(h2*(1+ar)*(1+ac)+r*s)
         +y*A2*B*(h2*ar*ac+r*s)-y*A*B2*(ar*ac+h2*r*s))
    W = h*(A+B)*(r+s+(1-y*A*B)*(ac*r+ar*s))
    Z = y*(V*V-W*W)/(N1*N2)
    assert Z[0] == 0
    coefficients = [arb(1)]
    for j in range(1, degree+1):
        coefficients.append(coefficients[-1]*rat(F((2*j-1)**2, 4*j*(j+1))))
    HG = arb_series([coefficients[-1]])
    for c in reversed(coefficients[:-1]):
        HG = HG*Z+c
    return arb.pi()/4*(V+W)/D/(N1*N2).sqrt()*HG


def polynomial_coefficients(data, key):
    vals = [arb(0)]*((max(map(int, data[key]))+1)//2)
    for n, coefficient in data[key].items():
        assert int(n) > 0 and int(n) % 2
        vals[(int(n)-1)//2] = rat(coefficient)
    return vals


def boundary(P, B, h, r, s, ar, ac, radius, panels):
    maximum = arb(0)
    minimum_D, minimum_NN = arb(10**12), arb(10**12)
    cap = rat("99/100")
    refinements = 0
    for j in range(panels):
        angle = arb.pi()*(arb(2*j+1)/panels+arb(0, arb(1)/panels))
        direction = acb(0, angle).exp()
        z = radius*direction
        a, b = z*P(z*z), z*B(z*z)
        D, N1, N2, pp, qq = polys(a, b, h, r, s, ar, ac)
        NN = N1*N2
        assert NN.real > 0, ("boundary N1N2", j, str(NN))
        Z = pp*qq/NN
        assert Z.real < cap, ("boundary Euler branch", j, str(Z))
        assert abs(D) > 0, ("boundary determinant", j, str(D))
        M = arb.pi()/4*abs(pp)/(abs(D)*abs(NN).sqrt())/(1-cap).sqrt()
        maximum = arb.max(maximum, M)
        minimum_D = arb.min(minimum_D, abs(D))
        minimum_NN = arb.min(minimum_NN, NN.real)
        stack = [(arb(1), radius, 0)]
        while stack:
            lo, hi, depth = stack.pop()
            radial = (lo+hi)/2+arb(0, (hi-lo)/2)
            zz = radial*direction
            aa, bb = zz*P(zz*zz), zz*B(zz*zz)
            DD = polys(aa, bb, h, r, s, ar, ac)[0]
            if abs(DD) > 0:
                continue
            assert depth < 12, ("annulus determinant", j, depth, str(DD))
            refinements += 1
            mid = (lo+hi)/2
            stack.extend([(lo, mid, depth+1), (mid, hi, depth+1)])
    return maximum, dict(boundary_modulus_bound=str(maximum),
                         minimum_D_boundary=str(minimum_D),
                         minimum_Re_N1N2_boundary=str(minimum_NN),
                         real_hypergeometric_argument_cap="99/100",
                         full_boundary_arcs=panels,
                         radial_subdivisions=refinements)


def one_mode(task):
    index, key, weight, data, settings, output, fingerprint = task
    path = Path(output)/f"mode_{index:04d}.json"
    if path.exists():
        old = json.loads(path.read_text())
        assert old["fingerprint"] == fingerprint
        assert old["key"] == key and old["weight"] == weight
        return old
    ctx.prec = settings["bits"]
    ctx.cap = settings["degree"]+1
    start = time.time()
    h = rat(data["theta"])
    r, ar, s, ac = map(rat, key)
    pc, bc = polynomial_coefficients(data, "P"), polynomial_coefficients(data, "B")
    cf = series(arb_series(pc), arb_series(bc), h, r, s, ar, ac, settings["degree"])
    radius = rat(settings["radius"])
    M, proof = boundary(acb_poly(pc), acb_poly(bc), h, r, s, ar, ac,
                        radius, settings["panels"])
    tail = M*radius**(-2*settings["degree"]-3)/(1-radius**(-2))
    out = dict(index=index, fingerprint=fingerprint, key=key, weight=weight,
               coefficients=[str(cf[j]) for j in range(settings["degree"]+1)],
               unweighted_scalar_tail=str(tail), proof=proof,
               elapsed_seconds=time.time()-start, status="PASS")
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(out, indent=2)+"\n")
    os.replace(temporary, path)
    print(json.dumps(dict(mode=index, elapsed=out["elapsed_seconds"], status="PASS")), flush=True)
    return out


def run(args):
    start = time.time()
    data = json.loads(Path(args.candidate).read_text())
    pnorm = sum(abs(F(x)) for x in data["P"].values())
    bnorm = sum(abs(F(x)) for x in data["B"].values())
    assert pnorm < 1 and bnorm <= 1
    assert F(data["target_gamma"]) > 0
    order = int(data["phase_order"])
    settings = dict(bits=args.bits, degree=args.degree, panels=args.panels,
                    radius=args.radius, phase_order=order)
    fingerprint = sha256(json.dumps(dict(data=data, settings=settings),
                                    sort_keys=True).encode()).hexdigest()
    dest = Path(args.output)
    dest.mkdir(parents=True, exist_ok=True)
    modes = exact_modes(data, order)
    tasks = [(i, list(map(str, key)), str(weight), data, settings,
              str(dest), fingerprint) for i, (key, weight) in enumerate(modes)]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        outputs = list(pool.map(one_mode, tasks, chunksize=1))
    ctx.prec = args.bits
    ctx.cap = args.degree+1
    total = [arb(0)]*(args.degree+1)
    tail = arb(0)
    for mode in outputs:
        weight = rat(mode["weight"])
        for j, value in enumerate(mode["coefficients"]):
            total[j] += weight*arb(value)
        tail += abs(weight)*arb(mode["unweighted_scalar_tail"])
    eta, lam = F(data["eta"]), abs(F(data["radial_lam"]))
    # |sin(2 lambda q)| <= lambda E, and e > 5/2.
    radial_bound = F(1) if eta == 0 else min(F(1), 2*lam/(5*eta))
    phase_bound = abs(F(data["epsilon"]))+abs(F(data["radial_epsilon"]))*radial_bound
    exact_phase_tail = (2*phase_bound)**(order+1)/factorial(order+1)
    nonlinear = sum((abs(x) for x in total[1:]), arb(0))
    gamma = total[0]-nonlinear-tail-rat(exact_phase_tail)
    target = rat(data["target_gamma"])
    result = dict(status="PASS" if gamma > target else "FAIL",
                  fingerprint=fingerprint, settings=settings,
                  exact_P_norm=str(pnorm), exact_B_norm=str(bnorm),
                  retained_primitive_count=len(modes),
                  linear_coefficient=str(total[0]),
                  finite_nonlinear_norm=str(nonlinear),
                  complete_scalar_degree_tail=str(tail),
                  exact_phase_multiplier_bound=str(phase_bound),
                  complete_phase_Taylor_tail=str(exact_phase_tail),
                  guaranteed_gamma=str(gamma), target_gamma=data["target_gamma"],
                  margin=str(gamma-target), reciprocal_target=str(1/target),
                  elapsed_seconds=time.time()-start,
                  arithmetic="Exact rational input and mode weights; outward Arb intervals; "
                             "complete arc/annulus coverings and both complete tails.")
    (dest/"certificate.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2), flush=True)
    assert gamma > target


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", default="radial_low_candidate.json")
    parser.add_argument("--output", default="radial_certificate")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--bits", type=int, default=256)
    parser.add_argument("--degree", type=int, default=420)
    parser.add_argument("--panels", type=int, default=4096)
    parser.add_argument("--radius", default="103/100")
    run(parser.parse_args())
