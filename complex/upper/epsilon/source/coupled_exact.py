"""Exact six-coordinate primitives and finite phase modes; standard library only."""
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from math import factorial
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "source_snapshot"))
from winding_exact_core import load, require

ZERO = (0, 0, 0)
NAMES = ("D", "N1", "N2", "P0", "Q0")


def candidate():
    c = load(HERE / "source_snapshot/rational_candidate.json")
    source = load(HERE / "source_snapshot/inputs.json")
    expected = sha256(json.dumps(dict(data=source["data"], epsilon=source["epsilon"],
        settings=source["settings"], hashes=source["source_hashes"]), sort_keys=True).encode()).hexdigest()
    require(c["fingerprint"] == source["fingerprint"] == expected, "Candidate binding differs")
    d = {k: source["data"][k] for k in ("theta", "k", "damping_w", "damping_x", "profiles")}
    d.update({k: c[k] for k in ("P", "B", "C", "epsilon", "mixing", "lambda_value")})
    d.update(old_order=4, new_order=4, damping_lower="999/1000")
    require(Q(d["epsilon"]) == Q(source["epsilon"]), "Epsilon binding differs")
    require(Q(d["mixing"]) == Q(source["settings"]["mix_ratio"]) * Q(d["epsilon"]),
            "Mixing binding differs")
    require(source["settings"]["old_order"] == source["settings"]["new_order"] == 4,
            "Phase truncation binding differs")
    require((Q(d["epsilon"]), Q(d["mixing"]), Q(d["lambda_value"])) ==
            (Q(1,100), Q(-1,50), Q(1,4)), "Different finite-amplitude candidate")
    return d


def input_facts(d):
    norms = {}
    for name in ("P", "B", "C"):
        require(all(int(n) > 0 and int(n) % 2 for n in d[name]), "Nonodd correlation map")
        norms[name] = sum(abs(Q(c)) for c in d[name].values())
        require(norms[name] <= 1, "Inadmissible correlation map")
    require(norms["P"] < 1 and d["C"] == {"7": "1"}, "Different endpoint hypotheses")
    lower, dw, dx = map(Q, (d["damping_lower"], d["damping_w"], d["damping_x"]))
    require(lower > 0 and dw > 0 and dx > 0 and dw*dx >= lower*lower, "Invalid damping lower bound")
    M = Q(0)
    for e, l, u in d["profiles"]:
        e, l, u = Q(e), abs(Q(l)), Q(u)
        require(u >= 0, "Negative damping")
        M += abs(e) * (1 if u == 0 else min(Q(1), 2*l/(5*u*lower)))
    L, N = d["old_order"], d["new_order"]
    error = (2*M)**(L+1)/factorial(L+1)
    error += sum((2*M)**j/factorial(j) for j in range(L+1)) * (2*abs(Q(d["epsilon"])))**(N+1)/factorial(N+1)
    b1, b3, b11 = (Q(d["B"][str(n)]) for n in (1, 3, 11))
    require(b1 < 0 and b3 > 0 and b11 < 0 and b3 >= 25*abs(b11), "B endpoint lemma does not apply")
    b_deficit = 2*abs(b1*b11)
    return dict(norms={k: str(v) for k, v in norms.items()}, M=str(M),
                phase_error=str(error), B_unit_circle_modulus_upper=str(1-b_deficit),
                B_unit_circle_squared_deficit=str(2*b_deficit))


def scalar(c):
    c = Q(c)
    return {ZERO: c} if c else {}


def add(a, b):
    out = dict(a)
    for k, c in b.items():
        out[k] = out.get(k, Q(0)) + c
        if not out[k]:
            del out[k]
    return out


def scale(a, c):
    c = Q(c)
    return {k: c*v for k, v in a.items()} if c else {}


def mul(a, b):
    out = defaultdict(Q)
    for i, c in a.items():
        for j, e in b.items():
            out[tuple(x+y for x, y in zip(i, j))] += c*e
    return {k: c for k, c in out.items() if c}


def determinant(M):
    """Division-free subset expansion of the Leibniz determinant."""
    n = len(M)
    dp = {0: scalar(1)}
    for row in range(n):
        following = {}
        for mask, value in dp.items():
            for column in range(n):
                if mask & (1 << column) or not M[row][column]:
                    continue
                inversions = (mask >> (column+1)).bit_count()
                term = scale(mul(value, M[row][column]), (-1)**inversions)
                key = mask | (1 << column)
                following[key] = add(following.get(key, {}), term)
        dp = following
    return dp.get((1 << n)-1, {})


def matrices(d, key):
    """Coordinate order W,X,W',X',Y,Y'; key r,u,tau,s,v,sigma."""
    r, u, tau, s, v, sigma = map(Q, key)
    h, m, dw, dx = (Q(d[k]) for k in ("theta", "mixing", "damping_w", "damping_x"))
    S = [[scalar(int(i == j)) for j in range(6)] for i in range(6)]
    for (i, j), powers in zip(((0, 2), (1, 3), (4, 5)), ((1, 0, 0), (0, 1, 0), (0, 0, 1))):
        S[i][j] = S[j][i] = {powers: Q(1)}
    J = [[Q(0) for j in range(6)] for i in range(6)]
    for i, damping in ((0, u*dw), (1, u*dx), (2, v*dw), (3, v*dx)):
        J[i][i] = -damping
    for i, j, frequency in ((0, 1, r), (2, 3, s), (0, 4, tau), (2, 5, sigma)):
        J[i][j], J[j][i] = -frequency, frequency
    M = []
    for i in range(6):
        row = []
        for j in range(6):
            entry = scalar(int(i == j))
            for k in range(6):
                entry = add(entry, scale(S[i][k], -J[k][j]))
            row.append(entry)
        M.append(row)
    U = [[Q(1), h, Q(0), Q(0), m, Q(0)], [Q(0), Q(0), Q(1), -h, Q(0), -m]]
    return S, J, M, U


def primitive(d, key):
    S, J, M, U = matrices(d, key)
    out = {"D": determinant(M)}
    for name, i, j in (("N1", 0, 0), ("N2", 1, 1), ("P0", 0, 1), ("Q0", 1, 0)):
        Su = []
        for row in S:
            value = {}
            for p, c in zip(row, U[j]):
                value = add(value, scale(p, c))
            Su.append(value)
        bordered = [row + [value] for row, value in zip(M, Su)]
        bordered.append([scalar(c) for c in U[i]] + [{}])
        out[name] = scale(determinant(bordered), -1)
    for name, polynomial in out.items():
        parity = int(name in ("P0", "Q0"))
        require(all(sum(k) % 2 == parity for k in polynomial), "Primitive parity differs")
    require(all(out[name].get(ZERO, 0) > 0 for name in ("D", "N1", "N2")), "Incorrect origin branch")
    return out


def evaluate(poly, a, b, c):
    return sum(coefficient*a**i*b**j*c**k for (i, j, k), coefficient in poly.items())


def phase_powers(profiles, order):
    generators = [(Q(e)*sign/2, Q(l)*sign, Q(u))
                  for e, l, u in profiles for sign in (-1, 1)]
    powers = [{(Q(0), Q(0)): Q(1)}]
    for n in range(1, order+1):
        row = defaultdict(Q)
        for (frequency, damping), weight in powers[-1].items():
            for w, f, u in generators:
                row[frequency+f, damping+u] += weight*w/n
        powers.append({k: w for k, w in row.items() if w})
    return powers


def phase_pairs(profiles, k, order, unordered):
    """Retain TOTAL row-plus-column order, with exact rational weights."""
    powers = phase_powers(profiles, order)
    out = defaultdict(Q)
    for i in range(order+1):
        for j in range(order+1-i):
            for (f, u), w in powers[i].items():
                for (g, v), z in powers[j].items():
                    left, right = (Q(k)-f, u), (Q(k)-g, v)
                    if unordered:
                        left, right = sorted((left, right))
                    out[(*left, *right)] += w*z
    return {key: value for key, value in out.items() if value}


def old_and_new(d):
    old = phase_pairs(d["profiles"], d["k"], d["old_order"], True)
    new = phase_pairs([[d["epsilon"], d["lambda_value"], "0"]], "0", d["new_order"], False)
    return old, new


def canonical(key):
    key = tuple(map(Q, key))
    left, right = sorted((key[:3], key[3:]))
    return (*left, *right)


def new_for_old(old_key, old_weight, new):
    r, u, s, v = old_key
    out = defaultdict(Q)
    for (tau, _, sigma, __), weight in new.items():
        out[canonical((r, u, tau, s, v, sigma))] += old_weight*weight
    return {k: w for k, w in out.items() if w}


def complete_modes(d):
    old, new = old_and_new(d)
    # Different old keys cannot give the same full canonical key.
    for key, weight in sorted(old.items()):
        yield from sorted(new_for_old(key, weight, new).items())


def omission_cost(d, key, weight):
    lower = Q(d["damping_lower"])
    return abs(weight)/((1+2*Q(key[1])*lower)*(1+2*Q(key[4])*lower))


def selected_modes(d, labels):
    k, lam = Q(d["k"]), Q(d["lambda_value"])
    specifications = {
        "base": (k, 0, 0, k, 0, 0),
        "new_one_side": (k, 0, -lam, k, 0, 0),
        "new_opposite": (k, 0, -lam, k, 0, lam),
        "new_same": (k, 0, -lam, k, 0, -lam),
        "damped_asymmetric": (k-Q(1, 8), 1, -lam, k, 0, 0),
        "damped_exchanged": (k-Q(1, 8), 1, 0, k, 0, -lam),
        "frequency_two": (k-2, 0, -lam, k, 0, 0),
        "frequency_four": (k-4, 0, 0, k, 0, lam),
    }
    old, new = old_and_new(d)
    rows = []
    for label in labels:
        require(label in specifications, "Unknown representative mode")
        key = canonical(specifications[label])
        r, u, tau, s, v, sigma = key
        old_key = (r, u, s, v)
        require(old_key in old, "Representative old mode absent")
        weight = new_for_old(old_key, old[old_key], new)[key]
        rows.append(dict(label=label, key=list(map(str, key)), weight=str(weight),
                         omission_cost=str(omission_cost(d, key, weight))))
    return rows


def encode_primitives(polynomials):
    return {name: [[list(k), str(v)] for k, v in sorted(p.items())] for name, p in polynomials.items()}


def decode_primitives(encoded):
    require(set(encoded) == set(NAMES), "Missing primitive")
    for name, rows in encoded.items():
        require(len(rows) == len({tuple(k) for k, v in rows}), "Duplicate polynomial monomial")
        require(all(len(k) == 3 and all(type(j) is int and j >= 0 for j in k)
                    and Q(v) != 0 for k,v in rows), "Invalid polynomial monomial")
    return {name: {tuple(k): Q(v) for k, v in rows} for name, rows in encoded.items()}
