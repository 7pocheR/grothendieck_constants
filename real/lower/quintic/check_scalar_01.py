"""Rational exponential, agreement-support, envelope and finite-ratio checks.

Only the Python standard library is used. Producer integral values and envelope
choices are not used to establish an inequality.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import argparse
import hashlib
import json
import time

A = F(33, 20)
LAMBDA = F(3, 20)
C = F(3681, 2500)
BITS = 224
SCALE = 2**BITS


def read(p):
    return json.loads(Path(p).read_text())


def save(p, obj):
    assert not Path(p).exists(), p
    Path(p).write_text(json.dumps(obj, indent=2)+'\n')


def down(q):
    return F(q.numerator*SCALE//q.denominator, SCALE)


def up(q):
    return F(-((-q.numerator*SCALE)//q.denominator), SCALE)


@lru_cache(maxsize=None)
def exp_negative(q):
    """A rational enclosure of exp(-q), q >= 0, including very large q."""
    q = F(q)
    assert q >= 0
    z = q
    squarings = 0
    while z > F(1, 2):
        z /= 2
        squarings += 1
    # Alternating series on [0,1/2]: S_39 <= exp(-z) <= S_40.
    term = total = F(1)
    lo = None
    for k in range(1, 41):
        term *= -z/k
        total += term
        if k == 39:
            lo = total
    hi = total
    assert 0 < lo <= hi <= 1
    lo, hi = down(lo), up(hi)
    for _ in range(squarings):
        lo, hi = down(lo*lo), up(hi*hi)
    return lo, hi


def atan_bounds(q, terms):
    q = F(q)
    terms_list = [(-1)**k*q**(2*k+1)/(2*k+1) for k in range(terms+1)]
    s = sum(terms_list[:-1], F(0))
    t = s+terms_list[-1]
    return min(s, t), max(s, t)


def pi_certificate():
    lo1, hi1 = atan_bounds(F(1, 5), 32)
    lo2, hi2 = atan_bounds(F(1, 239), 8)
    lo, hi = 16*lo1-4*hi2, 16*hi1-4*lo2
    lower = F(3141592653589793, 10**15)
    upper = F(1570796326794897, 5*10**14)
    assert lower < lo < hi < upper
    return lower, upper, {'machin_lower': str(lo), 'machin_upper': str(hi),
                         'pi_lower': str(lower), 'pi_upper': str(upper)}


def support(row):
    r, u, v = map(F, (row['r'], row['u'], row['v']))
    assert r >= 0
    ell = -u*(r-3)-v*(r*r-10*r+15)
    assert ell == F(row['ell'])
    penalty = F(3, 2)*u*u+30*v*v/LAMBDA
    assert penalty == F(row.get('penalty', 0))
    constant = F(row['constant_upper'])
    if u == v == 0:
        assert ell == constant == 0 and row['pieces'] == []
        return {'r': '0', 'u': '0', 'v': '0', 'ell': '0', 'penalty': '0',
                'constant_upper': '0', 'integral_lower': '0', 'integral_upper': '0',
                'constant_gap': '0', 'pieces': [], 'exponential_coefficients': []}
    aa, bb = ell-3*u+15*v, u-10*v
    # Q(y) = (y-r)(v*y + u + v*(r-10)); all its nonnegative roots are known.
    assert aa == -r*(u+v*(r-10)) and bb == u+v*(r-10)-v*r
    roots = {F(0)}
    if r > 0:
        roots.add(r)
    if v:
        r2 = 10-r-u/v
        if r2 > 0:
            roots.add(r2)
    cuts = sorted(roots)+[None]
    coefficients = {}
    pieces = []

    def add(y, coefficient):
        if y is not None:
            coefficients[y] = coefficients.get(y, F(0))+coefficient

    def multiplier(y):
        return aa+bb*(y+2)+v*(y*y+4*y+8)

    for lo, hi in zip(cuts, cuts[1:]):
        test = lo+1 if hi is None else (lo+hi)/2
        value = aa+bb*test+v*test*test
        assert value != 0
        sign = 1 if value > 0 else -1
        add(lo, sign*multiplier(lo))
        if hi is not None:
            add(hi, -sign*multiplier(hi))
        pieces.append({'lo': str(lo), 'hi': 'infinity' if hi is None else str(hi), 'sign': sign})
    assert pieces == [{k: p[k] for k in ('lo', 'hi', 'sign')} for p in row['pieces']]
    ilo = ihi = F(0)
    exponentials = []
    for y, coefficient in sorted(coefficients.items()):
        low, high = exp_negative(y/2)
        if coefficient >= 0:
            ilo += coefficient*low
            ihi += coefficient*high
        else:
            ilo += coefficient*high
            ihi += coefficient*low
        exponentials.append({'argument': str(y/2), 'coefficient': str(coefficient),
                             'exp_lower': str(low), 'exp_upper': str(high)})
    assert 0 <= ilo <= ihi
    assert penalty+ihi <= constant, ('invalid agreement constant', r, u, v, penalty+ihi, constant)
    return {'r': str(r), 'u': str(u), 'v': str(v), 'ell': str(ell),
            'penalty': str(penalty), 'constant_upper': str(constant),
            'integral_lower': str(ilo), 'integral_upper': str(ihi),
            'constant_gap': str(constant-penalty-ihi),
            'pieces': pieces, 'exponential_coefficients': exponentials}


def envelope(lines, lo, hi):
    """Sweep by the first future intersection, certify against every line."""
    x = lo
    intervals = []
    while x < hi:
        index = min(range(len(lines)), key=lambda i: (lines[i][0]*x+lines[i][1], lines[i][0], i))
        m, b = lines[index]
        next_x = hi
        for mm, bb in lines:
            if mm < m:
                z = (bb-b)/(m-mm)
                if z > x:
                    next_x = min(next_x, z)
        assert next_x > x
        # Every competing difference is affine: its two endpoint signs prove
        # global dominance on this entire closed interval.
        assert all(m*z+b <= mm*z+bb for z in (x, next_x) for mm, bb in lines)
        intervals.append({'lo': str(x), 'hi': str(next_x), 'index': index,
                          'slope': str(m), 'intercept': str(b)})
        x = next_x
    assert intervals[0]['lo'] == str(lo) and intervals[-1]['hi'] == str(hi)
    assert all(a['hi'] == b['lo'] for a, b in zip(intervals, intervals[1:]))
    return intervals


def check(source_dir, output_dir):
    started = time.monotonic()
    output_dir.mkdir(parents=True, exist_ok=True)
    source = read(source_dir/'agreement_quintic_closed_supports_01.json')
    rows = [op['supports'] for op in source['operators'] if op['lambda_'] == str(LAMBDA)]
    assert len(rows) == 1 and len(rows[0]) == 1007
    proved = [support(row) for row in rows[0]]
    agreement_lines = [(F(p['ell']), F(p['constant_upper'])) for p in proved]
    fibers = []
    for label in 'ABC':
        t = read(source_dir/f'quintic_176_target_{label}_05.json')
        assert F(t['a']) == A and F(t['lambda']) == LAMBDA
        assert F(t['D_bound']) == F(t['t'])+F(t['kappa'])
        fibers.append((-F(t['t']), F(t['D_bound'])))
    lo, hi = F(7, 25), F(1)
    agreement_envelope = envelope(agreement_lines, lo, hi)
    fiber_envelope = envelope(fibers, lo, hi)
    cuts = sorted({lo, hi} | {F(row[k]) for env in (agreement_envelope, fiber_envelope)
                             for row in env for k in ('lo', 'hi')})
    endpoints = []
    for x in cuts:
        av = min(m*x+b for m, b in agreement_lines)
        fv = min(m*x+b for m, b in fibers)
        value = A*x*x+av+fv
        assert value <= C
        endpoints.append({'x': str(x), 'value': str(value), 'gap': str(C-value)})
    old = read(source_dir/'quintic_176_exact_assembly_05.json')
    # Comparison occurs after independent construction and exact proof.
    assert [(r['x'], r['value']) for r in endpoints] == [(r['x'], r['value']) for r in old['assembly']['endpoints']]
    max_endpoint = max(F(r['value']) for r in endpoints)
    assert max_endpoint == F(old['assembly']['maximum'])
    pi_lo, pi_hi, pi_proof = pi_certificate()
    circle_upper = 3*pi_hi/40+F(119, 32)/pi_lo+F(97, 144)*lo*lo
    assert circle_upper < C
    # Complete finite-matrix lower estimate for n=100000, R=12, q=40.
    n, R, q = 100000, 12, 40
    epsilon_squared = F(n, n-2)*F(1, 2**(2*q))+4*n*F(3, 8)**72
    epsilon = F(1, 5*10**11)
    assert epsilon_squared < epsilon*epsilon
    vector_lower = (A+1)*F(n-1, n)**2-1-2*A*epsilon
    denominator_upper = 2*C/pi_lo
    ratio_lower = vector_lower/denominator_upper
    assert ratio_lower > F(44, 25) > F(7, 4)
    assert A/(2*C) == F(1375, 2454)
    data = {'a': str(A), 'lambda': str(LAMBDA), 'C': str(C),
            'agreement_supports': proved, 'agreement_envelope': agreement_envelope,
            'fiber_envelope': fiber_envelope, 'endpoints': endpoints,
            'maximum_endpoint': str(max_endpoint), 'minimum_endpoint_gap': str(C-max_endpoint),
            'pi_certificate': pi_proof,
            'circle': {'interval': ['0', str(lo)], 'upper': str(circle_upper), 'gap': str(C-circle_upper)},
            'finite_matrix': {'dimension': n, 'cube_radius': R, 'mesh_exponent': q,
                              'epsilon_squared_upper': str(epsilon_squared), 'epsilon_upper': str(epsilon),
                              'vector_lower': str(vector_lower), 'denominator_upper': str(denominator_upper),
                              'ratio_lower': str(ratio_lower), 'gap_above_44_over_25': str(ratio_lower-F(44, 25)),
                              'gap_above_7_over_4': str(ratio_lower-F(7, 4))}}
    save(output_dir/'SCALAR_MATHEMATICAL_DATA_01.json', data)
    names = ['agreement_quintic_closed_supports_01.json', 'quintic_176_exact_assembly_05.json'] + \
            [f'quintic_176_target_{label}_05.json' for label in 'ABC']
    report = {'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'inputs': {name: hashlib.sha256((source_dir/name).read_bytes()).hexdigest() for name in names},
              'supports_checked': len(proved), 'agreement_intervals': len(agreement_envelope),
              'fiber_intervals': len(fiber_envelope), 'endpoints_checked': len(endpoints),
              'all_support_constants_valid': True, 'exact_envelope_match': True,
              'fiber_universal_inequalities_assumed_only_in_final_assembly': True,
              'elapsed_seconds': time.monotonic()-started}
    save(output_dir/'SCALAR_AUDIT_01.json', report)
    print(json.dumps(report, indent=2))


def main():
    assert __debug__
    p = argparse.ArgumentParser()
    p.add_argument('--sources', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    check(args.sources, args.output)


if __name__ == '__main__':
    main()
