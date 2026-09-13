"""Standard-library exact reweighting of six supplied signed coefficient groups.

The group enclosures and original omission bound are numerical premises.
This program does not rerun the primitive coefficients or circle certificates.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent
CLASSES = ('00', '10', '20', '11', '3', '4')
DEN = 1 << 320
BASE = Q(1, 100)
SELECTED = BASE + Q(1, 2**28)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read(name):
    return json.loads((HERE/name).read_text())


def ceiling_grid(x):
    return Q(-((-x.numerator*DEN)//x.denominator), DEN)


def weight_polynomials(x):
    return {'00':1-x*x/2+3*x**4/32, '10':x-3*x**3/8,
            '20':x*x-x**4/3, '11':x*x-x**4/4, '3':x**3, '4':x**4}


def ratios(e):
    require(Q(9, 1000) <= e <= Q(11, 1000), 'Outside the proved epsilon interval')
    values, old = weight_polynomials(e), weight_polynomials(BASE)
    require(min(values.values()) > 0 and min(old.values()) > 0, 'Weight support changes')
    return {c:values[c]/old[c] for c in CLASSES}


def classify(i, j):
    a, b = sorted((abs(i), abs(j)))
    if a+b >= 3:
        return str(a+b)
    return {(0, 0):'00', (0, 1):'10', (0, 2):'20', (1, 1):'11'}[a,b]


def new_weight(i, j, e):
    a, b = abs(i), abs(j)
    sign = (-1)**((a if i < 0 else 0)+(b if j < 0 else 0))
    total = Q()
    for p in range(3):
        for q in range(3):
            degree = a+b+2*p+2*q
            if degree <= 4:
                total += sign*(-1)**(p+q)*(e/2)**degree / (
                    factorial(p)*factorial(p+a)*factorial(q)*factorial(q+b))
    return total


def phase_bound(candidate, e):
    g = Q(candidate['damping_lower'])
    dw, dx = Q(candidate['damping_w']), Q(candidate['damping_x'])
    require(g > 0 and dw > 0 and dx > 0 and dw*dx >= g*g, 'Damping inequality failed')
    M = Q()
    for amplitude, frequency, damping in candidate['profiles']:
        a, l, u = Q(amplitude), abs(Q(frequency)), Q(damping)
        require(u >= 0, 'Negative damping')
        M += abs(a)*(1 if u == 0 else min(Q(1), 2*l/(5*u*g)))
    return (2*M)**5/factorial(5) + sum(((2*M)**j/factorial(j) for j in range(5)), Q())*(2*e)**5/factorial(5)


def evaluate(groups, candidate, omitted, e):
    rr = ratios(e)
    nmax = max(len(g['lower']) for g in groups.values())
    lower, upper = [Q() for _ in range(nmax)], [Q() for _ in range(nmax)]
    for c in CLASSES:
        for n, (lo, hi) in enumerate(zip(groups[c]['lower'], groups[c]['upper'])):
            lower[n] += rr[c]*Q(lo, DEN)
            upper[n] += rr[c]*Q(hi, DEN)
    linear = lower[0]
    nonlinear = sum((max(abs(lo), abs(hi)) for lo, hi in zip(lower[1:], upper[1:])), Q())
    tails = sum((rr[c]*Q(groups[c]['tail'], DEN) for c in CLASSES), Q())
    phase = ceiling_grid(phase_bound(candidate, e))
    omission = ceiling_grid(max(rr.values())*omitted)
    gamma = linear-nonlinear-tails-phase-omission
    return {'epsilon':str(e), 'ratios':{c:str(rr[c]) for c in CLASSES},
            'linear_lower':str(linear), 'nonlinear_upper':str(nonlinear),
            'full_scalar_tail_upper':str(tails), 'full_phase_error_upper':str(phase),
            'full_omission_cost_upper':str(omission), 'gamma_lower':str(gamma),
            'coefficient_interval_width_sum':str(sum((hi-lo for lo, hi in zip(lower, upper)), Q()))}


def verify():
    manifest = read('immutable_manifest.json')
    actual = {p.relative_to(HERE).as_posix() for p in HERE.rglob('*')
              if p.is_file() and p.name != 'immutable_manifest.json' and '__pycache__' not in p.parts}
    require(actual == set(manifest['files']), 'Missing or additional package files')
    for name, record in manifest['files'].items():
        path = HERE/name
        require(not path.is_symlink() and path.resolve().is_relative_to(HERE), 'File escapes package')
        raw = path.read_bytes()
        require(len(raw) == record['bytes'] and sha256(raw).hexdigest() == record['sha256'],
                'File identity differs: '+name)
    candidate = read('certificate/baseline_candidate.json')
    point = read('certificate/endpoint.json')
    require(sha256((HERE/'certificate/baseline_candidate.json').read_bytes()).hexdigest() == point['candidate_sha256'],
            'Candidate identity differs')
    require(candidate == read('source/candidate_definition.json'), 'Scientific source candidate differs')
    require(candidate['epsilon'] == '1/100' and candidate['mixing'] == '-1/50'
            and candidate['lambda_value'] == '1/4' and candidate['old_order'] == candidate['new_order'] == 4,
            'Different physical family')
    for name in ('P', 'B', 'C'):
        p = {int(k):Q(v) for k, v in candidate[name].items()}
        require(all(k > 0 and k % 2 for k in p) and sum(map(abs, p.values()), Q()) <= 1,
                'Inadmissible preprocessing')
    require(candidate['C'] == {'7':'1'}, 'Different third preprocessing')
    require(Q(point['epsilon']) == SELECTED and Q(point['base_epsilon']) == BASE, 'Different epsilon')
    rows = read('certificate/six_groups.json')
    require(rows['dyadic_bits'] == point['dyadic_bits'] == 320 and Q(rows['base_epsilon']) == BASE,
            'Different coefficient grid')
    groups = rows['groups']
    require(set(groups) == set(CLASSES), 'Incomplete group inventory')
    expected = {'00':(1011,2693,720103), '10':(685,5076,1188433),
                '20':(543,852,184729), '11':(567,1110,282020),
                '3':(414,210,39053), '4':(161,7,1080)}
    for c, g in groups.items():
        require((len(g['lower']), g['modes'], g['coefficients']) == expected[c], 'Different group coverage')
        require(len(g['lower']) == len(g['upper']), 'Mismatched coefficient arrays')
        require(all(type(x) is int for k in ('lower', 'upper') for x in g[k]), 'Noninteger endpoint')
        require(all(lo <= hi for lo, hi in zip(g['lower'], g['upper'])), 'Reversed interval')
        require(type(g['tail']) is int and g['tail'] >= 0, 'Invalid whole scalar tail')
    require(sum(g['modes'] for g in groups.values()) == point['retained_modes'] == 9948
            and sum(g['coefficients'] for g in groups.values()) == point['original_coefficients'] == 2415418
            and point['omitted_modes'] == 110522, 'Incomplete physical support')

    # Direct finite Taylor sums check every ordered new-phase frequency pair.
    checked = 0
    rr = ratios(SELECTED)
    for i in range(-4, 5):
        for j in range(-4, 5):
            if abs(i)+abs(j) <= 4:
                old = new_weight(i,j,BASE)
                require(old != 0 and new_weight(i,j,SELECTED) == rr[classify(i,j)]*old,
                        'Six-class identity failed')
                checked += 1
    require(checked == 41, 'Incomplete ordered phase support')
    omitted = Q(point['original_omission_cost'])
    require(omitted >= 0, 'Negative omission bound')
    baseline = evaluate(groups,candidate,omitted,BASE)
    selected = evaluate(groups,candidate,omitted,SELECTED)
    require(baseline == point['baseline_arithmetic'] and selected == point['selected_arithmetic'],
            'Signed objective or complete errors differ')
    gamma = Q(selected['gamma_lower'])
    require(gamma == Q(point['exact_gamma_lower']) > Q(baseline['gamma_lower']) > 0,
            'Incorrect or nonimproving coefficient objective')
    reciprocal = 1/gamma
    require(reciprocal == Q(point['exact_complex_upper']) < Q(point['strict_decimal_upper']) == Q('1.404698554831'),
            'Incorrect strict reciprocal endpoint')
    inventory = read('certificate/external_dependencies.json')
    cursor = 0
    for i, block in enumerate(inventory['coefficient_group_blocks']):
        require(block['path'] == f'blocks/block_{i:04d}.json'
                and block['indices'] == [cursor,min(cursor+128,9948)], 'Incomplete external block partition')
        cursor = block['indices'][1]
    require(cursor == 9948 and len(inventory['coefficient_group_blocks']) == 78, 'Missing external blocks')
    return {'status':'PASS_EXACT_SIX_GROUP_ENDPOINT', 'selected':selected,
            'exact_upper':str(reciprocal), 'strict_decimal_upper':'1.404698554831',
            'gain_in_gamma_over_base':str(gamma-Q(baseline['gamma_lower'])),
            'positive_endpoint_margin':str(Q('1.404698554831')-reciprocal),
            'group_interval_entries_checked':3381, 'ordered_weight_pairs_checked':41,
            'largest_radial_index':1010, 'external_block_contents_checked':False,
            'original_mode_coefficient_records_decoded':0, 'primitive_reexecutions':0,
            'scope':'Exact reweighting and endpoint consequences; original group, primitive, contour and omission executions remain premises.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output:
        require(not args.output.resolve().is_relative_to(HERE) and not args.output.exists(),
                'Output must be new and outside this package')
    result = verify()
    if args.output:
        with args.output.open('x') as f:
            json.dump(result, f, sort_keys=True, indent=2)
            f.write('\n')
    print('PASS: exact six-group consequence; K_G^C < 1.404698554831.')
    print('Original numerical enclosures and group construction remain premises; no physical calculation rerun.')


if __name__ == '__main__':
    main()
