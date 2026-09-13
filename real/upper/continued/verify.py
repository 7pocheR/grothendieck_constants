"""Exact endpoint verification from supplied numerical enclosures.

Uses only the Python standard library. It does not recompute the Gaussian
integrals, the Fourier pair products, or the localized derivative bounds.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PHYSICAL = '3077b5a8011868ff8a9200590bebb73cb41ad3ba10b09eb37860d5074c3b1653'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(name):
    return json.loads((HERE / name).read_text())


def objective(rows):
    require(len(rows) == 404, 'Expected all coefficients from 0 through 403')
    intervals = []
    for n, row in enumerate(rows):
        require(isinstance(row, list) and len(row) == 2, 'Malformed interval')
        lo, hi = map(Q, row)
        require(lo <= hi, 'Reversed interval')
        if n % 2 == 0:
            require(lo == hi == 0, 'Nonzero even coefficient')
        intervals.append((lo, hi))
    penalty = sum((max(abs(lo), abs(hi)) for lo, hi in intervals[3::2]), Q())
    return intervals[1][0] - penalty


def sqrt_sum_bounded(a, b, t):
    """Equivalent rational test for sqrt(a)+sqrt(b) <= t, a,b,t >= 0."""
    require(a >= 0 and b >= 0 and t >= 0, 'Negative square-root input')
    residual = t*t - a - b
    return residual >= 0 and residual*residual >= 4*a*b


def arctan_bounds(q, terms):
    require(q > 1 and terms > 0, 'Invalid alternating series')
    s = sum((Q((-1)**k, (2*k+1)*q**(2*k+1)) for k in range(terms)), Q())
    other = s + Q((-1)**terms, (2*terms+1)*q**(2*terms+1))
    return min(s, other), max(s, other)


def verify():
    manifest = read('immutable_manifest.json')
    actual = {f.relative_to(HERE).as_posix() for f in HERE.rglob('*')
              if f.is_file() and f.name != 'immutable_manifest.json'
              and '__pycache__' not in f.parts}
    require(actual == set(manifest['files']), 'Incomplete or additional package files')
    for name, expected in manifest['files'].items():
        path = HERE / name
        require(not path.is_symlink() and path.resolve().is_relative_to(HERE),
                'File escapes the package')
        raw = path.read_bytes()
        require(len(raw) == expected['bytes'] and sha256(raw).hexdigest() == expected['sha256'],
                'File identity mismatch: ' + name)

    raw = (HERE / 'certificate/physical.json').read_bytes()
    require(sha256(raw).hexdigest() == PHYSICAL, 'Different physical candidate')
    physical = json.loads(raw)
    p = {int(k): Q(v) for k, v in physical['P'].items()}
    require(set(p) == set(range(1, 82, 2)) and all(p.values()), 'Different primary support')
    require(sum(map(abs, p.values()), Q()) == Q(999999999999, 10**12), 'Primary norm differs')
    require(physical['dimension'] == 4 and Q(physical['collective_amplitude']) == 0,
            'Different sign function')
    require({int(k): Q(v) for k, v in physical['Q'].items()} == {1: Q(1)},
            'Different auxiliary preprocessing')
    require(Q(physical['frequency_step']) == Q(3, 4), 'Different profile frequencies')
    require([Q(x['frequency']) for x in physical['outer_profile']] == [2, 4, 6, 8],
            'Different outer frequencies')

    coeff = read('certificate/head_coefficients.json')
    bounds = read('certificate/scalar_bounds.json')
    plan = read('certificate/head_error_input.json')
    require(coeff['degree'] == bounds['degree'] == 403, 'Different degree')
    require(bounds['physical_sha256'] == PHYSICAL, 'Different bound candidate')
    require(bounds['tail_panels'] == 2048 and bounds['cutoff_indices'] == [592, 1248],
            'Different circle cutoff')
    require(plan['parameters']['scalar_degree'] == 403
            and plan['parameters']['phase_degree'] == 31
            and len(plan['rows']) == 1 and plan['rows'][0]['frequency_cutoff'] == 768,
            'Different finite approximation')
    choices = plan['rows'][0]['frequency_choices']
    require([r['j'] for r in choices] == list(range(404)), 'Incomplete frequency allowance inventory')
    require(all(0 < Q(r['t']) < 1 and Q(r['y']) > 0 and Q(r['error_upper']) >= 0
                for r in choices), 'Invalid frequency allowance parameters')

    error = Q(bounds['head_analytic_error_upper'])
    require(error > 0 and error == Q(plan['rows'][0]['scalar_coefficient_error_upper']),
            'Different analytic error')
    finite = objective(coeff['finite_coefficients'])
    head = finite - error
    require(head == Q(bounds['head_lower_after_analytic_error']), 'Incorrect head arithmetic')
    a = Q(bounds['localized_squared_norm4_upper']) / (14*403**7)
    b = Q(bounds['localized_squared_norm6_upper']) / (22*403**11)
    tail = Q(bounds['tail_upper'])
    require(sqrt_sum_bounded(a, b, tail), 'Incorrect complete tail square-root bound')
    gamma = head - tail
    gamma0 = Q(bounds['gamma_strict_lower'])
    require(gamma0 == Q(17651832366373, 20000000000000), 'Different stated objective')
    require(gamma > gamma0 > 0, 'Objective margin is not positive')

    lower5, upper5 = arctan_bounds(5, 60)
    lower239, upper239 = arctan_bounds(239, 16)
    pi_lower, pi_upper = 16*lower5-4*upper239, 16*upper5-4*lower239
    endpoint = Q(bounds['K_strict_upper'])
    require(endpoint == Q(bounds['K_decimal_strict_upper']) == Q('1.779754412112'),
            'Different stated endpoint')
    require(pi_upper < 2*gamma0*endpoint, 'Strict endpoint comparison failed')

    dependency = read('certificate/external_dependencies.json')
    require(set(dependency['head_primitives']) ==
            {f'{kind}/node_{i:04d}.json.gz' for kind in ('features', 'profiles') for i in range(1, 769)},
            'Incomplete external primitive inventory')
    require(set(dependency['head_rows']) == {f'rows/row_{i:04d}.json.gz' for i in range(1, 769)},
            'Incomplete external row inventory')
    require(set(dependency['localized_tail_panels']) == {f'panel_{i:04d}.json' for i in range(2048)},
            'Incomplete external panel inventory')
    return {'status': 'PASS_EXACT_ENDPOINT_CONSEQUENCES',
            'scope': 'The supplied numerical enclosures are premises; no primitive, pair, or derivative-panel reexecution.',
            'package_files_checked': len(manifest['files']), 'coefficient_intervals_checked': 404,
            'finite_objective_lower': str(finite), 'analytic_error_subtracted_once': str(error),
            'head_lower': str(head), 'tail_upper': str(tail),
            'complete_objective_lower': str(gamma), 'strict_rational_objective_lower': str(gamma0),
            'positive_objective_margin': str(gamma-gamma0),
            'pi_enclosure': [str(pi_lower), str(pi_upper)],
            'strict_K_upper': str(endpoint), 'decimal_strict_K_upper': '1.779754412112',
            'positive_endpoint_margin': str(2*gamma0*endpoint-pi_upper),
            'external_data_contents_checked': False, 'primitive_integrals_reexecuted': 0,
            'pair_products_reexecuted': 0, 'localized_panels_reexecuted': 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output:
        require(not args.output.resolve().is_relative_to(HERE), 'Output must be outside the package')
        require(not args.output.exists(), 'Output already exists')
    result = verify()
    if args.output:
        with args.output.open('x') as f:
            json.dump(result, f, indent=2, sort_keys=True)
            f.write('\n')
    print('PASS: exact endpoint consequences; K_G^R < 1.779754412112.')
    print('Numerical enclosure premises retained; no primitive integrals, pair products or circle panels rerun.')


if __name__ == '__main__':
    main()
