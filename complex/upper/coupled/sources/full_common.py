"""Exact inventory, source binding, and dyadic coefficient accounting."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement
from math import factorial
from pathlib import Path
import json
import sys

# The explicitly bounded rational Cauchy tails can have over 8,000 decimal digits.
if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(50000)

from coupled_exact import HERE, candidate, complete_modes, input_facts, omission_cost
from winding_exact_core import load, require, modes as old_multinomial_modes, binary_ball


def digest(path):
    h = sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def source_binding():
    manifest = load(HERE / 'source_manifest.json')
    for name, value in manifest.items():
        path = (HERE / name).resolve()
        require(path.is_relative_to(HERE) and digest(path) == value,
                'Frozen source mismatch: ' + name)
    require(candidate() == load(HERE / 'candidate_definition.json'), 'Candidate differs')
    return digest(HERE / 'source_manifest.json')


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def independent_modes(d):
    """Multinomial expansion separate from the side-power recurrence."""
    old = old_multinomial_modes(dict(d, phase_order=d['old_order']))
    lam, epsilon = Q(d['lambda_value']), Q(d['epsilon'])
    generators = [(side, sign * lam, sign * epsilon / 2)
                  for side in (0, 1) for sign in (-1, 1)]
    new = defaultdict(Q)
    for n in range(d['new_order'] + 1):
        for choice in combinations_with_replacement(range(4), n):
            shifts, weight = [Q(0), Q(0)], Q(1)
            for i, count in Counter(choice).items():
                side, shift, coefficient = generators[i]
                shifts[side] -= count * shift
                weight *= coefficient**count / factorial(count)
            new[tuple(shifts)] += weight
    new = {k: w for k, w in new.items() if w}
    result = {}
    for (r, u, s, v), w in old:
        collected = defaultdict(Q)
        for (tau, sigma), z in new.items():
            first, second = sorted(((r, u, tau), (s, v, sigma)))
            collected[first + second] += w * z
        for key, weight in collected.items():
            if weight:
                require(key not in result, 'Independent inventory collided across old keys')
                result[key] = weight
    require(len(result) == 120470, 'Independent complete mode count differs')
    return result


def inventory(d, threshold, independently=False):
    stream = sorted(independent_modes(d).items()) if independently else list(complete_modes(d))
    # Sorting is only for durable identity, never separate sorting of side coordinates.
    stream.sort()
    require(len(stream) == 120470 and len({k for k, w in stream}) == 120470,
            'Missing or duplicate full mode')
    retained, omitted, h = [], Q(0), sha256()
    for index, (key, weight) in enumerate(stream):
        cost = omission_cost(d, key, weight)
        h.update(json.dumps([list(map(str, key)), str(weight)], separators=(',', ':')).encode() + b'\n')
        if cost <= threshold:
            omitted += cost
        else:
            retained.append(dict(label=f'mode_{index:06d}', key=list(map(str, key)),
                                 weight=str(weight), omission_cost=str(cost)))
    require(len(retained) == 9948 and threshold == Q(1, 10**10), 'Different retained set')
    return dict(mode_count=len(stream), retained=retained, omitted_count=len(stream)-len(retained),
                threshold=str(threshold), omission_cost=str(omitted), sorted_mode_digest=h.hexdigest())


def conditioning_rows(rows):
    """A fixed deterministic spectrum, including largest frequency and damping."""
    features = [lambda r: abs(Q(r['weight'])),
                lambda r: max(abs(Q(r['key'][i])) for i in (0, 3)),
                lambda r: sum(abs(Q(r['key'][i])) for i in (0, 3)),
                lambda r: sum(Q(r['key'][i]) for i in (1, 4)),
                lambda r: sum(abs(Q(r['key'][i])) for i in (2, 5)),
                lambda r: Q(r['omission_cost'])]
    selected = {}
    for feature in features:
        for r in sorted(rows, key=lambda r: (feature(r), r['label']), reverse=True)[:3]:
            selected[r['label']] = r
    # Preserve all retained representatives used in the previous pilot as well.
    from coupled_exact import selected_modes
    labels = ['base', 'new_one_side', 'new_opposite', 'new_same', 'damped_asymmetric',
              'damped_exchanged', 'frequency_two', 'frequency_four']
    by_key = {tuple(r['key']): r for r in rows}
    for r in selected_modes(candidate(), labels):
        if tuple(r['key']) in by_key:
            value = by_key[tuple(r['key'])]
            selected[value['label']] = value
    require(len(selected) <= 26, 'Conditioning spectrum exceeds its bound')
    return sorted(selected.values(), key=lambda r: r['label'])


def down(q, bits):
    q = Q(q)
    return (q.numerator << bits) // q.denominator


def up(q, bits):
    return -down(-Q(q), bits)


class DyadicAggregate:
    """Round each signed contribution outward before integer accumulation."""
    def __init__(self, bits=320):
        self.bits, self.lower, self.upper, self.tail = bits, [], [], 0
        self.count = 0

    def add(self, record):
        weight = Q(record['weight'])
        while len(self.lower) < len(record['coefficients']):
            self.lower.append(0)
            self.upper.append(0)
        for i, ball in enumerate(record['coefficients']):
            lo, hi = binary_ball(ball)
            a, b = sorted((weight * lo, weight * hi))
            self.lower[i] += down(a, self.bits)
            self.upper[i] += up(b, self.bits)
        self.tail += up(abs(weight) * Q(record['geometry']['scalar_tail_upper']), self.bits)
        self.count += 1

    def finish(self, inv, d, target):
        require(self.count == len(inv['retained']), 'Cannot accept incomplete retained set')
        nonlinear = sum(max(abs(a), abs(b)) for a, b in zip(self.lower[1:], self.upper[1:]))
        phase = up(Q(input_facts(d)['phase_error']), self.bits)
        omitted = up(Q(inv['omission_cost']), self.bits)
        gamma = self.lower[0] - nonlinear - self.tail - phase - omitted
        denominator = 1 << self.bits
        result = dict(total_modes=inv['mode_count'], retained_modes=self.count,
                      omitted_modes=inv['omitted_count'], dyadic_bits=self.bits,
                      linear_lower=str(Q(self.lower[0], denominator)),
                      nonlinear_upper=str(Q(nonlinear, denominator)),
                      full_scalar_tail_upper=str(Q(self.tail, denominator)),
                      full_phase_error_upper=str(Q(phase, denominator)),
                      full_omission_cost_upper=str(Q(omitted, denominator)),
                      gamma_lower=str(Q(gamma, denominator)), target=str(target),
                      sufficient_criterion_passed=Q(gamma, denominator) > target,
                      maximum_odd_index=len(self.lower)-1)
        return result
