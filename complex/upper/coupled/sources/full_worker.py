"""Bounded production and fresh outward replay of complete retained modes."""
from fractions import Fraction as Q
import gzip
import json
from pathlib import Path
import signal
import time

from coupled_exact import candidate, primitive, encode_primitives, decode_primitives
from coupled_arb import (formal_coefficients, geometry, suggested_degree, map_polynomials,
                         compiled, evaluate, map_arc, point, rectangle)
from full_common import atomic_json, digest
from winding_exact_core import load, require, binary_ball, arc_witness, dyad
from validate_records import validate_mode
from flint import ctx


def finite_width(coefficients, weight):
    return abs(Q(weight)) * sum(hi-lo for lo, hi in map(binary_ball, coefficients))


def records_overlap(a, b):
    for x, y in zip(a, b):
        lo, hi = binary_ball(x)
        left, right = binary_ball(y)
        require(max(lo, left) <= min(hi, right), 'Independent coefficient enclosures are disjoint')


def produce(row, d, config, binding, directory):
    directory.mkdir(parents=True, exist_ok=True)
    record = dict(row, source_binding=binding, status='UNRESOLVED', radius_attempts=[],
                  arithmetic_attempts=[])
    polynomials = primitive(d, row['key'])
    record['primitives'] = encode_primitives(polynomials)
    for radius in config['radii']:
        ctx.prec = config['precision_bits'][0]
        geometry_settings = dict(config, bits=ctx.prec, degree=0)
        name = 'geometry_' + radius.replace('/', '_') + '.jsonl.gz'
        try:
            with gzip.open(directory/name, 'wt', encoding='utf-8') as stream:
                proof = geometry(d, polynomials, Q(radius), 0, geometry_settings, stream)
            degree = suggested_degree(Q(row['weight']), Q(proof['boundary_modulus_upper']),
                                      Q(radius), Q(config['per_mode_tail']), cap=config['maximum_degree'])
            require(degree is not None, 'Required coefficient degree exceeds configured cap')
            degree = max(degree, config['minimum_degree'])
            proof['scalar_tail_upper'] = str(Q(proof['boundary_modulus_upper']) *
                Q(radius)**(-2*degree-3)/(1-Q(radius)**-2))
            for bits in config['precision_bits']:
                ctx.prec = bits
                coefficients = formal_coefficients(d, polynomials, degree)
                width = finite_width(coefficients, row['weight'])
                record['arithmetic_attempts'].append(dict(radius=radius, degree=degree,
                    bits=bits, weighted_width=str(width)))
                if width <= Q(config['per_mode_width']):
                    record.update(status='MODE_ARITHMETIC_PASS', degree=degree,
                        coefficients=coefficients, geometry=proof, geometry_file=name,
                        weighted_finite_coefficient_interval_width=str(width),
                        settings=dict(config, degree=degree, bits=bits))
                    atomic_json(directory/'mode.json', record)
                    # Independently check every rational implication immediately.
                    validate_mode(record, directory, d, binding)
                    return record
            raise ValueError('Coefficient width unresolved at precision cap')
        except (ValueError, ZeroDivisionError) as error:
            record['radius_attempts'].append(dict(radius=radius, reason=str(error)))
            atomic_json(directory/'unresolved.json', record)
    raise ValueError('No configured radius and arithmetic precision passed')


def fresh_geometry(record, source, output, d, polynomials, degree, settings):
    """Recompute all whole-arc enclosures on the production partition."""
    radius = Q(record['geometry']['radius'])
    maps, polys = map_polynomials(d, radius), compiled(polynomials)
    endpoints = {}
    def endpoint(angle):
        if angle == 2:
            angle = Q(0)
        if angle not in endpoints:
            z = point(angle, radius)
            endpoints[angle] = evaluate(polys, *(z*P(z*z) for P, _, __ in maps))
        return endpoints[angle]
    maximum, cursor, count = Q(0), Q(0), 0
    minima, sums = [None]*4, [[Q(0), Q(0)] for _ in range(3)]
    with gzip.open(source/record['geometry_file'], 'rt', encoding='utf-8') as incoming, \
         gzip.open(output/'geometry_fresh.jsonl.gz', 'wt', encoding='utf-8') as outgoing:
        for line in incoming:
            old = json.loads(line)
            left, right = Q(old['left']), Q(old['right'])
            require(left == cursor and left < right <= 2, 'Production contour gap or overlap')
            cursor = right
            values = evaluate(polys, *(map_arc(P, left, right, radius) for P in maps))
            z = values[3]*values[4]/(values[1]*values[2])
            first, last = endpoint(left), endpoint(right)
            ratios = [b/a for a, b in zip(first[:3], last[:3])]
            encoded = list(map(rectangle, values))
            encoded_z, encoded_ratios = rectangle(z), list(map(rectangle, ratios))
            witness = arc_witness(encoded, encoded_z, encoded_ratios)
            outgoing.write(json.dumps(dict(left=str(left), right=str(right), depth=old['depth'],
                values=encoded, Z=encoded_z, endpoint_ratios=encoded_ratios,
                exact_witness=witness), separators=(',', ':'))+'\n')
            count += 1
            maximum = max(maximum, dyad(witness['modulus_upper']))
            for i, pair in enumerate(witness['lower']):
                value = dyad(pair)
                minima[i] = value if minima[i] is None else min(minima[i], value)
            for i, (lo, hi) in enumerate(witness['argument_increments']):
                sums[i][0] += dyad(lo)
                sums[i][1] += dyad(hi)
    require(cursor == 2 and count == record['geometry']['leaves'], 'Incomplete production partition')
    require(all(-3 < lo <= hi < 3 for lo, hi in sums), 'Fresh zero winding unresolved')
    return dict(radius=str(radius), boundary_modulus_upper=str(maximum),
        scalar_tail_upper=str(maximum*radius**(-2*degree-3)/(1-radius**-2)),
        minima=list(map(str, minima)), argument_sums=[[str(lo), str(hi)] for lo, hi in sums],
        leaves=count, refinements=count-settings['panels'], evaluations=count,
        map_second_angular_derivative_bounds=[str(p[2]) for p in maps])


def replay(row, d, config, binding, source, output):
    output.mkdir(parents=True, exist_ok=True)
    original = load(source/'mode.json')
    for name in ('label', 'key', 'weight', 'omission_cost'):
        require(original[name] == row[name], 'Production mode identity differs')
    validate_mode(original, source, d, binding)
    # Reconstruct the exact determinant from the canonical key, not stored polynomials.
    polynomials = primitive(d, row['key'])
    require(polynomials == decode_primitives(original['primitives']), 'Rebuilt determinant differs')
    degree = original['degree'] + config['replay_extra_degree']
    bits = original['settings']['bits'] + config['replay_extra_bits']
    ctx.prec = bits
    settings = dict(config, degree=degree, bits=bits)
    proof = fresh_geometry(original, source, output, d, polynomials, degree, settings)
    coefficients = formal_coefficients(d, polynomials, degree)
    records_overlap(original['coefficients'], coefficients)
    width = finite_width(coefficients, row['weight'])
    require(width <= Q(config['per_mode_width']), 'Fresh coefficient width unresolved')
    require(abs(Q(row['weight']))*Q(proof['scalar_tail_upper']) <= Q(config['per_mode_tail']),
            'Fresh scalar tail exceeds configured budget')
    record = dict(row, status='MODE_ARITHMETIC_PASS', source_binding=binding,
        degree=degree, settings=settings, primitives=encode_primitives(polynomials),
        geometry=proof, geometry_file='geometry_fresh.jsonl.gz', coefficients=coefficients,
        weighted_finite_coefficient_interval_width=str(width),
        production_mode_sha256=digest(source/'mode.json'),
        production_geometry_sha256=digest(source/original['geometry_file']),
        qualification='Fresh outward coefficients and full-arc enclosures; same immutable analytic primitives.')
    atomic_json(output/'mode.json', record)
    validate_mode(record, output, d, binding)
    return record


def valid_resume(directory, row, d, binding, config, phase):
    path = directory/'mode.json'
    if not path.exists():
        return None
    record = load(path)
    require(record['source_binding'] == binding, 'Cannot reuse a record from different sources')
    for name in ('label', 'key', 'weight', 'omission_cost'):
        require(record[name] == row[name], 'Resumed mode identity differs')
    for name, value in config.items():
        require(record['settings'][name] == value, 'Resumed mode configuration differs')
    # Do not mistake a stale pass label for a verified record.
    validate_mode(record, directory, d, binding)
    require(finite_width(record['coefficients'], row['weight']) <= Q(config['per_mode_width']),
            'Resumed finite width exceeds budget')
    require(abs(Q(row['weight']))*Q(record['geometry']['scalar_tail_upper']) <= Q(config['per_mode_tail']),
            'Resumed scalar tail exceeds budget')
    return record


def work(task):
    row, config, binding, output, phase = task
    import resource
    resource.setrlimit(resource.RLIMIT_AS, (config['worker_memory_bytes'], config['worker_memory_bytes']))
    resource.setrlimit(resource.RLIMIT_FSIZE, (config['maximum_file_bytes'], config['maximum_file_bytes']))
    output, d = Path(output), candidate()
    started = time.monotonic()
    def alarm(signum, frame):
        raise TimeoutError('Per-mode wall-clock budget reached')
    signal.signal(signal.SIGALRM, alarm)
    signal.alarm(config['mode_seconds'])
    directory = output/phase/row['label']
    try:
        record = valid_resume(directory, row, d, binding, config, phase)
        if record is None:
            if directory.exists():
                # Preserve interrupted/failed attempts under a distinct non-mode name.
                saved = directory.with_name('incomplete_' + row['label'] + '_' + str(time.time_ns()))
                directory.rename(saved)
            record = (produce(row, d, config, binding, directory) if phase == 'production' else
                      replay(row, d, config, binding, output/'production'/row['label'], directory))
        result = dict(label=row['label'], phase=phase, status='PASS', degree=record['degree'],
            bits=record['settings']['bits'], radius=record['geometry']['radius'],
            leaves=record['geometry']['leaves'], mode_sha256=digest(directory/'mode.json'),
            geometry_sha256=digest(directory/record['geometry_file']))
    except Exception as error:
        result = dict(label=row['label'], phase=phase, status='UNRESOLVED',
                      error=type(error).__name__+': '+str(error))
    finally:
        signal.alarm(0)
    result['elapsed_seconds'] = time.monotonic()-started
    atomic_json(directory/'worker_result.json', result)
    return result
