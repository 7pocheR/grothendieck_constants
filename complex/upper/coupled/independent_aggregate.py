"""Independent integer reconstruction of complete, already validated mode sums.

Contour and primitive validity remains a dependency on the bound complete
validator execution. This program recomputes every coefficient contribution
without importing its coefficient decoder or accumulator.
"""
import argparse
import ctypes
import fcntl
from fractions import Fraction as Q
import hashlib
import json
from math import factorial
import os
from pathlib import Path
import resource
import signal
import sys
import time

sys.set_int_max_str_digits(50000)
BITS = 320
DENOMINATOR = 1 << BITS


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda: stream.read(1 << 20), b''):
            h.update(data)
    return h.hexdigest()


def object_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def load(path):
    def unique(items):
        out = {}
        for key, value in items:
            assert key not in out
            out[key] = value
        return out
    return json.loads(path.read_text(), object_pairs_hook=unique)


def bound(path, expected):
    assert digest(path) == expected, str(path)
    return load(path)


def rounded_ball(value, weight):
    assert set(value) == {'mid', 'rad'}
    m, me = value['mid']
    r, re = value['rad']
    assert all(type(x) is int for x in (m, me, r, re)) and r >= 0
    e = min(me, re)
    center, radius = m << (me-e), r << (re-e)
    first, second = sorted(((center-radius)*weight.numerator,
                            (center+radius)*weight.numerator))
    divisor = weight.denominator
    shift = e+BITS
    if shift >= 0:
        first, second = first << shift, second << shift
    else:
        divisor <<= -shift
    return first//divisor, -((-second)//divisor)


def up(value):
    value = Q(value)
    return -((-value.numerator*DENOMINATOR)//value.denominator)


def phase_error(candidate):
    lower = Q(candidate['damping_lower'])
    amplitude = Q(0)
    for eps, frequency, damping in candidate['profiles']:
        eps, frequency, damping = Q(eps), abs(Q(frequency)), Q(damping)
        assert damping >= 0
        factor = Q(1) if not damping else min(Q(1), 2*frequency/(5*damping*lower))
        amplitude += abs(eps)*factor
    old, new = candidate['old_order'], candidate['new_order']
    return (2*amplitude)**(old+1)/factorial(old+1) + sum(
        ((2*amplitude)**j/factorial(j) for j in range(old+1)), Q(0)
    )*(2*abs(Q(candidate['epsilon'])))**(new+1)/factorial(new+1)


ROOT = Path(__file__).resolve().parent


def write_atomic(path, value):
    import uuid
    assert not path.exists()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name+'.tmp.'+uuid.uuid4().hex)
    with temporary.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def package_files():
    manifest = load(ROOT/'SHA256SUMS.json')
    for name, expected in manifest.items():
        path = (ROOT/name).resolve()
        assert path.is_relative_to(ROOT) and digest(path) == expected
    return manifest


def relative(raw, name):
    assert not Path(name).is_absolute()
    path = (raw/name).resolve()
    assert path.is_relative_to(raw) and path.is_file()
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, required=True)
    parser.add_argument('--validation', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--stage', choices=('pilot', 'complete'), default='complete')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--cpu', type=int)
    args = parser.parse_args()
    if not (__debug__ and sys.flags.optimize == 0 and os.environ.get('PYTHONOPTIMIZE', '') in ('', '0')):
        raise ValueError('Run without -O and with PYTHONOPTIMIZE unset or zero')
    assert sys.platform == 'linux' and sys.dont_write_bytecode
    parent = os.getppid()
    assert ctypes.CDLL(None, use_errno=True).prctl(1, signal.SIGKILL) == 0
    assert os.getppid() == parent
    if args.cpu is not None:
        assert args.cpu in os.sched_getaffinity(0)
        os.sched_setaffinity(0, {args.cpu})
    seconds = 120 if args.stage == 'pilot' else 2400
    resource.setrlimit(resource.RLIMIT_AS, (2 << 30, 2 << 30))
    resource.setrlimit(resource.RLIMIT_CPU, (seconds+5, seconds+10))
    resource.setrlimit(resource.RLIMIT_FSIZE, (2 << 20, 2 << 20))
    signal.alarm(seconds)
    start = time.monotonic()
    package = package_files()
    source, raw, checked, output = ROOT/'sources', args.raw.resolve(), args.validation.resolve(), args.output.resolve()
    for other in (ROOT, raw, checked):
        assert not output.is_relative_to(other) and not other.is_relative_to(output)
    report_path = checked/'complete.json'
    expected_report = digest(report_path)
    report = bound(report_path, expected_report)
    binding = bound(checked/'binding.json', report['binding_sha256'])
    assert binding['package_sha256'] == object_hash(package)
    assert report['stage'] == 'complete' and report['validated'] == report['required'] == 9948
    assert not report['unresolved'] and report['unvisited_or_interrupted'] == 0
    assert report['status'] == 'COMPLETE_CRITERION_PASS' and report['full_candidate_certified'] is True
    manifest = bound(source/'source_manifest.json', binding['source_manifest_sha256'])
    for name, expected in manifest.items():
        assert digest(relative(source, name)) == expected
    assert binding['validator_manifest_sha256'] == digest(ROOT/'validation/manifest.json')
    assert binding['validator_source_sha256'] == digest(ROOT/'validation/checkpoint_validator.py')
    for name, expected in binding['evidence_dependency_sha256'].items():
        assert digest(relative(raw, name)) == expected
    inventory = load(raw/'complete_inventory.json')
    assert object_hash(inventory) == binding['complete_inventory_sha256']
    candidate = load(source/'candidate_definition.json')
    assert object_hash(candidate) == binding['candidate_sha256']
    rows = inventory['retained']
    assert len(rows) == len({r['label'] for r in rows}) == 9948
    assert {r['label'] for r in rows} == set(report['checkpoint_sha256'])
    assert {p.stem for p in (checked/'checkpoints').glob('mode_*.json')} == set(report['checkpoint_sha256'])
    selected = list(range(9948)) if args.stage == 'complete' else sorted({i*9947//31 for i in range(32)})
    execution = dict(package_sha256=object_hash(package), report_sha256=expected_report,
                     validation_binding_sha256=report['binding_sha256'], source_sha256=digest(Path(__file__)),
                     python=sys.version, platform=sys.platform, bits=BITS)
    fingerprint = object_hash(execution)
    output.mkdir(parents=True, exist_ok=True)
    with (output/'process.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        saved_binding = output/'binding.json'
        if saved_binding.exists():
            assert args.resume and load(saved_binding) == execution
        else:
            assert all(p.name == 'process.lock' or p.name.startswith('binding.json.tmp.') for p in output.iterdir())
            write_atomic(saved_binding, execution)
        lower, upper, tail, coefficients, times = [], [], 0, 0, []
        raw_bindings = hashlib.sha256()
        for ordinal in selected:
            started_mode = time.monotonic()
            row = rows[ordinal]
            checkpoint = bound(checked/'checkpoints'/(row['label']+'.json'), report['checkpoint_sha256'][row['label']])
            assert checkpoint['validator_source_sha256'] == binding['validator_source_sha256']
            descriptor = checkpoint['validated_descriptor']
            assert descriptor['row'] == row and descriptor['binding_sha256'] == report['binding_sha256']
            assert len(descriptor['dependency_sha256']) == 4
            for name, expected in sorted(descriptor['dependency_sha256'].items()):
                assert digest(relative(raw, name)) == expected
                raw_bindings.update(json.dumps([name, expected], separators=(',', ':')).encode()+b'\n')
            mode = load(raw/'replay'/row['label']/'mode.json')
            assert mode['label'] == row['label'] and mode['weight'] == row['weight'] and mode['key'] == row['key']
            weight = Q(row['weight'])
            assert len(mode['coefficients']) == mode['degree']+1 == descriptor['exact_summary']['coefficient_count']
            assert mode['geometry']['scalar_tail_upper'] == descriptor['exact_summary']['exact_unweighted_tail']
            dependency = dict(fingerprint=fingerprint, row=row,
                              validator_checkpoint_sha256=report['checkpoint_sha256'][row['label']],
                              raw_dependency_sha256=descriptor['dependency_sha256'])
            piece_path = output/'modes'/(row['label']+'.json')
            if piece_path.exists():
                saved = load(piece_path)
                assert saved['dependency'] == dependency
                piece = saved['contribution']
            else:
                pairs = [rounded_ball(value, weight) for value in mode['coefficients']]
                piece = dict(lower=[x[0] for x in pairs], upper=[x[1] for x in pairs],
                             tail=up(abs(weight)*Q(mode['geometry']['scalar_tail_upper'])), count=1, bits=BITS)
            assert object_hash(piece) == descriptor['exact_summary']['dyadic_contribution_sha256']
            assert piece['count'] == 1 and piece['bits'] == BITS
            assert len(piece['lower']) == len(piece['upper']) == len(mode['coefficients'])
            assert all(type(a) is int and type(b) is int and a <= b for a, b in zip(piece['lower'], piece['upper']))
            assert type(piece['tail']) is int and piece['tail'] >= 0
            # Bind both new and reused contributions to raw contents throughout this visit.
            for name, expected in descriptor['dependency_sha256'].items():
                assert digest(relative(raw, name)) == expected
            if not piece_path.exists():
                write_atomic(piece_path, dict(dependency=dependency, contribution=piece))
            while len(lower) < len(piece['lower']):
                lower.append(0)
                upper.append(0)
            for i, (a, b) in enumerate(zip(piece['lower'], piece['upper'])):
                lower[i] += a
                upper[i] += b
            tail += piece['tail']
            coefficients += len(piece['lower'])
            times.append(time.monotonic()-started_mode)
            if len(times) % 128 == 0:
                print(json.dumps({'modes': len(times), 'required': len(selected), 'seconds': time.monotonic()-start}), flush=True)
        result = dict(status='PASS_INDEPENDENT_INTEGER_AGGREGATION_PILOT', source_sha256=digest(Path(__file__)),
                      complete_report_sha256=expected_report, bound_source=binding,
                      selected_indices=selected, mode_count=len(selected), coefficient_count=coefficients,
                      raw_dependency_list_sha256=raw_bindings.hexdigest(), mode_seconds=times,
                      elapsed_seconds=time.monotonic()-start, independent_coefficient_decoder=True,
                      prior_complete_contour_verifier_execution_required=True,
                      new_physical_or_contour_evaluations=0, complete_aggregation_verified=False,
                      resume_dependency='A reused integer contribution is bound to this decoder source, its complete raw dependencies, and its exact validated contribution hash.')
        if args.stage == 'complete':
            assert {p.stem for p in (output/'modes').glob('mode_*.json')} == {r['label'] for r in rows}
            nonlinear = sum(max(abs(a), abs(b)) for a, b in zip(lower[1:], upper[1:]))
            phase, omission = up(phase_error(candidate)), up(Q(inventory['omission_cost']))
            terms = dict(linear_lower=lower[0], nonlinear_upper=nonlinear,
                         full_scalar_tail_upper=tail, full_phase_error_upper=phase,
                         full_omission_cost_upper=omission,
                         gamma_lower=lower[0]-nonlinear-tail-phase-omission)
            exact = {name: str(Q(value, DENOMINATOR)) for name, value in terms.items()}
            assert all(value == report['implication'][name] for name, value in exact.items())
            assert Q(exact['gamma_lower']) > Q(report['implication']['target'])
            result.update(status='PASS_COMPLETE_INDEPENDENT_INTEGER_AGGREGATION',
                          complete_aggregation_verified=True, exact_implication=exact,
                          rational_target=report['implication']['target'])
        assert digest(report_path) == expected_report and package_files() == package
        for name, expected in binding['evidence_dependency_sha256'].items():
            assert digest(relative(raw, name)) == expected
        destination = output/('complete.json' if args.stage == 'complete' else 'pilot.json')
        if destination.exists():
            previous = load(destination)
            assert previous['complete_report_sha256'] == expected_report and previous['source_sha256'] == result['source_sha256']
            assert previous.get('exact_implication') == result.get('exact_implication')
        else:
            write_atomic(destination, result)
        print(json.dumps({k:v for k,v in result.items() if k not in ('selected_indices', 'mode_seconds', 'bound_source')}, indent=2))


if __name__ == '__main__':
    main()
