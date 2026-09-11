"""Strict metadata and final rational comparison; no physical evaluations."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
from package_support import (ROOT, SOURCE, candidate, compressed, digest, input_facts,
                             load, object_hash, require, verify_package)

TERMS = ('linear_lower', 'nonlinear_upper', 'full_scalar_tail_upper',
         'full_phase_error_upper', 'full_omission_cost_upper', 'gamma_lower')
DEN = 1 << 320


def upward(value):
    value = Q(value)
    return Q(-((-value.numerator*DEN)//value.denominator), DEN)


def arithmetic(implication, endpoint):
    terms = {name:Q(implication[name]) for name in TERMS}
    require(all(DEN % q.denominator == 0 for q in terms.values()), 'A term is not on the 320-bit grid')
    require(all(terms[k] >= 0 for k in TERMS if k != 'gamma_lower'), 'Negative cost or linear lower bound')
    gamma = terms['linear_lower']-sum(terms[k] for k in TERMS[1:-1])
    require(gamma == terms['gamma_lower'], 'Exact subtraction differs')
    target, upper = Q(endpoint['gamma']), Q(endpoint['upper_bound'])
    require(target > 0 and upper == 1/target and gamma > target, 'Displayed bound lacks a strict rational margin')
    # The decimal below is only a checked outward display of the exact fraction.
    require(upper < Q('1.404698554869'), 'Decimal upper display is not outward')
    return gamma


def quick():
    verify_package()
    endpoint = load(ROOT/'data/endpoint.json')
    report_path, aggregate_path = ROOT/'data/reference_validation.json', ROOT/'data/reference_aggregation.json'
    require(digest(report_path) == endpoint['reference_validation_sha256'], 'Reference validation hash differs')
    require(digest(aggregate_path) == endpoint['reference_aggregation_sha256'], 'Reference aggregation hash differs')
    report, aggregate = load(report_path), load(aggregate_path)
    binding = load(ROOT/'data/reference_binding.json')
    require(digest(ROOT/'data/reference_binding.json') == report['binding_sha256'], 'Reference binding hash differs')
    require(aggregate['bound_source'] == binding and aggregate['complete_report_sha256'] == digest(report_path),
            'Reference aggregation binds a different complete validation')
    inv = compressed('inventory.json.gz')
    rows = inv['retained']
    labels = [r['label'] for r in rows]
    require(len(rows) == len(set(labels)) == 9948 and labels == sorted(labels), 'Retained inventory is incomplete or unordered')
    require(inv['mode_count'] == 120470 and inv['omitted_count'] == 110522 and inv['mode_count'] == len(rows)+inv['omitted_count'],
            'Complete inventory counts differ')
    require(object_hash(inv) == binding['complete_inventory_sha256'], 'Reference complete inventory hash differs')
    import gzip
    require(sha256(gzip.decompress((ROOT/'data/inventory.json.gz').read_bytes())).hexdigest() ==
            binding['evidence_dependency_sha256']['complete_inventory.json'], 'Reference inventory bytes differ')
    require(report['status'] == 'COMPLETE_CRITERION_PASS' and report['stage'] == 'complete' and
            report['validated'] == report['required'] == 9948 and report['full_candidate_certified'] is True and
            not report['unresolved'] and report['unvisited_or_interrupted'] == 0 and
            set(report['checkpoint_sha256']) == set(labels), 'Reference validation is partial')
    require(aggregate['status'] == 'PASS_COMPLETE_INDEPENDENT_INTEGER_AGGREGATION' and
            aggregate['complete_aggregation_verified'] is True and aggregate['selected_indices'] == list(range(9948)) and
            aggregate['mode_count'] == 9948 and aggregate['coefficient_count'] == 2415418,
            'Reference aggregation is partial')
    require(all(report['implication'][name] == aggregate['exact_implication'][name] for name in TERMS), 'Reference rational terms differ')
    d = candidate()
    facts = input_facts(d)
    require(object_hash(d) == binding['candidate_sha256'], 'Reference candidate differs')
    task = load(SOURCE/'full_task.json')
    reference_task = load(ROOT/'data/reference_task.json')
    require(task == {k:v for k,v in reference_task.items() if k != 'output'}, 'Mathematical task settings differ')
    require(digest(ROOT/'data/reference_task.json') == binding['task_sha256'], 'Reference task hash differs')
    require(digest(ROOT/'data/reference_execution_inputs.json') == binding['evidence_dependency_sha256']['execution_inputs.json'],
            'Reference outward execution metadata differs')
    threshold, g = Q(task['threshold']), Q(d['damping_lower'])
    require(threshold == Q(1,10**10) and inv['threshold'] == str(threshold), 'Omission threshold differs')
    for row in rows:
        key = tuple(map(Q, row['key']))
        require(len(key) == 6 and key[:3] <= key[3:] and key[1] >= 0 and key[4] >= 0, 'Noncanonical mode')
        cost = abs(Q(row['weight']))/((1+2*key[1]*g)*(1+2*key[4]*g))
        require(str(cost) == row['omission_cost'] and cost > threshold, 'Retained mode cost differs')
    implication = report['implication']
    require(Q(implication['full_phase_error_upper']) == upward(facts['phase_error']), 'Complete phase error differs')
    require(Q(implication['full_omission_cost_upper']) == upward(inv['omission_cost']), 'Complete omission summary differs')
    require(Q(implication['full_scalar_tail_upper']) <= 9948*Q(task['settings']['per_mode_tail'])+Q(9948,DEN),
            'Reported total tail exceeds the fixed budgets and outward rounding allowance')
    raw = compressed('reference_raw_sha256.json.gz')
    schedule = compressed('reference_schedule.json.gz')
    require([r['label'] for r in schedule] == labels and len(raw) == 4*9948, 'Reference raw inventory differs')
    require(sum(r['coefficient_count'] for r in schedule) == 2415418 and
            all(r['coefficient_count'] == r['degree']+1 for r in schedule), 'Reference coefficient counts differ')
    settings = task['settings']
    require(all(settings['minimum_degree'] <= r['degree']-settings['replay_extra_degree'] <= settings['maximum_degree'] and
                r['bits']-settings['replay_extra_bits'] in settings['precision_bits'] and
                r['radius'] in settings['radii'] for r in schedule), 'Reference schedule exceeds the unchanged settings')
    expected_raw = set()
    aggregate_hash = sha256()
    for row in schedule:
        label = row['label']
        names = [f'production/{label}/mode.json', f'production/{label}/geometry_'+row['radius'].replace('/','_')+'.jsonl.gz',
                 f'replay/{label}/mode.json', f'replay/{label}/geometry_fresh.jsonl.gz']
        expected_raw.update(names)
        for name in sorted(names):
            require(name in raw and len(raw[name]) == 64 and all(c in '0123456789abcdef' for c in raw[name]), 'Malformed raw hash')
            aggregate_hash.update(json.dumps([name,raw[name]],separators=(',',':')).encode()+b'\n')
    require(set(raw) == expected_raw and aggregate_hash.hexdigest() == aggregate['raw_dependency_list_sha256'],
            'Complete original raw dependency manifest differs')
    gamma = arithmetic(implication, endpoint)
    require(str(gamma) == endpoint['exact_gamma_lower'], 'Exact reference gamma differs')
    return dict(status='PASS_METADATA_AND_EXACT_ARITHMETIC', upper_bound=endpoint['upper_bound'],
                strict_gamma_margin=str(gamma-Q(endpoint['gamma'])), retained_modes=9948,
                complete_inventory_modes=120470, reference_coefficient_contributions=2415418,
                raw_files_evaluated=0, contour_checks_executed=0, complete_inventory_regenerated=False,
                scope='Package integrity and exact implications of reference metadata; outward execution and raw validation remain separate requirements.')


def final(validation, aggregation):
    package = verify_package()
    endpoint = load(ROOT/'data/endpoint.json')
    report_path = validation/'complete.json'
    report = load(report_path)
    binding = load(validation/'binding.json')
    aggregate = load(aggregation/'complete.json')
    aggregate_binding = load(aggregation/'binding.json')
    require(report['status'] == 'COMPLETE_CRITERION_PASS' and report['stage'] == 'complete' and
            report['validated'] == report['required'] == 9948 and not report['unresolved'] and
            report['unvisited_or_interrupted'] == 0 and report['full_candidate_certified'] is True,
            'Complete contour validation is required')
    require(aggregate['status'] == 'PASS_COMPLETE_INDEPENDENT_INTEGER_AGGREGATION' and
            aggregate['complete_aggregation_verified'] is True and aggregate['mode_count'] == 9948 and
            aggregate['selected_indices'] == list(range(9948)), 'Complete independent aggregation is required')
    require(aggregate['complete_report_sha256'] == digest(report_path) == aggregate_binding['report_sha256'],
            'Aggregation binds a different validation report')
    require(digest(validation/'binding.json') == report['binding_sha256'] and aggregate['bound_source'] == binding,
            'Validation binding differs')
    require(binding['package_sha256'] == aggregate_binding['package_sha256'] == object_hash(package) and
            aggregate['source_sha256'] == aggregate_binding['source_sha256'] == digest(ROOT/'independent_aggregate.py'),
            'Execution source package differs')
    require(set(report['checkpoint_sha256']) == {r['label'] for r in compressed('inventory.json.gz')['retained']},
            'Complete retained label set differs')
    require(all(report['implication'][k] == aggregate['exact_implication'][k] for k in TERMS), 'The six exact terms differ')
    gamma = arithmetic(aggregate['exact_implication'], endpoint)
    if binding['reference_records']:
        reference = load(ROOT/'data/reference_validation.json')['implication']
        require(all(reference[k] == aggregate['exact_implication'][k] for k in TERMS) and
                aggregate['coefficient_count'] == 2415418, 'Original raw reconstruction differs')
    return dict(status='PASS_COMPLETE_REPRODUCTION_ARITHMETIC', upper_bound=endpoint['upper_bound'],
                gamma_lower=str(gamma), strict_gamma_margin=str(gamma-Q(endpoint['gamma'])),
                coefficient_contributions=aggregate['coefficient_count'],
                equals_reference_exact_gamma=str(gamma) == endpoint['exact_gamma_lower'],
                scope='Final arithmetic of the bound complete executions; this command does not reevaluate raw enclosures.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--validation', type=Path)
    parser.add_argument('--aggregation', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    require(bool(args.validation) == bool(args.aggregation), 'Supply both execution directories or neither')
    result = final(args.validation, args.aggregation) if args.validation else quick()
    if args.output:
        from package_support import atomic_new, external
        atomic_new(external(args.output), result)
    print(json.dumps(result, indent=2))
