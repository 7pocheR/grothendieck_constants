"""Exact geometry and supplied-result checks; no Gaussian integral replay.

Only the standard library is used. The supplied Gaussian endpoints are
outputs of the source-bound outward integration program, not axioms proved
by JSON aggregation. Use run_replay.py for fresh Gaussian containment checks.
"""
import os
import sys
if sys.flags.optimize or os.environ.get('PYTHONOPTIMIZE'):
    raise RuntimeError('Python optimization is unsupported; remove -O/-OO and PYTHONOPTIMIZE.')

import argparse
from fractions import Fraction as Q
from functools import lru_cache
import gzip
import hashlib
from itertools import combinations, product
import json
from math import lcm
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parent
GRID = 2**52
COUNTS = {'A': (54542, 11093), 'B': (90401, 17460), 'C': (107078, 18846)}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def digest(data):
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path, data):
    path = Path(path)
    resolved = path.resolve()
    manifest_path = ROOT/'certificate_manifest.json'
    protected = {manifest_path.resolve()}
    if manifest_path.exists():
        manifest = read(manifest_path)
        protected.update((ROOT/name).resolve() for name in manifest['files'])
        protected.update((ROOT/record['path']).resolve() for record in manifest['artifacts'].values())
    require(resolved not in protected, 'Output would overwrite a distributed file')
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name+'.new')
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')
    temporary.replace(path)


def package_path(relative):
    path = (ROOT/relative).resolve()
    require(path.is_relative_to(ROOT) and path != ROOT, 'Invalid package path')
    return path


def manifest_check():
    manifest = read(ROOT/'certificate_manifest.json')
    require(manifest['format_version'] == 1, 'Unsupported manifest')
    require((manifest['vertices'], manifest['leaves'], manifest['parts']) == (252021, 47399, 4), 'Wrong certificate counts')
    for name, wanted in manifest['files'].items():
        require(sha(package_path(name)) == wanted, 'File hash mismatch: '+name)
    require(manifest['integration_checker_sha256'] == sha(ROOT/'check_quintic_01.py'), 'Integration source mismatch')
    for name, record in manifest['artifacts'].items():
        path = package_path(record['path'])
        require(path.stat().st_size == record['bytes'] and sha(path) == record['sha256'], 'Compressed hash mismatch: '+name)
    return manifest


def unpack_artifact(manifest, name, destination):
    record = manifest['artifacts'][name]
    destination = Path(destination)
    if destination.exists():
        require(destination.stat().st_size == record['uncompressed_bytes'] and
                sha(destination) == record['uncompressed_sha256'], 'Existing artifact differs: '+name)
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name+'.new')
    h, total = hashlib.sha256(), 0
    with gzip.open(package_path(record['path']), 'rb') as source, temporary.open('wb') as output:
        for block in iter(lambda: source.read(1024**2), b''):
            total += len(block)
            require(total <= record['uncompressed_bytes'], 'Oversized artifact: '+name)
            h.update(block)
            output.write(block)
    require(total == record['uncompressed_bytes'] and h.hexdigest() == record['uncompressed_sha256'], 'Uncompressed hash mismatch: '+name)
    temporary.replace(destination)
    return destination


def load_inventory(manifest, path=None):
    if path is None:
        with gzip.open(package_path(manifest['artifacts']['inventory.json']['path']), 'rb') as stream:
            raw = stream.read()
        rec = manifest['artifacts']['inventory.json']
        require(len(raw) == rec['uncompressed_bytes'] and hashlib.sha256(raw).hexdigest() == rec['uncompressed_sha256'], 'Inventory byte mismatch')
        inv = json.loads(raw)
        del raw
    else:
        rec = manifest['artifacts']['inventory.json']
        require(sha(path) == rec['uncompressed_sha256'], 'Inventory byte mismatch')
        inv = read(path)
    wanted = inv.pop('content_sha256')
    require(digest(inv) == wanted == manifest['inventory_content_sha256'], 'Inventory content mismatch')
    inv['content_sha256'] = wanted
    require(inv['checker_sha256'] == manifest['integration_checker_sha256'], 'Inventory checker mismatch')
    require(inv['targets_sha256'] == sha(ROOT/'targets_01.json') and
            inv['bindings_sha256'] == sha(ROOT/'bindings_01.json'), 'Inventory input mismatch')
    require(inv['targets'] == read(ROOT/'targets_01.json'), 'Inventory target mismatch')
    require(len(inv['vertices']) == 252021 and len(inv['leaves']) == 47399, 'Incomplete inventory')
    return inv


def pi_bounds():
    def atan(inv, n):
        x = Q(1, inv)
        s = sum(((-1)**j*x**(2*j+1)/(2*j+1) for j in range(n)), Q(0))
        other = s+(-1)**n*x**(2*n+1)/(2*n+1)
        return min(s, other), max(s, other)
    a, b = atan(5, 64)
    c, d = atan(239, 16)
    # Machin's identity follows from this tangent identity and its angle branch.
    require((Q(120, 119)-Q(1, 239))/(1+Q(120, 119)*Q(1, 239)) == 1, 'Machin identity')
    lo, hi, den = 16*a-4*d, 16*b-4*c, 2**256
    return Q(lo.numerator*den//lo.denominator, den), -Q((-hi).numerator*den//hi.denominator, den)


def determinant(m):
    if not m:
        return Q(1)
    return sum(((-1)**j*m[0][j]*determinant([row[:j]+row[j+1:] for row in m[1:]])
                for j in range(len(m))), Q(0))


def matrix_checks(target):
    require(Q(target['a']) == Q(33, 20) and Q(target['lambda']) == Q(3, 20), 'Operator mismatch')
    require(Q(target['t']) > 0 and Q(target['D_bound']) == Q(target['t'])+Q(target['kappa']), 'Target bound mismatch')
    require(list(map(Q, target['forcing_hermite'])) == [0, Q(target['t']), 0, 0, 0, 0], 'Forcing mismatch')
    require(len(target['matrices']) == 2, 'Parity count')
    count = 0
    for raw in target['matrices']:
        m = [list(map(Q, row)) for row in raw]
        require(len(m) == 3 and all(len(row) == 3 for row in m), 'Matrix shape')
        require(all(m[i][j] == m[j][i] for i in range(3) for j in range(3)), 'Matrix symmetry')
        for name, diagonal in [('positive', [0, 0, 0]), ('D', [1, Q(3, 20), 0]),
                               ('E', [-Q(33, 20), 1, Q(3, 20)]), ('J', [Q(3, 20), 0, 0])]:
            z = [[m[i][j]-(diagonal[i] if i == j else 0) for j in range(3)] for i in range(3)]
            for n in range(1, 4):
                for ix in combinations(range(3), n):
                    value = determinant([[z[i][j] for j in ix] for i in ix])
                    require(value > 0 if name == 'positive' else value >= 0, 'Principal minor fails')
                    count += 1
    return count


def geometry(inv):
    """Reconstruct every exact corner directly; no producer-history premise."""
    require(set(inv['targets']) == set('ABC'), 'Target labels')
    targets, tasks = inv['targets'], inv['vertices']
    leaves_by_label = {label: [] for label in 'ABC'}
    for leaf in inv['leaves']:
        require(leaf['label'] in leaves_by_label, 'Unknown leaf target')
        p = leaf['path']
        require(isinstance(p, str) and 4 <= len(p) <= 256 and set(p) <= {'0', '1'}, 'Invalid leaf path')
        leaves_by_label[leaf['label']].append(leaf)
    depth = max(len(leaf['path']) for leaf in inv['leaves'])
    denominators = [Q(q).denominator for t in targets.values() for pair in t['domain'] for q in pair]
    scale = lcm(*denominators)*2**depth
    encoded, seen, counts = [], set(), {label: 0 for label in 'ABC'}
    for task in tasks:
        require(task['label'] in targets and len(task['v']) == 6, 'Invalid vertex')
        coords = []
        for raw in task['v']:
            q = Q(raw)*scale
            require(q.denominator == 1, 'Vertex not on exact subdivision grid')
            coords.append(q.numerator)
        key = (task['label'], *coords)
        require(key not in seen, 'Duplicate mathematical vertex')
        seen.add(key)
        encoded.append(tuple(coords))
        counts[task['label']] += 1
    del seen
    used = bytearray(len(tasks))
    summaries, minor_count, corner_count = [], 0, 0
    for label in 'ABC':
        target, leaves = targets[label], leaves_by_label[label]
        minor_count += matrix_checks(target)
        require((counts[label], len(leaves)) == COUNTS[label], 'Target count mismatch')
        root = tuple(tuple(int(Q(v)*scale) for v in pair) for pair in target['domain'])
        require(len(root) == 6 and all(hi > 0 and lo == (0 if j < 2 else -hi)
                    for j, (lo, hi) in enumerate(root)), 'Invalid complete root domain')
        paths = sorted(leaf['path'] for leaf in leaves)
        require(len(set(paths)) == len(paths), 'Duplicate leaf path')
        require(all(not b.startswith(a) for a, b in zip(paths, paths[1:])), 'Leaf prefix overlap')
        for prefix in (f'{i:04b}' for i in range(16)):
            require(sum((Q(1, 2**(len(p)-4)) for p in paths if p.startswith(prefix)), Q(0)) == 1, 'Incomplete initial subtree')

        @lru_cache(maxsize=None)
        def box(path):
            if not path:
                return root
            parent = list(box(path[:-1]))
            axis = max(range(6), key=lambda j: parent[j][1]-parent[j][0])
            lo, hi = parent[axis]
            require((lo+hi) % 2 == 0, 'Insufficient exact grid')
            mid = (lo+hi)//2
            parent[axis] = (lo, mid) if path[-1] == '0' else (mid, hi)
            return tuple(parent)

        for leaf in leaves:
            rectangle = box(leaf['path'])
            correction = Q(sum((hi-lo)**2 for lo, hi in rectangle), 4*scale*scale)
            require(Q(leaf['correction']) == correction, 'Incorrect interpolation correction')
            indices = leaf['vertices']
            require(len(indices) == 64 and len(set(indices)) == 64, 'Wrong corner multiplicity')
            for i, corner in zip(indices, product(*rectangle)):
                require(type(i) is int and 0 <= i < len(tasks), 'Invalid corner index')
                require(tasks[i]['label'] == label and encoded[i] == corner, 'Corner coordinate mismatch')
                used[i] = 1
                corner_count += 1
        box.cache_clear()
        summaries.append({'label': label, 'vertices': counts[label], 'leaves': len(leaves)})
    require(all(used), 'Unused or omitted inventory vertex')
    return {'targets': summaries, 'corners_checked': corner_count,
            'principal_minors_checked': minor_count, 'coverage': 'all closed boxes and shared boundaries'}


def integral_fields(data):
    lo, hi = Q(data['lower']), Q(data['upper'])
    require(lo <= hi, 'Reversed integral interval')
    require(type(data['splits']) is int and 0 <= data['splits'] <= 20000 and
            type(data['signed_cells']) is int and data['signed_cells'] >= 0, 'Integral counters')
    require(0 <= Q(data['uncertain_upper']) <= Q(1, 2**48)+Q(1, 2**100), 'Unresolved integral budget')
    require(Q(data['tail_upper']) >= 0, 'Negative tail bound')
    return lo, hi


def domains_check(inv, domains, pi_hi):
    require(set(domains) == set('ABC'), 'Incomplete domain certificates')
    count = 0
    for label in 'ABC':
        target = inv['targets'][label]
        columns = domains[label]['compact_domain']
        require(len(columns) == 6 and {c['coordinate'] for c in columns} == set(range(6)), 'Domain columns')
        for column in columns:
            bound = Q(target['domain'][column['coordinate']][1])
            high = Q(column['fresh_upper'])
            require(Q(column['bound']) == bound and 0 <= high <= bound and
                    Q(column['gap']) == bound-high, 'Compact domain bound')
            data = column['integral']
            lo, hi = integral_fields(data)
            require(2*high*high >= hi*hi*pi_hi, 'Domain normalization')
            signed, unresolved = data['signed_panels'], data['unresolved_panels']
            panels = sorted((Q(c['lo']), Q(c['hi'])) for c in signed+unresolved)
            require(panels and panels[0][0] == -12 and panels[-1][1] == 12 and
                    all(a < b for a, b in panels) and
                    all(a[1] == b[0] for a, b in zip(panels, panels[1:])), 'Domain panel coverage')
            require(len(signed) == data['signed_cells'] and len(panels) == 24+data['splits'], 'Domain panel counts')
            require(all(Q(c['lower']) <= Q(c['upper']) for c in signed) and
                    all(Q(c['upper']) >= 0 for c in unresolved), 'Domain panel ordering')
            # Outward rounding of individual recorded panels can add tiny
            # amounts. Require the coarser final endpoints to contain them.
            known_lo = sum((Q(c['lower']) for c in signed), Q(0))
            known_hi = sum((Q(c['upper']) for c in signed), Q(0))
            unknown = sum((Q(c['upper']) for c in unresolved), Q(0))
            require(lo <= known_lo and hi >= known_hi+unknown+Q(data['tail_upper']), 'Domain endpoint assembly')
            count += 1
    return count


def verify(manifest, inv, replay_dir=None):
    started = time.monotonic()
    geo = geometry(inv)
    if replay_dir is None:
        with gzip.open(package_path(manifest['artifacts']['verified.json']['path']), 'rb') as stream:
            raw = stream.read()
        require(hashlib.sha256(raw).hexdigest() == manifest['artifacts']['verified.json']['uncompressed_sha256'], 'Final report bytes')
        supplied = json.loads(raw)
        domains = supplied['domain_proofs']
    else:
        replay_dir = Path(replay_dir)
        supplied = None
        domain_file = read(replay_dir/'domains.json')
        require(domain_file['checker_sha256'] == inv['checker_sha256'] and
                domain_file['inventory_sha256'] == inv['content_sha256'], 'Fresh domain binding')
        domains = domain_file['domain_proofs']
    pi_lo, pi_hi = pi_bounds()
    domain_count = domains_check(inv, domains, pi_hi)
    highs, part_reports, exceeded = [None]*len(inv['vertices']), [], 0
    for part in range(4):
        name = f'vertices_{part:02d}.jsonl'
        path = package_path(manifest['artifacts'][name]['path']) if replay_dir is None else replay_dir/name
        opener = gzip.open if replay_dir is None else open
        binding = {'inventory_sha256': inv['content_sha256'], 'checker_sha256': inv['checker_sha256'],
                   'part': part, 'parts': 4, 'precision': 256, 'grid_bits': 52, 'tolerance': str(Q(1, 2**48))}
        h, count, size = hashlib.sha256(), 0, 0
        with opener(path, 'rb') as stream:
            for line in stream:
                h.update(line)
                size += len(line)
                row = json.loads(line)
                i = row['index']
                require(type(i) is int and 0 <= i < len(highs) and i % 4 == part and highs[i] is None, 'Duplicate/wrong vertex index')
                task = inv['vertices'][i]
                require(row['binding'] == binding and row['task_sha256'] == digest(task), 'Vertex source/task binding')
                require(len(row['integrals']) == 2, 'Wrong integral count')
                intervals = [integral_fields(z) for z in row['integrals']]
                ilow, ihigh = sum(z[0] for z in intervals), sum(z[1] for z in intervals)
                lo, hi = Q(row['lower']), Q(row['upper'])
                require(lo <= hi and GRID % hi.denominator == 0 and GRID % lo.denominator == 0, 'Vertex interval/grid')
                penalty = sum((Q(v)**2 for v in task['v']), Q(0))
                require(ilow >= 0 and hi+penalty >= 0 and
                        8*(hi+penalty)**2 >= ihigh*ihigh*pi_hi, 'Vertex upper normalization')
                require(lo+penalty <= 0 or 8*(lo+penalty)**2 <= ilow*ilow*pi_lo, 'Vertex lower normalization')
                goes_over = hi > Q(task['proposed_upper'])
                require(row['within_proposed_upper'] is (not goes_over), 'Proposed-bound comparison mismatch')
                exceeded += goes_over
                highs[i] = hi.numerator*(GRID//hi.denominator)
                count += 1
        require(count == len(range(part, 252021, 4)), 'Incomplete replay part')
        if replay_dir is None:
            rec = manifest['artifacts'][name]
            require(size == rec['uncompressed_bytes'] and h.hexdigest() == rec['uncompressed_sha256'], 'Uncompressed vertex hash')
        part_reports.append({'part': part, 'vertices': count, 'sha256': h.hexdigest()})
    require(all(v is not None for v in highs), 'Missing replay vertices')
    leaves, summary = [], []
    for leaf in inv['leaves']:
        upper = Q(max(highs[i] for i in leaf['vertices']), GRID)+Q(leaf['correction'])
        gap = Q(inv['targets'][leaf['label']]['D_bound'])-upper
        require(gap >= 0, 'Leaf target bound fails')
        leaves.append({'label': leaf['label'], 'path': leaf['path'], 'upper': str(upper), 'gap': str(gap)})
    for label in 'ABC':
        rows = [r for r in leaves if r['label'] == label]
        summary.append({'label': label, 'leaves': len(rows),
                        'minimum_gap': str(min(Q(r['gap']) for r in rows)),
                        'maximum_upper': str(max(Q(r['upper']) for r in rows))})
    result = {'complete': True, 'vertices': len(highs), 'leaves': len(leaves),
              'proposed_bounds_exceeded': exceeded, 'summary': summary, 'domain_proofs': domains,
              'leaf_results': leaves, 'inventory_sha256': inv['content_sha256'],
              'checker_sha256': inv['checker_sha256'],
              'output_hashes': {f'vertices_{r["part"]:02d}.jsonl': r['sha256'] for r in part_reports}}
    if supplied is not None:
        require(result == supplied, 'Reconstructed certificate differs from supplied complete report')
    report = {'status': 'COMPLETE_GEOMETRY_AND_STORED_DATA_PASS', 'complete': True,
              'vertices': len(highs), 'leaves': len(leaves), 'parts': part_reports, 'summary': summary,
              'geometry': geo, 'compact_domain_columns': domain_count, 'proposed_bounds_exceeded': exceeded,
              'checker_sha256': inv['checker_sha256'], 'inventory_sha256': inv['content_sha256'],
              'validator_sha256': sha(__file__), 'new_gaussian_integrals': 0,
              'trust': 'Gaussian endpoint containment depends on execution of the reviewed outward integration kernel. This check validates geometry, source bindings and rational assembly.',
              'elapsed_seconds': time.monotonic()-started}
    return report, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['integrity', 'geometry', 'stored'], default='stored')
    parser.add_argument('--output', type=Path, default=ROOT/'verification/stored_validation.json')
    args = parser.parse_args()
    manifest = manifest_check()
    if args.mode == 'integrity':
        report = {'status': 'COMPRESSED_FILES_AND_SOURCES_MATCH', 'mathematical_verification': False}
    else:
        inv = load_inventory(manifest)
        if args.mode == 'geometry':
            report = {'status': 'COMPLETE_EXACT_GEOMETRY_PASS', 'geometry': geometry(inv), 'new_gaussian_integrals': 0}
        else:
            report, _ = verify(manifest, inv)
    save(args.output, report)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
