"""Independent exact coverage and Gaussian polynomial integral verification.

No producer module is imported. Range bounds use interval Horner evaluation;
signed integrals use Gaussian monomial recurrences and complete analytic tails.
All assertions are mandatory: do not run Python with -O.
"""
from fractions import Fraction as F
from pathlib import Path
from functools import lru_cache
from itertools import combinations, product
import argparse
import hashlib
import heapq
import json
import math
import time

from flint import arb, ctx

HERE = Path(__file__).resolve().parent
PRECISION = 256
GRID_BITS = 52
TOLERANCE = F(1, 2**48)
ctx.prec = PRECISION
ctx.threads = 1


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()


def digest(obj):
    return hashlib.sha256(canonical(obj)).hexdigest()


def save(path, obj):
    path = Path(path)
    tmp = path.with_name(path.name + '.new')
    tmp.write_text(json.dumps(obj, sort_keys=True, indent=2) + '\n')
    tmp.replace(path)


def A(q):
    q = F(q)
    return arb(q.numerator) / q.denominator


def lower(q):
    assert q.is_finite()
    return F(str(q.lower().fmpq()))


def upper(q):
    assert q.is_finite()
    return F(str(q.upper().fmpq()))


def upward(q, bits=GRID_BITS):
    r = upper(q) * 2**bits
    return F(-((-r.numerator) // r.denominator), 2**bits)


def downward(q, bits=GRID_BITS):
    r = lower(q) * 2**bits
    return F(r.numerator // r.denominator, 2**bits)


def det(M):
    if not M:
        return F(1)
    return sum(((-1)**j * M[0][j] * det([r[:j] + r[j+1:] for r in M[1:]])
                for j in range(len(M))), F(0))


def principal_minors(M, strict=False):
    n = len(M)
    assert all(len(r) == n for r in M)
    assert all(M[i][j] == M[j][i] for i in range(n) for j in range(n))
    values = []
    for size in range(1, n+1):
        for inds in combinations(range(n), size):
            d = det([[M[i][j] for j in inds] for i in inds])
            assert d > 0 if strict else d >= 0
            values.append({'indices': list(inds), 'determinant': str(d)})
    return values


def ldl(M):
    """Exact rational LDL^T; its positive square root gives Cholesky."""
    n = len(M)
    L = [[F(int(i == j)) for j in range(n)] for i in range(n)]
    D = []
    for j in range(n):
        D.append(M[j][j] - sum(L[j][k]**2 * D[k] for k in range(j)))
        assert D[j] > 0
        for i in range(j+1, n):
            L[i][j] = (M[i][j] - sum(L[i][k]*L[j][k]*D[k]
                                     for k in range(j))) / D[j]
    assert all(sum(L[i][k]*D[k]*L[j][k] for k in range(n)) == M[i][j]
               for i in range(n) for j in range(n))
    return L, D


def hermite_monomials():
    rows = [[F(1)], [F(0), F(1)]]
    for j in range(1, 5):
        row = [F(0)] + rows[j]
        for k, v in enumerate(rows[j-1]):
            row[k] -= j*v
        rows.append(row)
    return [[A(v)/arb(math.factorial(j)).sqrt() for v in row] + [arb(0)]*(5-j)
            for j, row in enumerate(rows)]


HERMITE = hermite_monomials()
NU = (arb(2)/arb.pi()).sqrt()
INV_SQRT_2PI = 1/(2*arb.pi()).sqrt()


@lru_cache(maxsize=131072)
def endpoint(s):
    """Six antiderivatives of s^j phi(s), j=0,...,5."""
    x = A(s)
    phi = (-x*x/2).exp()*INV_SQRT_2PI
    J = [(1 + (x/arb(2).sqrt()).erf())/2, -phi]
    for j in range(2, 6):
        J.append(-x**(j-1)*phi + (j-1)*J[j-2])
    return tuple(J)


def horner(coeff, x):
    y = arb(0)
    for c in reversed(coeff):
        y = y*x + c
    return y


def cell(coeff, lo, hi):
    interval = arb(A((lo+hi)/2), A((hi-lo)/2))
    ran = horner(coeff, interval)
    jl, jh = endpoint(lo), endpoint(hi)
    mass = jh[0]-jl[0]
    if ran >= 0:
        return True, sum((c*(b-a) for c, a, b in zip(coeff, jl, jh)), arb(0))
    if ran <= 0:
        return True, -sum((c*(b-a) for c, a, b in zip(coeff, jl, jh)), arb(0))
    return False, abs(ran).upper()*mass.upper()


def tail_bound(coeff, R=12):
    """E[|p(S)| 1_{|S|>R}] <= sum |c_j| T_j(R), no truncation."""
    r = arb(R)
    exponential = (-r*r/2).exp()
    T = [(r/arb(2).sqrt()).erfc(), NU*exponential]
    for j in range(2, 6):
        T.append(NU*r**(j-1)*exponential + (j-1)*T[j-2])
    return sum((abs(c)*t for c, t in zip(coeff, T)), arb(0))


def absolute_integral(coeff, tolerance=TOLERANCE, retain=False):
    """Certified enclosure of E|p(S)|, with fail-closed subdivision limits."""
    assert len(coeff) == 6
    known = arb(0)
    unknown = arb(0)
    heap = []
    signed = []
    serial = 0
    splits = 0

    def add(lo, hi):
        nonlocal known, unknown, serial
        resolved, value = cell(coeff, lo, hi)
        if resolved:
            known += value
            if retain:
                signed.append({'lo': str(lo), 'hi': str(hi),
                               'lower': str(downward(value, 100)),
                               'upper': str(upward(value, 100))})
        else:
            bound = value.upper()
            assert bound >= 0
            unknown += bound
            # A floating priority changes only subdivision order, never an inequality.
            heapq.heappush(heap, (-float(bound), serial, lo, hi, bound))
            serial += 1

    for i in range(-12, 12):
        add(F(i), F(i+1))
    while heap and not (unknown <= A(tolerance)):
        assert splits < 20000, 'Unresolved integral; increase bounded precision/subdivision budget.'
        _, _, lo, hi, bound = heapq.heappop(heap)
        unknown -= bound
        mid = (lo+hi)/2
        assert mid not in (lo, hi)
        add(lo, mid)
        add(mid, hi)
        splits += 1
    tail = tail_bound(coeff)
    lo = known.lower()
    hi = (known + unknown + tail).upper()
    assert lo <= hi
    result = {'lower': str(downward(lo)), 'upper': str(upward(hi)),
              'splits': splits, 'signed_cells': 24 + 2*splits - serial,
              'uncertain_upper': str(upward(unknown, 100)),
              'tail_upper': str(upward(tail, 120))}
    if retain:
        unresolved = [{'lo': str(l), 'hi': str(h), 'upper': str(upward(b, 100))}
                      for _, _, l, h, b in heap]
        result['signed_panels'] = signed
        result['unresolved_panels'] = unresolved
        panels = sorted([(F(p['lo']), F(p['hi'])) for p in signed+unresolved])
        assert panels[0][0] == -12 and panels[-1][1] == 12
        assert all(a[1] == b[0] for a, b in zip(panels, panels[1:]))
    return result


def model(target, domains=False):
    a, lam = F(target['a']), F(target['lambda'])
    assert a == F(33, 20) and lam == F(3, 20)
    t, kappa, bound = map(F, (target['t'], target['kappa'], target['D_bound']))
    assert t > 0 and bound == t+kappa
    forcing = list(map(F, target['forcing_hermite']))
    assert forcing == [0, t, 0, 0, 0, 0]
    baseD = [F(1), lam, F(0)]
    baseE = [-a, F(1), lam]
    columns = [[arb(0) for _ in range(6)] for _ in range(6)]
    certificates = []
    assert len(target['matrices']) == 2
    for parity in range(2):
        M = [list(map(F, row)) for row in target['matrices'][parity]]
        cert = {'parity': parity, 'positive': principal_minors(M, True)}
        for label, base in [('D', baseD), ('E', baseE)]:
            delta = [[M[i][j]-(base[i] if i == j else 0) for j in range(3)]
                     for i in range(3)]
            cert[label] = principal_minors(delta)
        L, diag = ldl(M)
        cert['unit_lower'] = [[str(q) for q in row] for row in L]
        cert['diagonal'] = list(map(str, diag))
        for j in range(3):
            for degree in range(6):
                columns[2*j+parity][degree] = sum(
                    (A(L[i][j])*A(diag[j]).sqrt()*HERMITE[2*i+parity][degree]
                     for i in range(j, 3)), arb(0))
        certificates.append(cert)
    domain_proof = []
    box = [[F(lo), F(hi)] for lo, hi in target['domain']]
    assert len(box) == 6
    for j, (lo, hi) in enumerate(box):
        assert hi > 0 and lo == (0 if j < 2 else -hi)
        if domains:
            integral = absolute_integral(columns[j], retain=True)
            b = upward(A(F(integral['upper']))/NU)
            assert b <= hi, ('compact domain invalid', j, b, hi)
            domain_proof.append({'coordinate': j, 'bound': str(hi),
                                 'fresh_upper': str(b), 'gap': str(hi-b),
                                 'integral': integral})
    return columns, certificates, domain_proof


def vertex(target, columns, coords, retain=False):
    v = list(map(F, coords))
    p = [sum((2*A(q)*col[j] for q, col in zip(v, columns)), arb(0))
         for j in range(6)]
    t = A(target['t'])
    plus, minus = list(p), list(p)
    plus[1] += t
    minus[1] -= t
    integrals = [absolute_integral(c, retain=retain) for c in (plus, minus)]
    penalty = sum(q*q for q in v)
    high = upward(sum((A(F(i['upper'])) for i in integrals), arb(0))/(2*NU)-A(penalty))
    low = downward(sum((A(F(i['lower'])) for i in integrals), arb(0))/(2*NU)-A(penalty))
    return {'lower': str(low), 'upper': str(high), 'integrals': integrals}


def split(box):
    axis = max(range(6), key=lambda j: box[j][1]-box[j][0])
    lo, hi = box[axis]
    mid = (lo+hi)/2
    left, right = list(box), list(box)
    left[axis], right[axis] = (lo, mid), (mid, hi)
    return tuple(left), tuple(right)


def encode(v):
    return tuple(str(F(x)) for x in v)


def source_binding(tree, sources):
    bindings = read(HERE/'bindings_01.json')
    assert {str(p.relative_to(tree)) for p in tree.rglob('*') if p.is_file()} == set(bindings['tree_files'])
    for name, wanted in bindings['tree_files'].items():
        assert sha(tree/name) == wanted, ('tree file changed', name)
    for name, wanted in bindings['source_files'].items():
        assert sha(sources/name) == wanted, ('source file changed', name)
    return bindings


def reconstruct(tree, sources):
    """Rebuild every needed vertex and dyadic leaf from the closed root box."""
    bindings = source_binding(tree, sources)
    targets = read(HERE/'targets_01.json')
    paths = [f'{j:04b}' for j in range(16)]
    expected_sources = {n: bindings['source_files'][n] for n in
                        ('certify_transverse_full_04.py', 'gaussian_absolute_polynomial_04.py', 'dual_seed_library.json')}
    batch = read(tree/'batch_manifest.json')
    assert batch['sources'] == expected_sources
    assert batch['batch_sha256'] == bindings['source_files']['certify_transverse_full_batch_04.py']
    assert batch['tasks'] == [{'label': label, 'path': p} for label in 'ABC' for p in paths]
    assert batch['inputs'] == [{'path': f'quintic_176_target_{label}_05.json', 'label': label,
                                'sha256': bindings['source_files'][f'quintic_176_target_{label}_05.json']}
                               for label in 'ABC']
    task_vertices = []
    all_leaves = []
    summary = []
    for label in 'ABC':
        target = targets[label]
        full_input = read(sources/f'quintic_176_target_{label}_05.json')
        clean_input = {k: v for k, v in full_input.items() if k != 'seed_vertices'}
        assert clean_input == {k: v for k, v in target.items() if k != 'domain'}
        input_hash = bindings['source_files'][f'quintic_176_target_{label}_05.json']
        manifest = read(tree/label/'manifest.json')
        assert manifest['input'] == full_input and manifest['input_sha256'] == input_hash
        assert manifest['sources'] == expected_sources and manifest['paths'] == paths
        root_box = tuple(tuple(map(F, row)) for row in target['domain'])

        @lru_cache(maxsize=None)
        def box_at(path):
            if not path:
                return root_box
            assert path[-1] in '01'
            return split(box_at(path[:-1]))[int(path[-1])]

        needed = {}
        proposed_records = {}
        logged = {}
        target_leaves = []
        target_nodes = 0
        extra_log_records = 0
        for path in paths:
            folder = tree/label/path
            cp = read(folder/'checkpoint.json')
            assert cp['input'] == full_input
            assert cp['binding'] == {'sources': expected_sources, 'input_sha256': input_hash,
                                      'initial_path': path, 'domain': target['domain'],
                                      'integral_tolerance': '1e-10'}
            assert cp['stack'] == [], ('unfinished stack', label, path)
            assert [(row['coordinate'], F(row['bound'])) for row in cp['domain_certificate']] == \
                   [(i, row[1]) for i, row in enumerate(root_box)]
            leaf_paths = []
            with (folder/'leaves.jsonl').open() as stream:
                for line in stream:
                    row = json.loads(line)
                    name = row['path']
                    assert name.startswith(path) and set(name) <= {'0', '1'}
                    leaf_paths.append(name)
                    box = box_at(name)
                    correction = sum((hi-lo)**2 for lo, hi in box)/4
                    assert F(row['correction']) == correction
                    vs = list(product(*box))
                    assert len(vs) == 64 and len(row['vertex_upper']) == 64
                    bounds = list(map(F, row['vertex_upper']))
                    assert F(row['upper']) == max(bounds)+correction
                    assert F(row['upper']) <= F(target['D_bound'])
                    keys = []
                    for v, proposed in zip(vs, bounds):
                        key = encode(v)
                        if key not in needed:
                            needed[key] = proposed
                            proposed_records[key] = {proposed}
                        else:
                            needed[key] = min(needed[key], proposed)
                            proposed_records[key].add(proposed)
                        keys.append(key)
                    target_leaves.append({'label': label, 'path': name,
                                          'correction': str(correction), 'vertices': keys})
            ordered = sorted(leaf_paths)
            assert len(set(ordered)) == len(ordered), ('duplicate leaf', label, path)
            assert all(not b.startswith(a) for a, b in zip(ordered, ordered[1:]))
            assert sum((F(1, 2**(len(p)-4)) for p in ordered), F(0)) == 1
            assert cp['accepted'] == len(ordered)
            assert cp['nodes'] == 2*len(ordered)-1
            target_nodes += cp['nodes']
            with (folder/'vertices.jsonl').open() as stream:
                for line in stream:
                    row = json.loads(line)
                    key = encode(row['v'])
                    assert len(key) == 6
                    lo, hi = F(row['lower']), F(row['upper'])
                    assert lo <= hi
                    if key in logged:
                        logged[key].add(hi)
                    else:
                        logged[key] = {hi}
                    extra_log_records += 1
        # Every proposed leaf bound must occur in the bound vertex log. All
        # required vertices, including shared boundary vertices, are integrated.
        for key, bound in needed.items():
            assert key in logged and proposed_records[key] <= logged[key], ('missing vertex log', label, key)
        vertices = sorted(needed, key=lambda key: tuple(map(F, key)))
        ids = {key: len(task_vertices)+j for j, key in enumerate(vertices)}
        for key in vertices:
            task_vertices.append({'label': label, 'v': list(key), 'proposed_upper': str(needed[key])})
        for leaf in target_leaves:
            leaf['vertices'] = [ids[key] for key in leaf['vertices']]
            all_leaves.append(leaf)
        summary.append({'label': label, 'leaves': len(target_leaves), 'vertices': len(vertices),
                        'nodes': target_nodes, 'logged_vertices': len(logged),
                        'vertex_log_records': extra_log_records,
                        'unused_logged_vertices': len(set(logged)-set(needed))})
        box_at.cache_clear()
    result = {'targets': targets, 'vertices': task_vertices, 'leaves': all_leaves,
              'summary': summary, 'bindings_sha256': sha(HERE/'bindings_01.json'),
              'targets_sha256': sha(HERE/'targets_01.json'), 'checker_sha256': sha(__file__)}
    result['content_sha256'] = digest(result)
    return result


def check_inventory(inv):
    d = dict(inv)
    wanted = d.pop('content_sha256')
    assert digest(d) == wanted
    assert inv['checker_sha256'] == sha(__file__)
    assert inv['targets_sha256'] == sha(HERE/'targets_01.json')
    assert inv['bindings_sha256'] == sha(HERE/'bindings_01.json')
    assert inv['targets'] == read(HERE/'targets_01.json')


def replay(inv, out, part, parts, seconds, max_new):
    check_inventory(inv)
    assert 0 <= part < parts and parts == 4
    out.mkdir(parents=True, exist_ok=True)
    path = out/f'vertices_{part:02d}.jsonl'
    binding = {'inventory_sha256': inv['content_sha256'], 'checker_sha256': sha(__file__),
               'part': part, 'parts': parts, 'precision': PRECISION,
               'grid_bits': GRID_BITS, 'tolerance': str(TOLERANCE)}
    completed = {}
    if path.exists():
        with path.open() as stream:
            for line in stream:
                row = json.loads(line)
                assert row['binding'] == binding
                i = row['index']
                assert i % parts == part and i not in completed and 0 <= i < len(inv['vertices'])
                assert row['task_sha256'] == digest(inv['vertices'][i])
                completed[i] = row
    started = time.monotonic()
    models = {label: model(target)[0] for label, target in inv['targets'].items()}
    count = 0
    with path.open('a') as stream:
        for i in range(part, len(inv['vertices']), parts):
            if i in completed:
                continue
            if count >= max_new or time.monotonic()-started >= seconds:
                break
            task = inv['vertices'][i]
            label = task['label']
            cert = vertex(inv['targets'][label], models[label], task['v'])
            cert.update({'index': i, 'task_sha256': digest(task), 'binding': binding})
            cert['within_proposed_upper'] = F(cert['upper']) <= F(task['proposed_upper'])
            stream.write(json.dumps(cert, sort_keys=True, separators=(',', ':'))+'\n')
            stream.flush()
            completed[i] = cert
            count += 1
    expected = len(range(part, len(inv['vertices']), parts))
    status = {'binding': binding, 'completed': len(completed), 'expected': expected,
              'new_vertices': count, 'elapsed_seconds': time.monotonic()-started,
              'complete': len(completed) == expected, 'output_sha256': sha(path)}
    save(out/f'progress_{part:02d}.json', status)
    return status


def finish(inv, out, tree, sources):
    check_inventory(inv)
    rebuilt = reconstruct(tree, sources)
    assert rebuilt == inv, 'Inventory must reconstruct exactly from the original bound files.'
    domains = {}
    for label, target in inv['targets'].items():
        _, minors, compact = model(target, domains=True)
        domains[label] = {'matrices': minors, 'compact_domain': compact}
    returned = {}
    output_hashes = {}
    for part in range(4):
        path = out/f'vertices_{part:02d}.jsonl'
        output_hashes[path.name] = sha(path)
        for line in path.read_text().splitlines():
            row = json.loads(line)
            i = row['index']
            assert i not in returned and i % 4 == part and 0 <= i < len(inv['vertices'])
            binding = {'inventory_sha256': inv['content_sha256'], 'checker_sha256': sha(__file__),
                       'part': part, 'parts': 4, 'precision': PRECISION,
                       'grid_bits': GRID_BITS, 'tolerance': str(TOLERANCE)}
            assert row['binding'] == binding and row['task_sha256'] == digest(inv['vertices'][i])
            assert F(row['lower']) <= F(row['upper'])
            integrals = row['integrals']
            assert len(integrals) == 2
            for integral in integrals:
                assert F(integral['lower']) <= F(integral['upper'])
                assert 0 <= integral['splits'] <= 20000 and integral['signed_cells'] >= 0
                assert 0 <= F(integral['uncertain_upper']) <= TOLERANCE+F(1, 2**100)
                assert F(integral['tail_upper']) >= 0
            penalty = sum(F(q)**2 for q in inv['vertices'][i]['v'])
            recomputed_upper = upward(sum((A(F(r['upper'])) for r in integrals), arb(0))/(2*NU)-A(penalty))
            recomputed_lower = downward(sum((A(F(r['lower'])) for r in integrals), arb(0))/(2*NU)-A(penalty))
            assert F(row['upper']) == recomputed_upper and F(row['lower']) == recomputed_lower
            returned[i] = row
    assert set(returned) == set(range(len(inv['vertices'])))
    leaf_results = []
    for leaf in inv['leaves']:
        maximum = max(F(returned[i]['upper']) for i in leaf['vertices'])
        bound = maximum+F(leaf['correction'])
        target = F(inv['targets'][leaf['label']]['D_bound'])
        assert bound <= target, ('fresh leaf inequality fails', leaf['label'], leaf['path'], bound)
        leaf_results.append({'label': leaf['label'], 'path': leaf['path'],
                             'upper': str(bound), 'gap': str(target-bound)})
    summary = []
    for label in 'ABC':
        rows = [r for r in leaf_results if r['label'] == label]
        summary.append({'label': label, 'leaves': len(rows),
                        'minimum_gap': str(min(F(r['gap']) for r in rows)),
                        'maximum_upper': str(max(F(r['upper']) for r in rows))})
    result = {'complete': True, 'vertices': len(returned), 'leaves': len(leaf_results),
              'proposed_bounds_exceeded': sum(F(returned[i]['upper']) > F(t['proposed_upper'])
                                              for i, t in enumerate(inv['vertices'])),
              'summary': summary, 'domain_proofs': domains, 'leaf_results': leaf_results,
              'inventory_sha256': inv['content_sha256'], 'checker_sha256': sha(__file__),
              'output_hashes': output_hashes}
    save(out/'verified.json', result)
    return {k: v for k, v in result.items() if k not in ('leaf_results', 'domain_proofs')}


def main():
    assert __debug__, 'Assertions are required.'
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=('prepare', 'replay', 'finish', 'smoke'))
    ap.add_argument('--tree', type=Path)
    ap.add_argument('--sources', type=Path)
    ap.add_argument('--inventory', type=Path)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--part', type=int, default=0)
    ap.add_argument('--parts', type=int, default=4)
    ap.add_argument('--seconds', type=float, default=240)
    ap.add_argument('--max-new', type=int, default=4000)
    ap.add_argument('--memory-gib', type=float)
    args = ap.parse_args()
    if args.memory_gib is not None:
        import resource
        assert 1 <= args.memory_gib <= 8
        limit = int(args.memory_gib*1024**3)
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    if args.mode == 'prepare':
        assert not args.output.exists()
        inv = reconstruct(args.tree, args.sources)
        save(args.output, inv)
        print(json.dumps({'summary': inv['summary'], 'content_sha256': inv['content_sha256']}))
    elif args.mode == 'replay':
        print(json.dumps(replay(read(args.inventory), args.output, args.part, args.parts,
                                args.seconds, args.max_new)))
    elif args.mode == 'finish':
        print(json.dumps(finish(read(args.inventory), args.output, args.tree, args.sources)))
    else:
        result = {}
        targets = read(HERE/'targets_01.json')
        for label, target in targets.items():
            columns, matrices, compact = model(target, domains=True)
            z = vertex(target, columns, ['0']*6, retain=True)
            assert F(z['lower']) <= F(target['t']) <= F(z['upper'])
            corner = [row[1] for row in target['domain']]
            result[label] = {'matrices': matrices, 'compact_domain': compact,
                             'zero': z, 'corner': vertex(target, columns, corner, retain=True)}
        save(args.output, result)
        print(json.dumps({'smoke_completed': True, 'domain_columns': 18, 'vertices': 6,
                          'output_sha256': sha(args.output)}))


if __name__ == '__main__':
    main()
