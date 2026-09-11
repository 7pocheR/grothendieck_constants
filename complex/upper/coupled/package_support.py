"""Package integrity, exact metadata and external output-path checks."""
from pathlib import Path
from hashlib import sha256
import gzip
import json
import os
import sys
import uuid

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT/'sources'
sys.path.insert(0, str(SOURCE))
from coupled_exact import candidate, input_facts
from full_common import digest, source_binding
from winding_exact_core import load, require


def unoptimized():
    require(__debug__ and sys.flags.optimize == 0 and os.environ.get('PYTHONOPTIMIZE', '') in ('', '0'),
            'Run without -O and with PYTHONOPTIMIZE unset or zero')


def object_hash(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def safe_path(base, name):
    require(isinstance(name, str) and not Path(name).is_absolute(), 'Expected a relative file name')
    path = (base/name).resolve()
    require(path.is_relative_to(base.resolve()) and path.is_file(), 'File is missing or outside its directory')
    return path


def external(path):
    path = Path(path).resolve()
    require(not path.is_relative_to(ROOT) and not ROOT.is_relative_to(path),
            'Choose an explicit output directory outside the package')
    return path


def disjoint(a, b):
    require(not a.is_relative_to(b) and not b.is_relative_to(a), 'Output and raw evidence must be disjoint')


def compressed(name):
    with gzip.open(ROOT/'data'/name, 'rt', encoding='utf8') as stream:
        # Same duplicate-key policy as the ordinary JSON loader.
        def unique(items):
            result = {}
            for key, value in items:
                require(key not in result, 'Duplicate compressed JSON key')
                result[key] = value
            return result
        return json.load(stream, object_pairs_hook=unique)


def verify_package():
    unoptimized()
    manifest = load(ROOT/'SHA256SUMS.json')
    actual = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file()
              and '__pycache__' not in p.parts and p.name != 'SHA256SUMS.json'}
    require(actual == set(manifest), 'Missing or unexpected package files')
    for name, expected in manifest.items():
        require(digest(safe_path(ROOT, name)) == expected, 'Package hash differs: '+name)
    source_binding()
    return manifest


def code_binding():
    manifest = verify_package()
    return object_hash(manifest)


def atomic_new(path, value, replace=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    require(replace or not path.exists(), 'Refusing to replace completed evidence: '+str(path))
    temporary = path.with_name(path.name+'.tmp.'+uuid.uuid4().hex)
    with temporary.open('x', encoding='utf8') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def environment(native=False):
    import platform
    out = dict(python=sys.version, platform=platform.platform(), optimized=sys.flags.optimize)
    if native:
        import importlib.metadata
        import flint
        out.update(python_flint=importlib.metadata.version('python-flint'),
                   flint_version=str(getattr(flint, '__FLINT_VERSION__', 'not exposed')),
                   flint_release=str(getattr(flint, '__FLINT_RELEASE__', 'not exposed')))
        require(out['python_flint'] == '0.8.0', 'Use python-flint==0.8.0')
    return out
