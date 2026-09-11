"""Materialize the two exact reference input files if a raw copy lacks them."""
import argparse
import gzip
from hashlib import sha256
from pathlib import Path
import os
import uuid
from package_support import ROOT, external, load, require, verify_package

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, required=True)
    args = parser.parse_args()
    verify_package()
    raw = external(args.raw)
    binding = load(ROOT/'data/reference_binding.json')
    contents = {'complete_inventory.json':gzip.decompress((ROOT/'data/inventory.json.gz').read_bytes()),
                'execution_inputs.json':(ROOT/'data/reference_execution_inputs.json').read_bytes()}
    raw.mkdir(parents=True, exist_ok=True)
    for name, data in contents.items():
        require(sha256(data).hexdigest() == binding['evidence_dependency_sha256'][name], 'Reference input bytes differ')
        path = raw/name
        if path.exists():
            require(path.read_bytes() == data, 'Existing raw metadata differs; it will not be replaced')
        else:
            temporary = path.with_name(name+'.tmp.'+uuid.uuid4().hex)
            with temporary.open('xb') as stream:
                stream.write(data)
            os.link(temporary, path)
            temporary.unlink()
    print('PASS_REFERENCE_INPUT_BYTES: no coefficient or contour files were generated')
