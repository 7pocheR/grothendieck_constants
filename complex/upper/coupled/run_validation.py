"""Run the bounded complete contour validator, preserving resumable checkpoints."""
import argparse
import os
from pathlib import Path
import sys
from package_support import ROOT, SOURCE, disjoint, external, require, verify_package

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--stage', required=True, choices=('conditioning', 'complete'))
    parser.add_argument('--cpus', type=int, choices=(1, 2, 3, 4), default=4)
    parser.add_argument('--reference', action='store_true', help='Validate the original raw files against the supplied reference hashes')
    args = parser.parse_args()
    verify_package()
    require(sys.platform.startswith('linux'), 'The bounded exact validator requires Linux')
    evidence, output = args.evidence.resolve(), external(args.output)
    disjoint(evidence, output)
    command = [sys.executable, '-B', str(ROOT/'validation/driver.py'),
               '--frozen-source', str(SOURCE), '--evidence', str(evidence), '--output', str(output),
               '--stage', args.stage, '--cpus', str(args.cpus)]
    if args.reference:
        command.append('--reference')
    os.execv(sys.executable, command)
