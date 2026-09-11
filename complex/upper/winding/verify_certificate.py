"""Exact package aggregation, optional full stored-geometry or replay validation.

None of these modes reexecutes Arb's polynomial or transcendental operations.
Run run_replay.py before --replay to obtain newly computed native enclosures.
"""
import argparse
from itertools import zip_longest
from pathlib import Path
import json
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parent/'sources'))
from exact_core import load,digest,require,binary_ball
from package_core import (HERE,SETTINGS,Aggregate,check_inputs,compact_interval,compact_records,
                          fingerprint,geometry_summary,safe_path,source_hashes,unoptimized)
from geometry_check import fresh_geometry
from check_algebra import check_polynomials,expand_modes

def check_manifest():
    manifest=load(HERE/'SHA256SUMS.json')
    require(manifest['schema']==1,'Unknown manifest schema')
    for name,wanted in manifest['files'].items():
        require(digest(safe_path(HERE,name))==wanted,'Package hash differs: '+name)
    return len(manifest['files'])

def validate_native(directory,row,payload):
    folder=directory/f"mode_{row['index']:04d}"
    completion=load(folder/'complete.json')
    require(completion['fingerprint']==fingerprint(payload),'Checkpoint input binding differs')
    attempt=safe_path(folder,completion['attempt'])
    record_path,geometry_path=attempt/'mode.json',attempt/'geometry.jsonl.gz'
    require(digest(record_path)==completion['mode_sha256'] and digest(geometry_path)==completion['geometry_sha256'],'Checkpoint record bytes differ')
    record=load(record_path)
    require(record['fingerprint']==fingerprint(payload) and all(record[k]==v for k,v in row.items()),'Native mode identity differs')
    return record,geometry_path

def main():
    unoptimized()
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--raw-records',type=Path,help='Optional complete supplied native records in the manifest layout')
    group.add_argument('--replay',type=Path,help='Output of run_replay.py; validates every retained mode')
    parser.add_argument('--output',type=Path,help='Write a new exact JSON result (refuses replacement)')
    args=parser.parse_args()
    started=time.monotonic();files=check_manifest()
    data,plan,formula,norms,amplitude,phase,omission=check_inputs()
    polynomial_terms=check_polynomials(formula)
    independent=expand_modes(data)
    require(len(independent)==2962 and sum(w for k,w in independent)==1,'Independent phase expansion')
    for row in plan['retained']:
        key,weight=independent[row['index']]
        require(row['key']==list(map(str,key)) and row['weight']==str(weight),'Independent retained mode expansion')
    from exact_core import modes
    require(independent==modes(data),'Independent complete mode expansion')
    aggregation=Aggregate(plan,phase,omission)
    raw_by_index={};payload=None
    if args.raw_records:
        raw=load(HERE/'data/raw_records_manifest.json')['modes']
        raw_by_index={row['index']:row for row in raw}
        require(len(raw_by_index)==len(raw)==2204 and set(raw_by_index)=={r['index'] for r in plan['retained']},'Raw manifest mode completeness')
    if args.replay:
        envelope=load(args.replay/'inputs.json');payload=envelope['payload']
        require(envelope['fingerprint']==fingerprint(payload),'Replay input fingerprint')
        require(payload['source_hashes']==source_hashes(),'Replay source bytes differ')
        require(payload['data']==data and payload['plan']==plan and payload['formula']==formula and payload['settings']==SETTINGS,'Replay mathematical input differs')
    for row,compact in zip_longest(plan['retained'],compact_records()):
        require(row is not None and compact is not None and row['index']==compact['index'],'Compact data mode completeness or order')
        if args.replay:
            record,geometry_path=validate_native(args.replay,row,payload)
            fresh_geometry(geometry_path,record,SETTINGS)
            aggregation.add(row,record,False)
        elif args.raw_records:
            spec=raw_by_index[row['index']]
            for kind in ('mode','geometry'):
                p=safe_path(args.raw_records,spec[kind]['path'])
                require(p.stat().st_size==spec[kind]['bytes'] and digest(p)==spec[kind]['sha256'],'Raw record bytes differ')
            record=load(safe_path(args.raw_records,spec['mode']['path']))
            require(all(record[k]==v for k,v in row.items()),'Raw mode schedule differs')
            require(len(record['coefficients'])==len(compact['coefficients'])==row['degree']+1,'Raw coefficient completeness')
            for native,packed in zip(record['coefficients'],compact['coefficients']):
                a,b=binary_ball(native);lo,hi=compact_interval(packed)
                require(lo<=a<=b<=hi,'Compact interval does not contain the native enclosure')
            for key,value in compact['geometry'].items():
                require(record['geometry'][key]==value,'Compact geometric summary differs')
            fresh_geometry(safe_path(args.raw_records,spec['geometry']['path']),record,SETTINGS)
            aggregation.add(row,compact,True)
        else:
            aggregation.add(row,compact,True)
    require(aggregation.coefficients==792700,'Incomplete retained coefficients')
    if not args.replay:require(aggregation.leaves==2305576,'Supplied leaf count differs')
    result=aggregation.result(plan['target_gamma'])
    result.update(status='PASS',scope=('All recorded rational arc implications and complete coverage, with exact aggregation; native operations were not reexecuted by this validator.' if args.replay or args.raw_records else 'Package integrity, Gaussian polynomial identities, complete phase modes, compact coefficient aggregation and rational tails. Arc rectangles and native operations were not checked by this command.'),manifest_files=files,polynomial_terms=polynomial_terms,complete_modes=2962,retained_modes=2204,omitted_modes=758,exact_norms={k:str(v) for k,v in norms.items()},phase_amplitude=str(amplitude),seconds=time.monotonic()-started)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
