"""Independent group-law audit for frontier extraction SAT/UNSAT outcomes."""
import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, EllipticCurve
from factor_base_yield_screen import bit_span


def point_sets(E,basis,offset,exponents,r):
    F=E.base_field();base=set()
    for x in bit_span(basis,offset):
        try:lifts=E.lift_x(F.from_integer(x),all=True)
        except ValueError:lifts=[]
        base.update(P for P in lifts if r*P==E(0))
    sets=[]
    for exponent in exponents:
        current=set(base)
        for _ in range(exponent):current={E(P[0]**2,P[1]**2) for P in current}
        sets.append(current)
    return sets


def reachable(E,sets):
    sums={E(0)}
    for points in sets:sums={S+P for S in sums for P in points}
    return sums-{E(0)}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--runs',help='comma-separated run directories')
    args=parser.parse_args()
    if args.output.exists():parser.error('Choose a new output')
    runs=args.runs.split(',') if args.runs else ['frontier-extraction-20260909-a','frontier-replication-20260909-b','frontier-replication-20260909-c']
    cache={};reports=[]
    for run in runs:
        report=json.loads((args.source/run/'result.json').read_text());r=int(report['subgroup_order']);n=int(report.get('n',17))
        F=GF(2**n,name='a');E=EllipticCurve(F,[1,1,0,0,1])
        P=E(F.from_integer(report['P'][0]),F.from_integer(report['P'][1]));Q=E(F.from_integer(report['Q'][0]),F.from_integer(report['Q'][1]))
        arms=[]
        for arm in report['arms']:
            key=(n,tuple(arm['basis']),arm['offset'],tuple(arm['exponents']))
            if key not in cache:cache[key]=reachable(E,point_sets(E,arm['basis'],arm['offset'],arm['exponents'],r))
            reach=cache[key]
            if len(reach)!=arm['predicted_coverage']:raise RuntimeError(f'Coverage mismatch for {run} {arm["name"]}: {len(reach)} != {arm["predicted_coverage"]}')
            counts=dict(verified_in_set=0,unsat_outside_set=0,timeouts_in_set=0,timeouts_outside_set=0)
            checks=[]
            for trial in arm['runs']:
                target=trial['a']*P+trial['b']*Q;inside=target in reach;status=trial['status']
                if status=='verified_relation':
                    if not inside:raise RuntimeError('Certified relation outside group reachability set')
                    counts['verified_in_set']+=1
                elif status=='exhausted_polynomial_models':
                    if inside:raise RuntimeError(f'Solver UNSAT contradicts group decomposition: {run} {arm["name"]} target {trial["target_index"]}')
                    counts['unsat_outside_set']+=1
                elif status=='timeout':counts['timeouts_in_set' if inside else 'timeouts_outside_set']+=1
                else:raise RuntimeError(f'Unexpected status {status}')
                checks.append(dict(target_index=trial['target_index'],status=status,group_reachable=inside))
            arms.append(dict(name=arm['name'],independent_group_coverage=len(reach),counts=counts,checks=checks))
        reports.append(dict(run=run,arms=arms))
    result=dict(scope='independent group-law audit of bounded n=17 batches',method='actual Sage curve additions; no Semaev/CNF/log oracle',
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reports=reports)
    args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps([{ 'run':r['run'],'arms':[{a['name']:a['counts']} for a in r['arms']]} for r in reports]))


if __name__=='__main__':main()
