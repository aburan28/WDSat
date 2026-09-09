"""Continue validated n19 WDSat streams through relation dependency extraction."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import time

from sage.all import GF, EllipticCurve
from factor_base_frobenius_schedule import score_sets
from factor_base_sat_rank import RelationRank, columns, point_key
from factor_base_frobenius_sat import rewritten_coefficients, frobenius_power
from factor_base_yield_screen import bit_span, subgroup_generator
import wdsat_adapter


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,required=True);parser.add_argument('--solver',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--target-seed',type=int,default=202609104);parser.add_argument('--no-reuse-initial',action='store_true')
    parser.add_argument('--run-protocol',type=Path)
    args=parser.parse_args();args.solver=args.solver.resolve();config_path=args.source/'n19-wdsat-extraction.json';config=json.loads(config_path.read_text());args.output.mkdir(parents=True,exist_ok=False)
    for path in [Path(__file__),Path(wdsat_adapter.__file__),args.solver,config_path]:(args.output/path.name).write_bytes(path.read_bytes())
    if args.run_protocol is not None:(args.output/args.run_protocol.name).write_bytes(args.run_protocol.read_bytes())
    frontier=json.loads((args.source/'dimension-frontier-n19-20260909-a'/'result.json').read_text());admission=json.loads((args.source/'n19-sat-admission-20260909-a'/'result.json').read_text());initial=json.loads((args.source/'n19-wdsat-20260909-c'/'result.json').read_text())
    F=GF(2**19,name='a');E=EllipticCurve(F,[1,1,0,0,1]);r=int(frontier['subgroup_order']);P=subgroup_generator(E,r);secret=31415;Q=secret*P
    piP=frobenius_power(E,P,1);walk=E(0);lam=None;log_by_point={}
    for j in range(r):
        log_by_point[point_key(walk)]=j
        if walk==piP:lam=j
        walk+=P
    if lam is None:raise RuntimeError('No Frobenius eigenvalue')
    arms=[];states={}
    for name in ['trace-k7','sparse-k6']:
        candidate=next(row for row in frontier['candidates'] if row['name']==name);control=next(c for c in candidate['controls'] if c['name']=='best');exponents=control['schedule']
        base=[F.from_integer(v) for v in candidate['basis']];offset=F.from_integer(candidate['offset']);arm=dict(name=name+'-best',basis=candidate['basis'],offset=candidate['offset'],exponents=exponents,
            bases=[[v**(2**e) for v in base] for e in exponents],offsets=[offset**(2**e) for e in exponents],cofactor_multiplier=1,predicted_coverage=control['coverage'])
        col=columns(E,r,arm);reps=[]
        for x in bit_span(arm['basis'],arm['offset']):
            try:lifts=E.lift_x(F.from_integer(x),all=True)
            except ValueError:lifts=[]
            for point in lifts:
                if r*point==E(0):reps.append(log_by_point[point_key(point)])
        _,coverage=score_sets([[(pow(lam,e,r)*s)%r for s in reps] for e in exponents],r)
        if coverage!=control['coverage']:raise RuntimeError('Coverage reconstruction failed')
        arms.append(arm);states[arm['name']]=dict(columns=col,rank=RelationRank(r),runs=[],seconds=0.0,extraction=None,coverage_mask=None)
        sets=[[(pow(lam,e,r)*s)%r for s in reps] for e in exponents];states[arm['name']]['coverage_mask']=score_sets(sets,r)[0]
    rng=random.Random(args.target_seed+19);targets=[]
    while len(targets)<config['maximum_targets']:
        a,b=rng.randrange(r),rng.randrange(r);R=a*P+b*Q
        if R!=E(0):targets.append((a,b,R))
    # Reconstruct and reverify admitted Gaussian rows only for their original stream.
    if not args.no_reuse_initial:
        if args.target_seed!=202609104:raise ValueError('Initial records belong only to target seed 202609104')
        for arm in arms:
            state=states[arm['name']];prior=next(row for row in initial['arms'] if row['name']==arm['name'])
            admission_arm=next(row for row in admission['arms'] if row['name']==arm['name'])
            for trial in admission_arm['runs']:
                a,b,_=targets[trial['target_index']]
                if trial['a']!=a or trial['b']!=b:raise RuntimeError('Regenerated target stream disagrees with admission artifacts')
            for record in prior['runs']:
                if record['mode']!='gaussian':continue
                idx=record['target_index'];a,b,R=targets[idx];reachable=bool(state['coverage_mask'] & (1 << ((a+b*secret)%r)))
                if record['status']=='unsat':
                    if reachable:raise RuntimeError('Initial UNSAT contradicts reachability')
                elif record['status']=='verified_relation':
                    if not reachable:raise RuntimeError('Initial relation outside reachability')
                    coefficients=rewritten_coefficients(E,record['witness'],arm['exponents'],lam,r,state['columns']);lhs=E(0)
                    for key,pos in state['columns'].items():lhs+=coefficients[pos]*E(F.from_integer(key[0]),F.from_integer(key[1]))
                    if lhs!=R:raise RuntimeError('Initial relation rewrite failed')
                    d=state['rank'].add(coefficients,a,b)
                    if d is not None:raise RuntimeError('Unexpected extraction in initial rows')
                    record=dict(record);record.update(a=a,b=b,coefficients=coefficients,reused=True)
                state['runs'].append(record);state['seconds']+=record['seconds']
            if len(state['rank'].rows)!=prior['summaries']['gaussian']['rank']:raise RuntimeError('Initial rank reconstruction mismatch')
    (args.output/'targets.json').write_text(json.dumps([dict(index=i,a=a,b=b,target=point_key(R)) for i,(a,b,R) in enumerate(targets)],indent=2)+'\n')
    def report():
        return dict(scope='n19 Gaussian WDSat extraction; implementation evidence',target_seed=args.target_seed,
                    reused_initial=not args.no_reuse_initial,subgroup_order=r,P=point_key(P),Q=point_key(Q),frobenius_eigenvalue=lam,arms=[dict(name=arm['name'],basis=arm['basis'],offset=arm['offset'],exponents=arm['exponents'],column_count=len(states[arm['name']]['columns']),rank=len(states[arm['name']]['rank'].rows),seconds=states[arm['name']]['seconds'],extraction=states[arm['name']]['extraction'],statuses={s:sum(row['status']==s for row in states[arm['name']]['runs']) for s in sorted(set(row['status'] for row in states[arm['name']]['runs']))},runs=states[arm['name']]['runs']) for arm in arms])
    order_rng=random.Random(args.target_seed*100+19)
    start_index=0 if args.no_reuse_initial else 12
    for index in range(start_index,len(targets)):
        active=[arm for arm in arms if states[arm['name']]['extraction'] is None]
        if not active:break
        order_rng.shuffle(active);a,b,R=targets[index]
        for arm in active:
            state=states[arm['name']];start=time.perf_counter();result=wdsat_adapter.solve(E,arm['bases'],arm['offsets'],R,r,args.solver,config['watchdog_seconds'],args.output/'raw'/arm['name']/str(index))
            reachable=bool(state['coverage_mask'] & (1 << ((a+b*secret)%r)))
            if result['status']=='unsat' and reachable:raise RuntimeError('WDSat UNSAT contradicts reachability')
            if result['status']=='verified_relation':
                if not reachable:raise RuntimeError('Relation outside reachability')
                coefficients=rewritten_coefficients(E,result['witness'],arm['exponents'],lam,r,state['columns']);lhs=E(0)
                for key,pos in state['columns'].items():lhs+=coefficients[pos]*E(F.from_integer(key[0]),F.from_integer(key[1]))
                if lhs!=R:raise RuntimeError('Rewritten relation failed')
                d=state['rank'].add(coefficients,a,b);result['coefficients']=coefficients
                if d is not None:
                    if d!=secret or d*P!=Q:raise RuntimeError('Extracted log failed')
                    state['extraction']=dict(log=d,target_index=index,relations=sum(row['status']=='verified_relation' for row in state['runs'])+1,seconds=state['seconds']+time.perf_counter()-start);result['verified_extracted_log']=d
            result.update(target_index=index,a=a,b=b,reused=False,reachable=reachable,rank_after=len(state['rank'].rows),measured_seconds=time.perf_counter()-start)
            state['runs'].append(result);state['seconds']+=result['measured_seconds'];(args.output/'progress.json').write_text(json.dumps(report(),indent=2)+'\n')
            print(json.dumps(dict(arm=arm['name'],target=index,status=result['status'],rank=result['rank_after'],extracted=state['extraction'] is not None,seconds=result['measured_seconds'])),flush=True)
    final=report();(args.output/'result.json').write_text(json.dumps(final,indent=2)+'\n');(args.output/'manifest.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),adapter_sha256=hashlib.sha256(Path(wdsat_adapter.__file__).read_bytes()).hexdigest(),config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest()),indent=2)+'\n')


if __name__=='__main__':main()
