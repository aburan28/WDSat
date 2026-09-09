"""End-to-end extraction for the n=17 dimension/schedule frontier."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import time

from sage.all import GF, EllipticCurve
import binary_factor_base as lab
from factor_base_sat_rank import RelationRank, columns, point_key
from factor_base_frobenius_sat import rewritten_coefficients, frobenius_power
from factor_base_yield_screen import subgroup_generator


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();config_path=args.source/'frontier-extraction-amendment.json';config=json.loads(config_path.read_text())
    extraction=json.loads((args.source/'extraction-replication.json').read_text());frontier=json.loads((args.source/'dimension-frontier-20260909-a'/'result.json').read_text())
    args.output.mkdir(parents=True,exist_ok=False)
    for path in [Path(__file__),Path(lab.__file__),config_path]:(args.output/path.name).write_bytes(path.read_bytes())
    (args.output/'manifest.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest()),indent=2)+'\n')
    F=GF(2**17,name='a');E=EllipticCurve(F,[1,1,0,0,1]);r=int(frontier['subgroup_order']);P=subgroup_generator(E,r);secret=extraction['simulation_secret']%r;Q=secret*P
    piP=frobenius_power(E,P,1);walk=E(0);lam=None
    for j in range(r):
        if walk==piP:lam=j;break
        walk+=P
    if lam is None:raise RuntimeError('Frobenius eigenvalue not found')
    requested=[('sparse-k5','best'),('sparse-k5','consecutive'),('trace-k6','best'),('polynomial-k6','best'),('sparse-k4','best')]
    arms=[];states={}
    for candidate_name,control_name in requested:
        candidate=next(row for row in frontier['candidates'] if row['name']==candidate_name)
        control=next(row for row in candidate['controls'] if row['name']==control_name)
        exponents=control['schedule'];base=[F.from_integer(v) for v in candidate['basis']];offset=F.from_integer(candidate['offset'])
        arm=dict(name=candidate_name+'-'+control_name,basis=candidate['basis'],offset=candidate['offset'],
                 exponents=exponents,bases=[[v**(2**e) for v in base] for e in exponents],
                 offsets=[offset**(2**e) for e in exponents],cofactor_multiplier=1,
                 selection_seconds=candidate['selection_seconds'] or 0.0,
                 predicted_coverage=control['coverage'],predicted_systems=control['expected_systems_for_column_count'])
        arms.append(arm);states[arm['name']]=dict(columns=columns(E,r,arm),rank=RelationRank(r),runs=[],seconds=0.0,extraction=None)
    seed=extraction['target_seeds'][0];rng=random.Random(seed+17);targets=[]
    while len(targets)<96:
        a,b=rng.randrange(r),rng.randrange(r);R=a*P+b*Q
        if R!=E(0):targets.append((a,b,R))
    (args.output/'targets.json').write_text(json.dumps([dict(index=i,a=a,b=b,target=point_key(R)) for i,(a,b,R) in enumerate(targets)],indent=2)+'\n')
    order_rng=random.Random(seed+1700)
    def report():
        rows=[]
        for arm in arms:
            state=states[arm['name']]
            rows.append(dict(name=arm['name'],basis=arm['basis'],offset=arm['offset'],exponents=arm['exponents'],
                column_count=len(state['columns']),predicted_coverage=arm['predicted_coverage'],predicted_systems=arm['predicted_systems'],
                selection_seconds=arm['selection_seconds'],rank=len(state['rank'].rows),measured_seconds=state['seconds'],
                charged_seconds=state['seconds']+arm['selection_seconds'],extraction=state['extraction'],
                statuses={s:sum(run['status']==s for run in state['runs']) for s in sorted(set(run['status'] for run in state['runs']))},runs=state['runs']))
        return dict(scope='n=17 frontier extraction; no broad breakthrough claim',subgroup_order=r,P=point_key(P),Q=point_key(Q),frobenius_eigenvalue=lam,arms=rows)
    for index,(a,b,R) in enumerate(targets):
        active=[arm for arm in arms if states[arm['name']]['extraction'] is None]
        if not active:break
        order_rng.shuffle(active)
        for arm in active:
            state=states[arm['name']];start=time.perf_counter();base=[F.from_integer(v) for v in arm['basis']]
            trial=lab.solve_verified(E,base,3,R,r,'cadical',10,F.from_integer(arm['offset']),trace_filter=True,
                exact_factor_filter=True,factor_bases=arm['bases'],factor_offsets=arm['offsets'],artifact_dir=args.output/'raw'/arm['name']/str(index))
            trial.update(target_index=index,a=a,b=b,rank_before=len(state['rank'].rows))
            if trial['rejected_projections']:raise RuntimeError('Invalid factor projection')
            if trial['status']=='verified_relation':
                coefficients=rewritten_coefficients(E,trial['witness'],arm['exponents'],lam,r,state['columns']);lhs=E(0)
                for key,pos in state['columns'].items():lhs+=coefficients[pos]*E(F.from_integer(key[0]),F.from_integer(key[1]))
                if lhs!=R:raise RuntimeError('Rewritten relation failed')
                d=state['rank'].add(coefficients,a,b);trial['coefficients']=coefficients
                if d is not None:
                    if d!=secret or d*P!=Q:raise RuntimeError('Extracted log failed')
                    state['extraction']=dict(log=d,target_index=index,relations=sum(run['status']=='verified_relation' for run in state['runs'])+1,
                                             measured_seconds=state['seconds']+time.perf_counter()-start);trial['verified_extracted_log']=d
            trial['rank_after']=len(state['rank'].rows);trial['measured_wall_seconds']=time.perf_counter()-start
            state['seconds']+=trial['measured_wall_seconds'];state['runs'].append(trial);(args.output/'progress.json').write_text(json.dumps(report(),indent=2)+'\n')
            print(json.dumps(dict(arm=arm['name'],target=index,status=trial['status'],rank=trial['rank_after'],extracted=state['extraction'] is not None,seconds=trial['measured_wall_seconds'])),flush=True)
    final=report();(args.output/'result.json').write_text(json.dumps(final,indent=2)+'\n')
    if not all(row['extraction'] for row in final['arms']):print('Some frontier arms did not extract within 96 targets')


if __name__=='__main__':main()
