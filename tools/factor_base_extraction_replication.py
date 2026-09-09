"""End-to-end toy relation collection through verified log extraction."""
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
    args=parser.parse_args();config_path=args.source/'extraction-replication.json';config=json.loads(config_path.read_text())
    args.output.mkdir(parents=True,exist_ok=False)
    for path in [Path(__file__),Path(lab.__file__),config_path]: (args.output/path.name).write_bytes(path.read_bytes())
    (args.output/'manifest.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest()),indent=2)+'\n')
    profile=config['profile'];n=profile['n'];F=GF(2**n,name='a');E=EllipticCurve(F,[1,1,0,0,1])
    screened=json.loads((args.source/'run-20260909-a'/f'n{n}-results.json').read_text())
    sparse_result=json.loads((args.source/'sparse-20260909-a'/f'n{n}-result.json').read_text())
    sparse_selection=json.loads((args.source/'sparse-20260909-a'/f'n{n}-selection.json').read_text())
    r=int(screened['subgroup_order']);P=subgroup_generator(E,r);Q=(config['simulation_secret']%r)*P
    piP=frobenius_power(E,P,1);walk=E(0);lam=None
    for j in range(r):
        if walk==piP:lam=j;break
        walk+=P
    if lam is None or lam*P!=piP:raise RuntimeError('Frobenius eigenvalue failed')
    trace=next(row for row in screened['results'] if row['name']=='trace-zero-affine')
    specs=[('polynomial',[1 << i for i in range(profile['dimension'])],0,[0,0,0],0.0),
           ('trace',trace['basis'],1,[0,0,0],0.0),
           ('sparse-same',sparse_result['basis'],1,[0,0,0],sparse_selection['selection_wall_seconds']),
           ('sparse-frobenius',sparse_result['basis'],1,[0,1,2],sparse_selection['selection_wall_seconds'])]
    arms=[];states={}
    for name,basis_int,offset_int,exponents,selection_seconds in specs:
        base=[F.from_integer(v) for v in basis_int];offset=F.from_integer(offset_int)
        bases=[[v**(2**e) for v in base] for e in exponents];offsets=[offset**(2**e) for e in exponents]
        candidate=dict(name=name,basis=basis_int,offset=offset_int,cofactor_multiplier=1,
                       exponents=exponents,bases=bases,offsets=offsets,selection_seconds=selection_seconds)
        arms.append(candidate);col=columns(E,r,candidate)
        states[name]=dict(columns=col,rank=RelationRank(r),runs=[],measured_seconds=0.0,extraction=None)
    targets=[]
    for seed in config['target_seeds']:
        rng=random.Random(seed+n)
        count=0
        while count<config['targets_per_seed']:
            a,b=rng.randrange(r),rng.randrange(r);R=a*P+b*Q
            if R==E(0):continue
            targets.append((seed,count,a,b,R));count+=1
    (args.output/'targets.json').write_text(json.dumps([dict(seed=s,index=i,a=a,b=b,target=point_key(R)) for s,i,a,b,R in targets],indent=2)+'\n')
    order_rng=random.Random(sum(config['target_seeds'])+n)
    def report():
        rows=[]
        for arm in arms:
            state=states[arm['name']]
            rows.append(dict(name=arm['name'],basis=arm['basis'],offset=arm['offset'],exponents=arm['exponents'],
                column_count=len(state['columns']),rank=len(state['rank'].rows),selection_seconds=arm['selection_seconds'],
                measured_seconds=state['measured_seconds'],total_charged_seconds=state['measured_seconds']+arm['selection_seconds'],
                extraction=state['extraction'],statuses={s:sum(run['status']==s for run in state['runs'])
                    for s in sorted(set(run['status'] for run in state['runs']))},runs=state['runs']))
        return dict(scope='toy end-to-end extraction replication; no broad cryptanalytic claim',profile=profile,
                    subgroup_order=r,P=point_key(P),Q=point_key(Q),frobenius_eigenvalue=lam,arms=rows)
    for global_index,(seed,index,a,b,R) in enumerate(targets):
        active=[arm for arm in arms if states[arm['name']]['extraction'] is None]
        if not active:break
        order_rng.shuffle(active)
        for arm in active:
            state=states[arm['name']];start=time.perf_counter();base=[F.from_integer(v) for v in arm['basis']]
            trial=lab.solve_verified(E,base,3,R,r,'cadical',config['solver_watchdog_seconds_per_target'],F.from_integer(arm['offset']),
                trace_filter=True,exact_factor_filter=True,factor_bases=arm['bases'],factor_offsets=arm['offsets'],
                artifact_dir=args.output/'raw'/arm['name']/str(global_index))
            trial.update(seed=seed,seed_index=index,global_target_index=global_index,a=a,b=b,rank_before=len(state['rank'].rows))
            if trial['rejected_projections']:raise RuntimeError('Exact factor predicate admitted invalid projection')
            if trial['status']=='verified_relation':
                coefficients=rewritten_coefficients(E,trial['witness'],arm['exponents'],lam,r,state['columns'])
                lhs=E(0)
                for key,pos in state['columns'].items():lhs+=coefficients[pos]*E(F.from_integer(key[0]),F.from_integer(key[1]))
                if lhs!=R:raise RuntimeError('Rewritten relation failed')
                d=state['rank'].add(coefficients,a,b);trial['coefficients']=coefficients
                if d is not None:
                    if d*P!=Q or d!=config['simulation_secret']%r:raise RuntimeError('Extracted log verification failed')
                    state['extraction']=dict(log=d,seed=seed,seed_index=index,global_target_index=global_index,
                                             relations=sum(run['status']=='verified_relation' for run in state['runs'])+1,
                                             measured_seconds=state['measured_seconds']+time.perf_counter()-start)
                    trial['verified_extracted_log']=d
            trial['rank_after']=len(state['rank'].rows);trial['measured_wall_seconds']=time.perf_counter()-start
            state['measured_seconds']+=trial['measured_wall_seconds'];state['runs'].append(trial)
            (args.output/'progress.json').write_text(json.dumps(report(),indent=2)+'\n')
            print(json.dumps(dict(arm=arm['name'],target=global_index,status=trial['status'],rank=trial['rank_after'],
                                  extracted=state['extraction'] is not None,seconds=trial['measured_wall_seconds'])),flush=True)
    final=report();(args.output/'result.json').write_text(json.dumps(final,indent=2)+'\n')
    if not all(row['extraction'] for row in final['arms']):print('One or more arms did not extract within the frozen target maximum',flush=True)


if __name__=='__main__':main()
