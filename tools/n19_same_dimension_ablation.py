"""Three-stream n19 same-dimension factor-base extraction ablation."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import subprocess
import time

from sage.all import GF, EllipticCurve
from sage.env import SAGE_VERSION
from factor_base_sat_rank import RelationRank, columns, point_key
from factor_base_frobenius_sat import rewritten_coefficients, koblitz_frobenius_eigenvalue
from factor_base_yield_screen import bit_span, subgroup_generator
import wdsat_adapter


def point_reachability(E,basis,offset,exponents,r):
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
    sums={E(0)}
    for points in sets:sums={S+P for S in sums for P in points}
    return sums-{E(0)}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,required=True);parser.add_argument('--solver',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.solver=args.solver.resolve();config_path=args.source/'n19-same-dimension-ablation.json';config=json.loads(config_path.read_text());frontier=json.loads((args.source/'dimension-frontier-n19-20260909-a'/'result.json').read_text())
    args.output.mkdir(parents=True,exist_ok=False)
    for path in [Path(__file__),Path(wdsat_adapter.__file__),args.solver,config_path]:(args.output/path.name).write_bytes(path.read_bytes())
    F=GF(2**19,name='a');E=EllipticCurve(F,[1,1,0,0,1]);r=int(frontier['subgroup_order']);P=subgroup_generator(E,r);lam=koblitz_frobenius_eigenvalue(E,P,r)
    arms=[]
    validation_start=time.perf_counter()
    for spec in config['arms']:
        candidate=next(row for row in frontier['candidates'] if row['name']==spec['candidate']);base=[F.from_integer(v) for v in candidate['basis']];offset=F.from_integer(candidate['offset']);exponents=spec['schedule']
        arm=dict(name=spec['name'],basis=candidate['basis'],offset=candidate['offset'],exponents=exponents,
                 bases=[[v**(2**e) for v in base] for e in exponents],offsets=[offset**(2**e) for e in exponents],
                 cofactor_multiplier=1,selection_seconds=candidate['selection_seconds'] or 0.0)
        arm['columns']=columns(E,r,arm);arm['reachable']=point_reachability(E,arm['basis'],arm['offset'],exponents,r);arms.append(arm)
    validation_setup_seconds=time.perf_counter()-validation_start
    streams=[]
    for stream_spec in config['streams']:
        seed=stream_spec['target_seed'];secret=stream_spec['secret'];Q=secret*P;rng=random.Random(seed+19);targets=[]
        while len(targets)<config['maximum_targets']:
            a,b=rng.randrange(r),rng.randrange(r);R=a*P+b*Q
            if R!=E(0):targets.append((a,b,R))
        states={arm['name']:dict(rank=RelationRank(r),runs=[],solver_seconds=0.0,pipeline_seconds=0.0,extraction=None) for arm in arms};order_rng=random.Random(seed*100+19)
        def stream_report():
            return dict(target_seed=seed,secret=secret,Q=point_key(Q),arms=[dict(name=arm['name'],basis=arm['basis'],offset=arm['offset'],exponents=arm['exponents'],column_count=len(arm['columns']),reachable_targets=len(arm['reachable']),selection_seconds=arm['selection_seconds'],rank=len(states[arm['name']]['rank'].rows),solver_seconds=states[arm['name']]['solver_seconds'],pipeline_seconds=states[arm['name']]['pipeline_seconds'],extraction=states[arm['name']]['extraction'],statuses={s:sum(row['status']==s for row in states[arm['name']]['runs']) for s in sorted(set(row['status'] for row in states[arm['name']]['runs']))},runs=states[arm['name']]['runs']) for arm in arms])
        for index,(a,b,R) in enumerate(targets):
            active=[arm for arm in arms if states[arm['name']]['extraction'] is None]
            if not active:break
            order_rng.shuffle(active)
            for arm in active:
                state=states[arm['name']];start=time.perf_counter();result=wdsat_adapter.solve(E,arm['bases'],arm['offsets'],R,r,args.solver,config['watchdog_seconds'],args.output/'raw'/str(seed)/arm['name']/str(index));reachable=R in arm['reachable']
                if result['status']=='unsat' and reachable:raise RuntimeError('WDSat UNSAT contradicts point-addition reachability')
                if result['status']=='verified_relation':
                    if not reachable:raise RuntimeError('Relation outside point-addition reachability')
                    coefficients=rewritten_coefficients(E,result['witness'],arm['exponents'],lam,r,arm['columns']);lhs=E(0)
                    for key,pos in arm['columns'].items():lhs+=coefficients[pos]*E(F.from_integer(key[0]),F.from_integer(key[1]))
                    if lhs!=R:raise RuntimeError('Rewritten relation failed')
                    d=state['rank'].add(coefficients,a,b);result['coefficients']=coefficients
                    if d is not None:
                        if d!=secret or d*P!=Q:raise RuntimeError('Extracted log failed')
                        state['extraction']=dict(log=d,target_index=index,relations=sum(row['status']=='verified_relation' for row in state['runs'])+1);result['verified_extracted_log']=d
                elapsed=time.perf_counter()-start;result.update(target_index=index,a=a,b=b,reachable=reachable,rank_after=len(state['rank'].rows),pipeline_seconds=elapsed)
                state['runs'].append(result);state['solver_seconds']+=result['solver_seconds'];state['pipeline_seconds']+=elapsed
                if state['extraction'] is not None:state['extraction'].update(solver_seconds=state['solver_seconds'],pipeline_seconds=state['pipeline_seconds'])
                (args.output/'progress.json').write_text(json.dumps(dict(streams=streams+[stream_report()]),indent=2)+'\n')
                print(json.dumps(dict(seed=seed,arm=arm['name'],target=index,status=result['status'],rank=result['rank_after'],extracted=state['extraction'] is not None,solver_seconds=result['solver_seconds'])),flush=True)
        streams.append(stream_report())
    environment=dict(platform=platform.platform(),machine=platform.machine(),cpu_count=os.cpu_count(),python=platform.python_version(),sage=SAGE_VERSION,
                     compiler=subprocess.run(['cc','--version'],capture_output=True,text=True,check=True).stdout.strip(),
                     solver_sha256=hashlib.sha256(args.solver.read_bytes()).hexdigest(),source_commit=subprocess.run(['git','rev-parse','HEAD'],capture_output=True,text=True,check=True).stdout.strip())
    result=dict(scope='n19 same-dimension factor-base implementation ablation',subgroup_order=r,P=point_key(P),frobenius_eigenvalue=lam,
                eigenvalue_method='roots of T^2-tT+2, selected by public point equality',validation_setup_seconds=validation_setup_seconds,
                environment=environment,streams=streams)
    (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n');(args.output/'manifest.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),adapter_sha256=hashlib.sha256(Path(wdsat_adapter.__file__).read_bytes()).hexdigest(),config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest()),indent=2)+'\n')


if __name__=='__main__':main()
