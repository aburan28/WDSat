"""Paired SAT/rank comparison for same versus Frobenius-image factor spaces."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import time

from sage.all import GF, EllipticCurve
import binary_factor_base as lab
from factor_base_sat_rank import RelationRank, columns, point_key
from factor_base_yield_screen import subgroup_generator


def frobenius_power(E, point, exponent):
    for _ in range(exponent):
        point=E(point[0]**2,point[1]**2)
    return point


def koblitz_frobenius_eigenvalue(E, point, subgroup_order):
    """Compute pi's subgroup eigenvalue from T^2-tT+2, then identify its root."""
    if E.a2()**2 != E.a2() or E.a6()**2 != E.a6():
        raise ValueError('Curve coefficients must lie in GF(2)')
    base_curve=EllipticCurve(GF(2),[1,int(E.a2()),0,0,int(E.a6())])
    trace=3-int(base_curve.cardinality())
    field=GF(subgroup_order);root=field(trace*trace-8).sqrt();inv2=field(2)**-1
    roots=[int((field(trace)+root)*inv2),int((field(trace)-root)*inv2)]
    image=frobenius_power(E,point,1)
    matches=[value for value in roots if value*point==image]
    if len(matches)!=1:raise RuntimeError('Could not identify a unique Frobenius eigenvalue')
    value=matches[0]
    if (value*value-trace*value+2)%subgroup_order:raise RuntimeError('Invalid Frobenius characteristic root')
    return value


def rewritten_coefficients(E,witness,exponents,lam,r,column_map):
    coefficients=[0]*len(column_map); n=E.base_field().degree()
    for pair,e in zip(witness,exponents):
        point=E(E.base_field().from_integer(pair[0]),E.base_field().from_integer(pair[1]))
        base_point=frobenius_power(E,point,(-e)%n)
        if frobenius_power(E,base_point,e)!=point:
            raise RuntimeError('Inverse Frobenius rewrite failed')
        pos,neg=point_key(base_point),point_key(-base_point);key=min(pos,neg)
        coefficients[column_map[key]] += (1 if pos==key else -1)*pow(lam,e,r)
    return coefficients


def run_profile(n,source,output,config,base_config):
    sparse=json.loads((source/'sparse-20260909-a'/f'n{n}-result.json').read_text())
    profile=sparse['profile'];F=GF(2**n,name='a');E=EllipticCurve(F,[1,1,0,0,1]);r=int(sparse['subgroup_order'])
    P=subgroup_generator(E,r);Q=(base_config['simulation_secret']%r)*P
    # Compute and verify lambda without using it in SAT.
    piP=frobenius_power(E,P,1);walk=E(0);lam=None
    for j in range(r):
        if walk==piP:lam=j;break
        walk+=P
    if lam is None:raise RuntimeError('Frobenius eigenvalue not found')
    base=[F.from_integer(v) for v in sparse['basis']];offset=F.from_integer(sparse['offset'])
    setting=next(row for row in config['profiles'] if row['n']==n)
    arms=[]
    for name,exponents in [('same',[0,0,0]),('consecutive',[0,1,2])]:
        bases=[];offsets=[]
        for e in exponents:
            bases.append([v**(2**e) for v in base]);offsets.append(offset**(2**e))
        arms.append(dict(name=name,exponents=exponents,bases=bases,offsets=offsets,
                         basis=[int(v.to_integer()) for v in base],offset=int(offset.to_integer()),
                         cofactor_multiplier=1))
    column_map=columns(E,r,arms[0])
    states={arm['name']:dict(rank=RelationRank(r),runs=[],seconds=0.0,first_extraction=None) for arm in arms}
    rng=random.Random(base_config['target_seed']+n);targets=[]
    while len(targets)<setting['targets']:
        a,b=rng.randrange(r),rng.randrange(r);R=a*P+b*Q
        if R!=E(0):targets.append((a,b,R))
    root=output/f'n{n}';root.mkdir()
    (root/'targets.json').write_text(json.dumps([dict(a=a,b=b,target=point_key(R)) for a,b,R in targets],indent=2)+'\n')
    def report():
        rows=[]
        for arm in arms:
            state=states[arm['name']]
            rows.append(dict(name=arm['name'],exponents=arm['exponents'],rank=len(state['rank'].rows),
                column_count=len(column_map),measured_seconds=state['seconds'],first_extraction=state['first_extraction'],
                statuses={s:sum(run['status']==s for run in state['runs'])
                          for s in sorted(set(run['status'] for run in state['runs']))},runs=state['runs']))
        return dict(scope='paired prior-art Frobenius SAT diagnostic; not a breakthrough claim',profile=profile,
                    subgroup_order=r,frobenius_eigenvalue=lam,P=point_key(P),Q=point_key(Q),arms=rows)
    for index,(a,b,R) in enumerate(targets):
        order=list(arms);rng.shuffle(order)
        for arm in order:
            state=states[arm['name']];start=time.perf_counter()
            trial=lab.solve_verified(E,base,3,R,r,'cadical',setting['watchdog_seconds'],offset,
                trace_filter=True,exact_factor_filter=True,factor_bases=arm['bases'],factor_offsets=arm['offsets'],
                artifact_dir=root/arm['name']/str(index))
            trial.update(target_index=index,a=a,b=b,rank_before=len(state['rank'].rows))
            if trial['rejected_projections']:raise RuntimeError('Exact filter admitted invalid projection')
            if trial['status']=='verified_relation':
                coefficients=rewritten_coefficients(E,trial['witness'],arm['exponents'],lam,r,column_map)
                lhs=E(0)
                for key,pos in column_map.items():lhs += coefficients[pos]*E(F.from_integer(key[0]),F.from_integer(key[1]))
                if lhs!=R:raise RuntimeError('Rewritten point relation failed')
                d=state['rank'].add(coefficients,a,b);trial['coefficients']=coefficients
                if d is not None:
                    if d*P!=Q:raise RuntimeError('Extracted log failed')
                    trial['verified_extracted_log']=d
                    if state['first_extraction'] is None:state['first_extraction']=dict(target_index=index,log=d)
            trial['rank_after']=len(state['rank'].rows);trial['measured_wall_seconds']=time.perf_counter()-start
            state['seconds']+=trial['measured_wall_seconds'];state['runs'].append(trial)
            (root/'progress.json').write_text(json.dumps(report(),indent=2)+'\n')
            print(json.dumps(dict(n=n,arm=arm['name'],target=index,status=trial['status'],
                                  rank=trial['rank_after'],seconds=trial['measured_wall_seconds'])),flush=True)
    (output/f'n{n}-result.json').write_text(json.dumps(report(),indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();config_path=args.source/'frobenius-sat-amendment.json';config=json.loads(config_path.read_text())
    base_config=json.loads((args.source/'sat-rank-amendment.json').read_text());args.output.mkdir(parents=True,exist_ok=False)
    for path in [Path(__file__),Path(lab.__file__),config_path]:(args.output/path.name).write_bytes(path.read_bytes())
    (args.output/'manifest.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest()),indent=2)+'\n')
    for setting in config['profiles']:run_profile(setting['n'],args.source,args.output,config,base_config)


if __name__=='__main__':main()
