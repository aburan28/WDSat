"""Exact Frobenius-image schedule screen with a group-law-only verification."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import random
import time

from sage.all import GF, EllipticCurve
from factor_base_yield_screen import bit_span, subgroup_generator
from factor_base_sat_rank import point_key


def add_mask(mask, values, r):
    result=0
    for value in values:
        result |= (mask << value) | (mask >> (r-value))
    return result & ((1 << r)-1)


def score_sets(sets, r):
    mask=1
    for values in sets:
        mask=add_mask(mask,values,r)
    return mask, (mask & ~1).bit_count()


def point_factor_base(E, basis, offset, r):
    F=E.base_field(); points=set()
    for x in bit_span(basis,offset):
        try: lifts=E.lift_x(F.from_integer(x),all=True)
        except ValueError: lifts=[]
        points.update(P for P in lifts if r*P==E(0))
    return points


def group_coverage(E, point_sets):
    sums={E(0)}
    for points in point_sets:
        sums={S+P for S in sums for P in points}
    return len(sums-{E(0)})


def screen(n, root, config, output):
    sparse=json.loads((root/'sparse-20260909-a'/f'n{n}-result.json').read_text())
    F=GF(2**n,name='a'); E=EllipticCurve(F,[1,1,0,0,1]); r=int(sparse['subgroup_order'])
    P=subgroup_generator(E,r); base=point_factor_base(E,sparse['basis'],sparse['offset'],r)
    log_by_point={}; Q=P
    oracle_start=time.perf_counter()
    for scalar in range(r):
        log_by_point[point_key(Q)]=scalar; Q+=P
    if len(log_by_point)!=r: raise RuntimeError('Subgroup log oracle incomplete')
    frob=E(P[0]**2,P[1]**2); lam=log_by_point[point_key(frob)]
    if frob!=lam*P: raise RuntimeError('Frobenius eigenvalue verification failed')
    scalars=sorted(log_by_point[point_key(Q)] for Q in base)
    scores=[]
    for e in range(n):
        for f in range(e,n):
            schedule=[0,e,f]
            sets=[[(pow(lam,j,r)*s)%r for s in scalars] for j in schedule]
            mask,coverage=score_sets(sets,r)
            scores.append(dict(schedule=schedule,coverage=coverage,coverage_mask_hex=hex(mask)))
    selected=min(scores,key=lambda row:(-row['coverage'],row['schedule']))
    controls=[]
    for name,schedule in [('same',[0,0,0]),('consecutive',[0,1,2]),('selected',selected['schedule'])]:
        sets=[[(pow(lam,j,r)*s)%r for s in scalars] for j in schedule]
        mask,coverage=score_sets(sets,r)
        point_sets=[]
        for j in schedule:
            current=set(base)
            for _ in range(j): current={E(Q[0]**2,Q[1]**2) for Q in current}
            point_sets.append(current)
        independent=group_coverage(E,point_sets)
        if independent!=coverage: raise RuntimeError('Group-law coverage disagrees with cyclic oracle')
        controls.append(dict(name=name,schedule=schedule,coverage=coverage,
                             independent_group_coverage=independent,coverage_mask_hex=hex(mask)))
    rng=random.Random(202609097+n); nulls=[]
    for _ in range(32):
        sets=[]
        for _ in range(3):
            half=rng.sample(range(1,(r+1)//2),len(base)//2)
            sets.append([v for s in half for v in (s,r-s)])
        nulls.append(score_sets(sets,r)[1])
    result=dict(scope='prior-art optimization screen; not a novelty or breakthrough claim',n=n,
                field_modulus=str(F.modulus()),subgroup_order=r,basis=sparse['basis'],offset=sparse['offset'],
                factor_points=len(base),frobenius_eigenvalue=lam,scores=scores,selected=selected,
                controls=controls,matched_random_coverages=nulls,
                oracle_and_screen_wall_seconds=time.perf_counter()-oracle_start)
    (output/f'n{n}-result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(n=n,factor_points=len(base),selected=selected,controls=controls,
                         null_mean=sum(nulls)/len(nulls))),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();config_path=args.source/'frobenius-schedule-amendment.json'
    config=json.loads(config_path.read_text());args.output.mkdir(parents=True,exist_ok=False)
    for path in [Path(__file__),config_path]: (args.output/path.name).write_bytes(path.read_bytes())
    (args.output/'manifest.json').write_text(json.dumps(dict(
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest()),indent=2)+'\n')
    for n in config['profiles']: screen(n,args.source,config,args.output)


if __name__=='__main__':main()
