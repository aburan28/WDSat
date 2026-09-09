"""Exact n=17 dimension/basis/Frobenius-schedule Pareto screen."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

from sage.all import GF, EllipticCurve
import binary_factor_base as lab
from factor_base_frobenius_schedule import score_sets, group_coverage
from factor_base_sat_rank import point_key
from factor_base_yield_screen import bit_rank, bit_span, subgroup_generator


def equation_shape(F,E,basis,offset,schedule,P,r):
    base=[F.from_integer(v) for v in basis];off=F.from_integer(offset)
    bases=[[v**(2**e) for v in base] for e in schedule];offsets=[off**(2**e) for e in schedule]
    start=time.perf_counter();equations,names=lab.descend_multi(F,bases,P[0],E.a6(),offsets,E.a2().trace())
    nvars,clauses=lab.cnf(equations,len(names));factor_clauses,valid,invalid=lab.exact_factor_clauses_multi(E,bases,r,offsets)
    nonlinear={mon for eq in equations for mon in eq if mon.bit_count()>1}
    return dict(boolean_variables=len(names),anf_equations=len(equations),anf_terms=sum(map(len,equations)),
                nonlinear_monomials=len(nonlinear),cnf_variables=nvars,cnf_clauses=len(clauses)+len(factor_clauses),
                exact_factor_clauses=len(factor_clauses),valid_assignments=[len(row) for row in valid],
                construction_wall_seconds=time.perf_counter()-start)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--n',type=int,default=17);parser.add_argument('--dimensions',help='comma-separated override')
    args=parser.parse_args();config_path=args.source/'dimension-frontier-amendment.json';config=json.loads(config_path.read_text())
    args.output.mkdir(parents=True,exist_ok=False)
    for path in [Path(__file__),Path(lab.__file__),config_path]:(args.output/path.name).write_bytes(path.read_bytes())
    n=args.n;dimensions=[int(x) for x in args.dimensions.split(',')] if args.dimensions else config['dimensions']
    F=GF(2**n,name='a');E=EllipticCurve(F,[1,1,0,0,1]);r=int(max(p for p,_ in E.cardinality().factor()));P=subgroup_generator(E,r)
    if E.cardinality()!=2*r:raise ValueError('Frontier screen currently requires cofactor two')
    logs_by_x={};log_by_point={};Q=E(0)
    for j in range(r):
        log_by_point[point_key(Q)]=j
        if j and j<=(r-1)//2:logs_by_x[int(Q[0].to_integer())]=j
        Q+=P
    lam=log_by_point[point_key(E(P[0]**2,P[1]**2))]
    directions={i:int((F.gen()**i+F((F.gen()**i).trace())).to_integer()) for i in range(1,n)}
    candidates=[]
    for k in dimensions:
        polynomial=[1 << i for i in range(k)]
        trace=[directions[i] for i in range(1,k+1)]
        candidates.extend([dict(name=f'polynomial-k{k}',family='polynomial',dimension=k,basis=polynomial,offset=0,selection_seconds=0.0),
                           dict(name=f'trace-k{k}',family='trace',dimension=k,basis=trace,offset=1,selection_seconds=0.0)])
        start=time.perf_counter();best=None;count=0
        for exponents in itertools.combinations(range(1,n),k):
            basis=[directions[i] for i in exponents];valid=sum(x in logs_by_x for x in bit_span(basis,1));count+=1
            item=(valid,-sum(exponents),tuple(-i for i in exponents),basis,exponents)
            if best is None or item[:3]>best[:3]:best=item
        candidates.append(dict(name=f'sparse-k{k}',family='sparse',dimension=k,basis=best[3],offset=1,
            selection_seconds=time.perf_counter()-start,selection_candidates=count,selection_exponents=list(best[4])))
    if n==17 and 5 in dimensions:
        candidates.append(dict(name='all-lift-flat-k5',family='exploratory-flat',dimension=5,
            basis=[70724,4642,8208,34008,2066],offset=48581,selection_seconds=None))
    rows=[]
    for candidate in candidates:
        if bit_rank(candidate['basis'])!=candidate['dimension']:raise RuntimeError('Dependent factor basis')
        reps=[logs_by_x[x] for x in bit_span(candidate['basis'],candidate['offset']) if x in logs_by_x]
        scalars=[v for s in reps for v in (s,r-s)];schedule_scores=[]
        for e in range(n):
            for f in range(e,n):
                schedule=[0,e,f];sets=[[(pow(lam,j,r)*s)%r for s in scalars] for j in schedule]
                mask,coverage=score_sets(sets,r);schedule_scores.append(dict(schedule=schedule,coverage=coverage))
        best_schedule=min(schedule_scores,key=lambda row:(-row['coverage'],row['schedule']))
        controls=[]
        for name,schedule in [('same',[0,0,0]),('consecutive',[0,1,2]),('best',best_schedule['schedule'])]:
            sets=[[(pow(lam,j,r)*s)%r for s in scalars] for j in schedule]
            _,coverage=score_sets(sets,r);fraction=coverage/(r-1)
            controls.append(dict(name=name,schedule=schedule,coverage=coverage,coverage_fraction=fraction,
                                 expected_systems_for_column_count=len(reps)/fraction if fraction else None,
                                 encoding=equation_shape(F,E,candidate['basis'],candidate['offset'],schedule,P,r)))
        rows.append(candidate|dict(actual_factor_points=len(scalars),sign_folded_columns=len(reps),
                                   schedule_scores=schedule_scores,controls=controls))
        print(json.dumps(dict(name=candidate['name'],columns=len(reps),
            controls=[{k:c[k] for k in ['name','schedule','coverage','expected_systems_for_column_count']} for c in controls])),flush=True)
    points=[]
    for row in rows:
        control=next(c for c in row['controls'] if c['name']=='best')
        point=dict(row=row,control=control)
        points.append(point)
    frontier=[]
    for point in points:
        metrics=(point['control']['expected_systems_for_column_count'],point['control']['encoding']['anf_terms'],
                 point['control']['encoding']['cnf_clauses'])
        dominated=False
        for other in points:
            if other is point:continue
            om=(other['control']['expected_systems_for_column_count'],other['control']['encoding']['anf_terms'],
                other['control']['encoding']['cnf_clauses'])
            if all(a<=b for a,b in zip(om,metrics)) and any(a<b for a,b in zip(om,metrics)):
                dominated=True;break
        if not dominated:frontier.append(point['row']['name'])
    # Independent curve-addition check for every reported frontier point.
    for name in frontier:
        row=next(row for row in rows if row['name']==name);control=next(c for c in row['controls'] if c['name']=='best')
        base_points=set()
        for x in bit_span(row['basis'],row['offset']):
            try:lifts=E.lift_x(F.from_integer(x),all=True)
            except ValueError:lifts=[]
            base_points.update(Q for Q in lifts if r*Q==E(0))
        point_sets=[]
        for e in control['schedule']:
            point_sets.append({E(Q[0]**(2**e),Q[1]**(2**e)) for Q in base_points})
        control['independent_group_coverage']=group_coverage(E,point_sets)
        if control['independent_group_coverage']!=control['coverage']:raise RuntimeError('Frontier group coverage mismatch')
    result=dict(scope=f'n={n} exact Pareto screen; not a breakthrough claim',n=n,field_modulus=str(F.modulus()),
                subgroup_order=r,frobenius_eigenvalue=lam,frontier=frontier,candidates=rows)
    (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    (args.output/'manifest.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest()),indent=2)+'\n')
    print(json.dumps(dict(frontier=frontier)))


if __name__=='__main__':main()
