"""Paired exact-lift CNF filter benchmark on frozen factor-base candidates."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import time

from sage.all import GF, EllipticCurve
import binary_factor_base as lab
from factor_base_sat_rank import RelationRank, columns, relation_coefficients, point_key
from factor_base_yield_screen import subgroup_generator


def run_profile(n, base, output, base_config, filter_config):
    screened = json.loads((base/'run-20260909-a'/f'n{n}-results.json').read_text())
    sparse = json.loads((base/'sparse-20260909-a'/f'n{n}-result.json').read_text())
    profile = screened['profile']
    F = GF(2**n, name='a')
    E = EllipticCurve(F, [1, 1, 0, 0, 1])
    r = int(screened['subgroup_order'])
    P = subgroup_generator(E, r)
    Q = (base_config['simulation_secret'] % r)*P
    trace = next(row for row in screened['results'] if row['name'] == 'trace-zero-affine')
    families = [('trace', trace['basis']), ('sparse', sparse['basis'])]
    arms = []
    for name, basis in families:
        for filtered in [False, True]:
            arms.append(dict(name=name+('-filtered' if filtered else '-unfiltered'),
                             family=name, basis=basis, offset=1, cofactor_multiplier=1,
                             exact_factor_filter=filtered))
    setting = next(row for row in filter_config['profiles'] if row['n'] == n)
    rng = random.Random(base_config['target_seed']+n)
    targets = []
    while len(targets) < setting['targets']:
        a, b = rng.randrange(r), rng.randrange(r)
        R = a*P+b*Q
        if R != E(0):
            targets.append((a, b, R))
    states = {}
    for arm in arms:
        states[arm['name']] = dict(columns=columns(E, r, arm), rank=RelationRank(r),
                                  runs=[], seconds=0.0, first_extraction=None)
    root = output/f'n{n}'
    root.mkdir()
    (root/'targets.json').write_text(json.dumps([dict(a=a,b=b,target=point_key(R)) for a,b,R in targets],indent=2)+'\n')

    def make_report():
        rows=[]
        for arm in arms:
            state=states[arm['name']]
            statuses={s:sum(run['status']==s for run in state['runs'])
                      for s in sorted(set(run['status'] for run in state['runs']))}
            rows.append(arm | dict(column_count=len(state['columns']), rank=len(state['rank'].rows),
                statuses=statuses, measured_seconds=state['seconds'], first_extraction=state['first_extraction'],
                runs=state['runs']))
        return dict(scope='paired exact-lift filter diagnostic; not a breakthrough claim', profile=profile,
                    subgroup_order=r, P=point_key(P), Q=point_key(Q), arms=rows)

    for index,(a,b,R) in enumerate(targets):
        order=list(arms); rng.shuffle(order)
        for arm in order:
            state=states[arm['name']]
            start=time.perf_counter()
            trial=lab.solve_verified(E,[F.from_integer(v) for v in arm['basis']],3,R,r,'cadical',
                setting['watchdog_seconds'],F(1),trace_filter=True,
                artifact_dir=root/arm['name']/str(index), exact_factor_filter=arm['exact_factor_filter'])
            trial.update(target_index=index,a=a,b=b,rank_before=len(state['rank'].rows))
            if arm['exact_factor_filter'] and trial['rejected_projections']:
                raise RuntimeError('Exact factor filter admitted an invalid coordinate projection')
            if trial['status']=='verified_relation':
                coefficients,killed=relation_coefficients(E,trial['witness'],1,state['columns'])
                lhs=E(0)
                for key,pos in state['columns'].items():
                    lhs += coefficients[pos]*E(F.from_integer(key[0]),F.from_integer(key[1]))
                if lhs != R:
                    raise RuntimeError('Certified relation does not sum to target')
                d=state['rank'].add(coefficients,a,b)
                trial.update(coefficients=coefficients,killed_torsion_points=killed)
                if d is not None:
                    if d*P != Q: raise RuntimeError('Extracted log failed point verification')
                    trial['verified_extracted_log']=d
                    if state['first_extraction'] is None:
                        state['first_extraction']=dict(target_index=index,log=d)
            trial['rank_after']=len(state['rank'].rows)
            trial['measured_wall_seconds']=time.perf_counter()-start
            state['seconds'] += trial['measured_wall_seconds']; state['runs'].append(trial)
            (root/'progress.json').write_text(json.dumps(make_report(),indent=2)+'\n')
            print(json.dumps(dict(n=n,arm=arm['name'],target=index,status=trial['status'],
                                  rejects=trial['rejected_projections'],rank=trial['rank_after'],
                                  seconds=trial['measured_wall_seconds'])),flush=True)
    (output/f'n{n}-result.json').write_text(json.dumps(make_report(),indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    base_config=json.loads((args.source/'sat-rank-amendment.json').read_text())
    filter_path=args.source/'exact-lift-amendment.json'
    filter_config=json.loads(filter_path.read_text())
    args.output.mkdir(parents=True,exist_ok=False)
    for path in [Path(__file__),Path(lab.__file__),filter_path]:
        (args.output/path.name).write_bytes(path.read_bytes())
    (args.output/'manifest.json').write_text(json.dumps(dict(
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        filter_config_sha256=hashlib.sha256(filter_path.read_bytes()).hexdigest()),indent=2)+'\n')
    for setting in filter_config['profiles']:
        run_profile(setting['n'],args.source,args.output,base_config,filter_config)


if __name__=='__main__': main()
