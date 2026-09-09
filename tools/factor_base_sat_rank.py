"""Paired SAT relation collection and rank validation on frozen selected bases."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import time

from sage.all import GF, EllipticCurve
import binary_factor_base as lab
from factor_base_yield_screen import subgroup_generator, bit_span


class RelationRank:
    """Eliminate factor-base columns, retaining a,b in sum(c_i log F_i)=a+b*d."""
    def __init__(self, r):
        self.r = r
        self.rows = {}

    def add(self, coefficients, a, b):
        r = self.r
        row = [c % r for c in coefficients]
        a, b = a % r, b % r
        for pivot in sorted(self.rows):
            if row[pivot]:
                factor = row[pivot]
                old, aa, bb = self.rows[pivot]
                row = [(x-factor*y) % r for x, y in zip(row, old)]
                a, b = (a-factor*aa) % r, (b-factor*bb) % r
        pivot = next((j for j, value in enumerate(row) if value), None)
        if pivot is not None:
            inv = pow(row[pivot], -1, r)
            self.rows[pivot] = ([x*inv % r for x in row], a*inv % r, b*inv % r)
            return None
        if b:
            return -a*pow(b, -1, r) % r
        if a:
            raise RuntimeError('Inconsistent relation dependency')
        return None


def point_key(P):
    return int(P[0].to_integer()), int(P[1].to_integer())


def columns(E, r, candidate):
    by_key = {}
    h = candidate['cofactor_multiplier']
    for x in bit_span(candidate['basis'], candidate['offset']):
        try:
            points = E.lift_x(E.base_field().from_integer(x), all=True)
        except ValueError:
            continue
        for P in points:
            if h == 1 and r*P != E(0):
                continue
            image = h*P
            if image == E(0):
                continue
            key = min(point_key(image), point_key(-image))
            by_key[key] = None
    return {key: i for i, key in enumerate(sorted(by_key))}


def relation_coefficients(E, witness, h, column_map):
    coefficients = [0]*len(column_map)
    killed = 0
    for x, y in witness:
        P = h*E(E.base_field().from_integer(x), E.base_field().from_integer(y))
        if P == E(0):
            killed += 1
            continue
        pos, neg = point_key(P), point_key(-P)
        key = min(pos, neg)
        coefficients[column_map[key]] += 1 if pos == key else -1
    return coefficients, killed


def run_profile(n, config, source_dir, output):
    screened = json.loads((source_dir/'run-20260909-a'/f'n{n}-results.json').read_text())
    sparse = json.loads((source_dir/'sparse-20260909-a'/f'n{n}-result.json').read_text())
    profile = screened['profile']
    F = GF(2**n, name='a')
    E = EllipticCurve(F, [1, 1, 0, 0, 1])
    r = int(screened['subgroup_order'])
    assert E.cardinality() == 2*r
    P = subgroup_generator(E, r)
    Q = (config['simulation_secret'] % r)*P
    trace = next(row for row in screened['results'] if row['name'] == 'trace-zero-affine')
    dense = max((row for row in screened['results'] if row['name'].startswith('density-selected')),
                key=lambda row: row['actual_subgroup_points'])
    arms = [dict(name='trace', basis=trace['basis'], offset=1, cofactor_multiplier=1),
            dict(name='sparse', basis=sparse['basis'], offset=1, cofactor_multiplier=1),
            dict(name='dense', basis=dense['basis'], offset=1, cofactor_multiplier=1),
            dict(name='full-polynomial', basis=[1 << i for i in range(profile['dimension'])],
                 offset=0, cofactor_multiplier=2)]
    settings = next(row for row in config['profiles'] if row['n'] == n)
    rng = random.Random(config['target_seed']+n)
    targets = []
    while len(targets) < settings['targets']:
        a, b = rng.randrange(r), rng.randrange(r)
        R = a*P+b*Q
        if R != E(0):
            targets.append((a, b, R))
    states = {}
    for arm in arms:
        col = columns(E, r, arm)
        states[arm['name']] = dict(columns=col, rank=RelationRank(r), runs=[],
                                  first_extraction=None, measured_seconds=0.0)
    def report():
        result = dict(scope='toy SAT relation/rank validation; not a cryptanalytic breakthrough',
                      profile=profile, subgroup_order=r, public_P=point_key(P), public_Q=point_key(Q),
                      solver='cadical', solver_watchdog_seconds=settings['solver_watchdog_seconds_per_target'],
                      peak_child_memory_bytes=None,
                      memory_note='Per-child peak memory is unavailable from this subprocess adapter.', arms=[])
        for arm in arms:
            state = states[arm['name']]
            result['arms'].append(arm | dict(
                columns=[list(key) for key in state['columns']], rank=len(state['rank'].rows),
                first_extraction=state['first_extraction'], measured_seconds=state['measured_seconds'],
                runs=state['runs']))
        return result
    nout = output/f'n{n}'
    nout.mkdir()
    (nout/'targets.json').write_text(json.dumps([dict(a=a, b=b, target=point_key(R)) for a,b,R in targets], indent=2)+'\n')
    for index, (a, b, R) in enumerate(targets):
        order = list(arms)
        rng.shuffle(order)
        for arm in order:
            state = states[arm['name']]
            h = arm['cofactor_multiplier']
            begin = time.perf_counter()
            trial = lab.solve_verified(E, [F.from_integer(v) for v in arm['basis']], 3,
                R, r if h == 1 else 2*r, 'cadical', settings['solver_watchdog_seconds_per_target'],
                F.from_integer(arm['offset']), trace_filter=h == 1,
                artifact_dir=nout/arm['name']/str(index))
            trial.update(target_index=index, a=a, b=b, rank_before=len(state['rank'].rows))
            if trial['status'] == 'verified_relation':
                coefficients, killed = relation_coefficients(E, trial['witness'], h, state['columns'])
                # Recheck the cleared equation using actual points, not known logs.
                lhs = E(0)
                for key, position in state['columns'].items():
                    lhs += coefficients[position]*E(F.from_integer(key[0]), F.from_integer(key[1]))
                if lhs != h*R:
                    raise RuntimeError('Cofactor-cleared relation certificate failed')
                d = state['rank'].add(coefficients, h*a, h*b)
                trial.update(coefficients=coefficients, killed_torsion_points=killed)
                if d is not None:
                    if d*P != Q:
                        raise RuntimeError('Extracted logarithm fails public-point verification')
                    trial['verified_extracted_log'] = d
                    if state['first_extraction'] is None:
                        state['first_extraction'] = dict(target_index=index, log=d,
                            measured_seconds=state['measured_seconds']+time.perf_counter()-begin)
            trial['rank_after'] = len(state['rank'].rows)
            trial['measured_wall_seconds'] = time.perf_counter()-begin
            state['measured_seconds'] += trial['measured_wall_seconds']
            state['runs'].append(trial)
            # A progress snapshot is operational, final result is written below.
            (nout/'progress.json').write_text(json.dumps(report(), indent=2)+'\n')
            print(json.dumps(dict(n=n, arm=arm['name'], target=index, status=trial['status'],
                                  rank=trial['rank_after'], seconds=trial['measured_wall_seconds'])), flush=True)
    (output/f'n{n}-result.json').write_text(json.dumps(report(), indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    config_path = args.source/'sat-rank-amendment.json'
    config = json.loads(config_path.read_text())
    args.output.mkdir(parents=True, exist_ok=False)
    for path in [Path(__file__), Path(lab.__file__), Path(__file__).with_name('factor_base_yield_screen.py'), config_path]:
        (args.output/path.name).write_bytes(path.read_bytes())
    manifest = dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest())
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    for profile in config['profiles']:
        run_profile(profile['n'], config, args.source, args.output)


if __name__ == '__main__':
    main()
