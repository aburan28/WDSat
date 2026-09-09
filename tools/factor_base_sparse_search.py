"""Enumerate sparse trace-compatible affine directions without target feedback."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import random
import time

from sage.all import GF, EllipticCurve
from factor_base_yield_screen import bit_span, shape, subgroup_generator, cyclic_sumset


def search(profile, m, output):
    n, k = profile['n'], profile['dimension']
    F = GF(2**n, name='a')
    E = EllipticCurve(F, [1, 1, 0, 0, 1])
    r = int(max(p for p, _ in E.cardinality().factor()))
    if profile['a'] != 1 or profile['b'] != 1 or E.cardinality() != 2*r or n % 2 == 0:
        raise ValueError('This screen is only for A=B=1, odd n and cofactor two')
    directions = {i: int((F.gen()**i + F((F.gen()**i).trace())).to_integer()) for i in range(1, n)}
    cache = {}
    scores = []
    start = time.perf_counter()
    queries = 0
    for exponents in itertools.combinations(range(1, n), k):
        basis = [directions[i] for i in exponents]
        count = 0
        for x in bit_span(basis, 1):
            queries += 1
            if x not in cache:
                value = F.from_integer(x)
                assert value.trace() == 1
                cache[x] = int(value != 0 and (1/value).trace() == 0)
            count += cache[x]
        scores.append(dict(exponents=list(exponents), valid_x=count))
    winner = min(scores, key=lambda row: (-row['valid_x'], sum(row['exponents']), row['exponents']))
    basis = [directions[i] for i in winner['exponents']]
    selection_seconds = time.perf_counter()-start
    selection = dict(profile=profile, directions=directions, scores=scores, winner=winner,
                     basis=basis, offset=1, membership_queries=queries,
                     unique_membership_queries=len(cache), selection_wall_seconds=selection_seconds)
    # Selection is durable before any group/log oracle or equation-size comparison.
    (output/f'n{n}-selection.json').write_text(json.dumps(selection, indent=2)+'\n')
    wanted = set(bit_span(basis, 1))
    valid = []
    for x in wanted:
        try:
            points = E.lift_x(F.from_integer(x), all=True)
        except ValueError:
            points = []
        subgroup_points = [P for P in points if r*P == E(0)]
        if bool(subgroup_points) != bool(cache[x]) or len(subgroup_points) not in (0, 2):
            raise RuntimeError('Trace membership and curve/subgroup checks disagree')
        if subgroup_points:
            valid.append(x)
    assert len(valid) == winner['valid_x']
    P = subgroup_generator(E, r)
    logs = []
    Q = P
    oracle_start = time.perf_counter()
    for j in range(1, (r+1)//2):
        if int(Q[0].to_integer()) in wanted:
            logs.extend([j, r-j])
        Q += P
    assert len(logs) == 2*len(valid)
    covered = (cyclic_sumset(logs, r, m)&~1).bit_count()
    rng = random.Random(202609091+n)
    nulls = []
    for _ in range(32):
        sample = rng.sample(range(1, (r+1)//2), len(valid))
        nulls.append((cyclic_sumset([s for j in sample for s in (j, r-j)], r, m)&~1).bit_count())
    result = dict(profile=profile, subgroup_order=r, field_modulus=str(F.modulus()),
                  basis=basis, offset=1, exponents=winner['exponents'],
                  actual_subgroup_points=len(logs), sign_folded_columns=len(valid),
                  exact_nonzero_target_coverage=covered, matched_cardinality_random_coverages=nulls,
                  oracle_wall_seconds=time.perf_counter()-oracle_start,
                  encoding=shape(F, E, basis, 1, m, P))
    (output/f'n{n}-result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text())
    args.output.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).read_bytes()
    (args.output/'factor_base_sparse_search.py').write_bytes(source)
    (args.output/'manifest.json').write_text(json.dumps(dict(
        scope='density selection and exact toy coverage; no breakthrough claim',
        source_sha256=hashlib.sha256(source).hexdigest(),
        protocol_sha256=hashlib.sha256(args.protocol.read_bytes()).hexdigest()), indent=2)+'\n')
    for profile in protocol['profiles']:
        search(profile, protocol['summands'], args.output)


if __name__ == '__main__':
    main()
