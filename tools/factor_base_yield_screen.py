"""Frozen factor-base screen with exact small-group coverage controls.

Run with sage -python. Logarithms are confined to the separately timed oracle;
the factor-space selector uses only actual point membership counts.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import random
import time

from sage.all import GF, EllipticCurve, set_random_seed
from sage.env import SAGE_VERSION
import binary_factor_base as lab


def bit_rank(values):
    pivots = {}
    for value in values:
        while value:
            bit = value.bit_length()-1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
    return len(pivots)


def bit_span(basis, offset):
    values = [offset]
    for b in basis:
        values += [v ^ b for v in values]
    return values


def random_trace_basis(F, k, rng):
    basis = []
    while len(basis) < k:
        value = rng.randrange(1 << F.degree())
        value ^= int(F.from_integer(value).trace())
        if bit_rank(basis + [value]) > len(basis):
            basis.append(value)
    return basis


def cyclic_sumset(scalars, r, m):
    """Exact integer bitset arithmetic, no floating-point convolution."""
    bits = 1
    mask = (1 << r)-1
    for _ in range(m):
        result = 0
        for s in scalars:
            result |= (bits << s) | (bits >> (r-s))
        bits = result & mask
    return bits


def subgroup_generator(E, r):
    F = E.base_field()
    h = E.cardinality() // r
    for i in range(1, F.order()):
        try:
            points = E.lift_x(F.from_integer(i), all=True)
        except ValueError:
            continue
        for Q in points:
            P = h*Q
            if P != E(0) and r*P == E(0):
                return P
    raise RuntimeError('No subgroup generator')


class Membership:
    def __init__(self, E, r):
        self.E, self.r = E, r
        self.cache = {}
        self.queries = self.misses = 0

    def test(self, x):
        self.queries += 1
        if x not in self.cache:
            self.misses += 1
            try:
                points = self.E.lift_x(self.E.base_field().from_integer(x), all=True)
            except ValueError:
                points = []
            self.cache[x] = sum(self.r*P == self.E(0) for P in points)
        return self.cache[x]


def shape(F, E, basis, offset, m, R):
    start = time.perf_counter()
    equations, names = lab.descend(F, [F.from_integer(v) for v in basis], m,
                                  R[0], E.a6(), F.from_integer(offset), E.a2().trace())
    nonlinear = {mon for eq in equations for mon in eq if mon.bit_count() > 1}
    nvars, clauses = lab.cnf(equations, len(names))
    return dict(boolean_variables=len(names), equations=len(equations),
                anf_terms=sum(map(len, equations)), nonlinear_monomials=len(nonlinear),
                cnf_variables=nvars, cnf_clauses=len(clauses),
                encoding_wall_seconds=time.perf_counter()-start)


def screen_profile(profile, protocol, output_dir):
    set_random_seed(protocol['seeds'][0])
    n, k = profile['n'], profile['dimension']
    F = GF(2**n, name='a')
    E = EllipticCurve(F, [1, F.from_integer(profile['a']), 0, 0, F.from_integer(profile['b'])])
    order = int(E.cardinality())
    r = int(max(p for p, _ in E.cardinality().factor()))
    if order != 2*r:
        raise ValueError('This frozen screen requires cofactor two')
    P = subgroup_generator(E, r)
    membership = Membership(E, r)
    generated, dimensions, factors = lab.candidates(F, k, protocol['seeds'][0])
    selected = []
    for name, basis in generated:
        basis = [int(b.to_integer()) for b in basis]
        if name == 'polynomial':
            selected.append(dict(name='polynomial', basis=basis, offset=0))
        elif name == 'trace-zero' or name.startswith('invariant'):
            selected.append(dict(name=name+'-affine', basis=basis, offset=1))
    pools = []
    selection_start = time.perf_counter()
    for seed in protocol['seeds']:
        rng = random.Random(seed)
        pool = []
        for index in range(protocol['pool_size_per_seed']):
            start = time.perf_counter()
            before = membership.queries
            before_misses = membership.misses
            basis = random_trace_basis(F, k, rng)
            values = bit_span(basis, 1)
            count = sum(membership.test(x) for x in values)
            pool.append(dict(index=index, basis=basis, offset=1, subgroup_points=count,
                             wall_seconds=time.perf_counter()-start,
                             membership_queries=membership.queries-before,
                             uncached_membership_queries=membership.misses-before_misses))
        winner = max(pool, key=lambda row: (row['subgroup_points'], -row['index']))
        pools.append(dict(seed=seed, candidates=pool, selected_index=winner['index']))
        for label, item in [('random-first', pool[0]), ('density-selected', winner)]:
            selected.append(dict(name=f'{label}-{seed}', basis=item['basis'], offset=1))
    selection_seconds = time.perf_counter()-selection_start
    # Freeze the selection before constructing any logarithm oracle.
    selection_record = dict(profile=profile, pools=pools, selected=selected,
                            selection_wall_seconds=selection_seconds,
                            membership_queries=membership.queries,
                            uncached_membership_queries=membership.misses,
                            invariant_dimensions=dimensions, frobenius_factors=factors)
    (output_dir / f'n{n}-selection.json').write_text(json.dumps(selection_record, indent=2)+'\n')

    oracle_start = time.perf_counter()
    x_logs = {}
    Q = P
    for scalar in range(1, (r+1)//2):
        x = int(Q[0].to_integer())
        if x in x_logs:
            raise RuntimeError('Subgroup oracle repeated an x-coordinate before sign quotient')
        x_logs[x] = scalar
        Q += P
    if len(x_logs) != (r-1)//2:
        raise RuntimeError('Incorrect oracle cardinality')
    # Check membership-table agreement and sampled independent scalar multiplication.
    for x, count in membership.cache.items():
        if (count == 2) != (x in x_logs):
            raise RuntimeError('Subgroup membership and group enumeration disagree')
    rng = random.Random(protocol['seeds'][-1])
    for x in rng.sample(list(x_logs), min(64, len(x_logs))):
        if int((x_logs[x]*P)[0].to_integer()) != x:
            raise RuntimeError('Oracle disagrees with scalar multiplication')
    oracle_seconds = time.perf_counter()-oracle_start
    results = []
    for i, candidate in enumerate(selected):
        values = bit_span(candidate['basis'], candidate['offset'])
        representative_logs = [x_logs[x] for x in values if x in x_logs]
        scalars = [v for s in representative_logs for v in (s, r-s)]
        coverage_mask = cyclic_sumset(scalars, r, protocol['summands'])
        coverage = (coverage_mask & ~1).bit_count()
        rng = random.Random(protocol['seeds'][0] + i + n*1000)
        matched_coverages = []
        for _ in range(protocol['pool_size_per_seed']):
            sample = rng.sample(range(1, (r+1)//2), len(representative_logs))
            null_scalars = [v for s in sample for v in (s, r-s)]
            matched_coverages.append((cyclic_sumset(null_scalars, r, protocol['summands']) & ~1).bit_count())
        sizes = shape(F, E, candidate['basis'], candidate['offset'], protocol['summands'], P)
        results.append(candidate | dict(
            actual_subgroup_points=len(scalars), sign_folded_columns=len(representative_logs),
            nonzero_target_coverage=coverage, coverage_fraction=coverage/(r-1),
            coverage_per_column=coverage/len(representative_logs) if representative_logs else None,
            coverage_mask_hex=hex(coverage_mask),
            matched_cardinality_random_coverages=matched_coverages, encoding=sizes))
    result = dict(scope='exhaustive small-field screening, not an index-calculus attack result',
                  profile=profile, field_modulus=str(F.modulus()), group_order=order, subgroup_order=r,
                  subgroup_generator=[int(P[0].to_integer()), int(P[1].to_integer())],
                  oracle_wall_seconds=oracle_seconds, oracle_sign_quotient_entries=len(x_logs),
                  selection_file=f'n{n}-selection.json', results=results)
    (output_dir / f'n{n}-results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(n=n, subgroup_order=r, oracle_seconds=oracle_seconds,
                         candidates=[dict(name=row['name'], points=row['actual_subgroup_points'],
                                          coverage=row['nonzero_target_coverage'],
                                          null_mean=sum(row['matched_cardinality_random_coverages'])/len(row['matched_cardinality_random_coverages']))
                                     for row in results])), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text())
    args.output.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).read_bytes()
    (args.output/'factor_base_yield_screen.py').write_bytes(source)
    (args.output/'binary_factor_base.py').write_bytes(Path(lab.__file__).read_bytes())
    (args.output/'protocol.json').write_bytes(args.protocol.read_bytes())
    manifest = dict(sage=SAGE_VERSION, python=platform.python_version(), platform=platform.platform(),
                    source_sha256=hashlib.sha256(source).hexdigest(),
                    protocol_sha256=hashlib.sha256(args.protocol.read_bytes()).hexdigest(),
                    oracle_uses_known_logs=True, selector_uses_known_logs=False)
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    for profile in protocol['profiles']:
        screen_profile(profile, protocol, args.output)


if __name__ == '__main__':
    main()
