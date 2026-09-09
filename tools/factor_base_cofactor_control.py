"""Exact cofactor-aware full-curve control for the frozen yield screen."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from sage.all import GF, EllipticCurve
from factor_base_yield_screen import subgroup_generator, cyclic_sumset, bit_span


def profile_control(profile, m):
    start = time.perf_counter()
    n, k = profile['n'], profile['dimension']
    F = GF(2**n, name='a')
    E = EllipticCurve(F, [1, F.from_integer(profile['a']), 0, 0, F.from_integer(profile['b'])])
    N = int(E.cardinality())
    r = int(max(p for p, _ in E.cardinality().factor()))
    if N != 2*r:
        raise ValueError('Control requires cofactor two')
    P = subgroup_generator(E, r)
    T = E.lift_x(F(0))
    assert T != E(0) and 2*T == E(0)
    G = P+T
    assert r*G == T and N*G == E(0) and 2*G == 2*P
    wanted_x = set(bit_span([1 << i for i in range(k)], 0))
    scalars = []
    subgroup_scalars = []
    representative_points = []
    Q = G
    for scalar in range(1, r+1):
        x = int(Q[0].to_integer())
        if x in wanted_x and scalar != r:  # T maps to zero: no unknown log column.
            assert 2*Q != E(0)
            scalars.extend([scalar, N-scalar])
            representative_points.append([x, int(Q[1].to_integer()), scalar])
            if scalar % 2 == 0:
                subgroup_scalars.extend([scalar, N-scalar])
        Q += G
    full_mask = cyclic_sumset(scalars, N, m)
    even_mask = ((1 << N)-1)//3
    subgroup_target_coverage = (full_mask & even_mask & ~1).bit_count()
    projected = sorted({2*s % r for s in scalars})
    projected_columns = {min(s, r-s) for s in projected if s}
    projected_coverage = (cyclic_sumset(projected, r, m) & ~1).bit_count()
    # Independently project the full sumset. Doubling kills its torsion component.
    projected_full_mask = 0
    remaining = full_mask
    while remaining:
        lowest = remaining & -remaining
        s = lowest.bit_length()-1
        projected_full_mask |= 1 << (2*s % r)
        remaining ^= lowest
    assert projected_full_mask == cyclic_sumset(projected, r, m)
    sub_only_coverage = (cyclic_sumset(subgroup_scalars, N, m) & even_mask & ~1).bit_count()
    return dict(profile=profile, field_modulus=str(F.modulus()), group_order=N, subgroup_order=r,
                curve_generator=[int(G[0].to_integer()), int(G[1].to_integer())],
                subgroup_generator=[int(P[0].to_integer()), int(P[1].to_integer())],
                points=representative_points, full_curve_factor_points=len(scalars),
                subgroup_only_factor_points=len(subgroup_scalars),
                original_sign_columns=len(scalars)//2, cofactor_image_sign_columns=len(projected_columns),
                subgroup_only_target_coverage=sub_only_coverage,
                full_curve_single_subgroup_target_coverage=subgroup_target_coverage,
                full_curve_two_preimage_target_coverage=projected_coverage,
                two_preimage_query_multiplier_upper_bound=2,
                wall_seconds=time.perf_counter()-start)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Choose a new output path')
    protocol = json.loads(args.protocol.read_text())
    result = dict(scope='small-field cofactor-control oracle, not an attack performance measurement',
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), profiles=[])
    for profile in protocol['profiles']:
        row = profile_control(profile, protocol['summands'])
        result['profiles'].append(row)
        print(json.dumps({k: v for k, v in row.items() if k != 'points'}), flush=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
