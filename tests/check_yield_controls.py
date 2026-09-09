"""Independent arithmetic tests for cofactor handling and relation-rank recovery."""
import itertools
from pathlib import Path
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from sage.all import GF, EllipticCurve
from factor_base_cofactor_control import profile_control
from factor_base_sat_rank import RelationRank, columns, relation_coefficients, point_key
from factor_base_yield_screen import bit_span, subgroup_generator
from factor_base_frobenius_schedule import add_mask, score_sets


class YieldControlTests(unittest.TestCase):
    def test_distinct_cyclic_sumsets(self):
        for r in [7,11,71]:
            sets=[[1,r-1],[2,r-2],[3,r-3]]
            mask,coverage=score_sets(sets,r)
            expected={sum(v)%r for v in itertools.product(*sets)}
            self.assertEqual(mask,sum(1 << v for v in expected))
            self.assertEqual(coverage,len(expected-{0}))
    def test_rank_recovery(self):
        rng = random.Random(345)
        r, secret = 991, 42
        logs = [14, 571, 52, 608]
        rank = RelationRank(r)
        extracted = []
        for _ in range(12):
            row = [rng.randrange(r) for _ in logs]
            b = rng.randrange(1, r)
            a = (sum(x*y for x, y in zip(row, logs))-b*secret) % r
            d = rank.add(row, a, b)
            if d is not None:
                extracted.append(d)
        self.assertEqual(len(rank.rows), len(logs))
        self.assertTrue(extracted)
        self.assertTrue(all(d == secret for d in extracted))

    def test_inconsistent_dependency(self):
        rank = RelationRank(71)
        with self.assertRaises(RuntimeError):
            rank.add([0, 0], 1, 0)

    def test_cofactor_coverage_against_group_addition(self):
        profile = dict(n=7, a=1, b=1, dimension=3)
        result = profile_control(profile, 3)
        F = GF(128, name='a')
        E = EllipticCurve(F, [1, 1, 0, 0, 1])
        points = []
        for x in bit_span([1, 2, 4], 0):
            try:
                points.extend(P for P in E.lift_x(F.from_integer(x), all=True) if 2*P != E(0))
            except ValueError:
                pass
        sums = {P+Q+R for P,Q,R in itertools.product(points, repeat=3)}
        subgroup_targets = {R for R in sums if R != E(0) and 71*R == E(0)}
        cleared_targets = {2*R for R in sums}-{E(0)}
        self.assertEqual(len(subgroup_targets), result['full_curve_single_subgroup_target_coverage'])
        self.assertEqual(len(cleared_targets), result['full_curve_two_preimage_target_coverage'])

    def test_cleared_coefficients_with_torsion(self):
        F = GF(128, name='a')
        E = EllipticCurve(F, [1, 1, 0, 0, 1])
        candidate = dict(basis=[1, 2, 4], offset=0, cofactor_multiplier=2)
        column_map = columns(E, 71, candidate)
        points = []
        for x in bit_span(candidate['basis'], 0):
            try:
                points.extend(E.lift_x(F.from_integer(x), all=True))
            except ValueError:
                pass
        T = E.lift_x(F(0))
        point = next(P for P in points if 2*P != E(0))
        witness = [point_key(T), point_key(point), point_key(point)]
        coeff, killed = relation_coefficients(E, witness, 2, column_map)
        total = E(0)
        for key, index in column_map.items():
            total += coeff[index]*E(F.from_integer(key[0]), F.from_integer(key[1]))
        self.assertEqual(killed, 1)
        self.assertEqual(total, 2*(T+point+point))


if __name__ == '__main__':
    unittest.main()
