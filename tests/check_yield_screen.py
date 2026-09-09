"""Run with sage -python; validates the exact coverage oracle and selector."""
import contextlib
import io
import itertools
import json
from pathlib import Path
import random
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
import factor_base_yield_screen as screen
from sage.all import GF, Matrix


class YieldScreenTests(unittest.TestCase):
    def test_cyclic_sumset(self):
        rng = random.Random(12871)
        for r in [3, 7, 11, 71]:
            for _ in range(20):
                values = rng.sample(range(r), rng.randrange(r+1))
                for m in [1, 2, 3]:
                    expected = {sum(items) % r for items in itertools.product(values, repeat=m)}
                    actual = screen.cyclic_sumset(values, r, m)
                    self.assertEqual(actual, sum(1 << v for v in expected))

    def test_basis_and_trace(self):
        F = GF(2**11, name='a')
        rng = random.Random(871)
        for _ in range(20):
            basis = screen.random_trace_basis(F, 4, rng)
            self.assertEqual(screen.bit_rank(basis), 4)
            self.assertEqual(Matrix(GF(2), [F.from_integer(b)._vector_() for b in basis]).rank(), 4)
            values = screen.bit_span(basis, 1)
            self.assertEqual(len(set(values)), 16)
            self.assertTrue(all(F.from_integer(v).trace() == 1 for v in values))

    def test_small_profile_integrity(self):
        protocol = dict(seeds=[821, 822], pool_size_per_seed=3, summands=3)
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            root = Path(directory)
            screen.screen_profile(dict(n=7, a=1, b=1, dimension=3), protocol, root)
            selection = json.loads((root/'n7-selection.json').read_text())
            result = json.loads((root/'n7-results.json').read_text())
            self.assertEqual(result['subgroup_order'], 71)
            for pool in selection['pools']:
                winner = max(pool['candidates'], key=lambda row: (row['subgroup_points'], -row['index']))
                self.assertEqual(winner['index'], pool['selected_index'])
            for row in result['results']:
                self.assertEqual(row['nonzero_target_coverage'], (int(row['coverage_mask_hex'], 16)&~1).bit_count())
                self.assertEqual(row['actual_subgroup_points'], 2*row['sign_folded_columns'])


if __name__ == '__main__':
    unittest.main()
