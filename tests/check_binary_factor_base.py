"""Run: sage -python tests/check_binary_factor_base.py"""
import itertools
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
import binary_factor_base as lab
from sage.all import GF, EllipticCurve


class FactorBaseTests(unittest.TestCase):
    def test_invariant_dimensions_and_closure(self):
        for n, k, expected in [(7, 3, [0, 1, 3, 4, 6, 7]),
                               (17, 3, [0, 1, 8, 9, 16, 17])]:
            F = GF(2**n, name='a')
            bases, dimensions, _ = lab.candidates(F, k, 1)
            self.assertEqual(dimensions, expected)
            for name, basis in bases:
                self.assertEqual(lab.basis_rank(basis), k)
                if name.startswith('invariant'):
                    values = lab.span(basis)
                    self.assertTrue(all(x*x in values for x in values))

    def test_descent_against_direct_field_arithmetic(self):
        F = GF(8, name='a')
        basis = [F(1), F.gen()]
        for m in [2, 3]:
            for target in F:
                equations, names = lab.descend(F, basis, m, target, F(1))
                def s3(x, y, z):
                    return (x*y + x*z + y*z)**2 + x*y*z + 1
                for bits in itertools.product((0, 1), repeat=len(names)):
                    xs = [sum((basis[j]*bits[i*2+j] for j in range(2)), F(0)) for i in range(m)]
                    if m == 2:
                        expected = s3(xs[0], xs[1], target) == 0
                    else:
                        u = sum((F.gen()**j*bits[m*2+j] for j in range(3)), F(0))
                        expected = s3(xs[0], xs[1], u) == 0 and s3(u, xs[2], target) == 0
                    self.assertEqual(lab.evaluate(equations, bits), expected)

    def test_cnf_projection_equivalence(self):
        # Shared products, constant-only contradiction, tautology, and mixed XOR.
        systems = [[{1, 2, 3}, {3, 4, 0}], [{0}], [set()], [{1, 2, 4, 0}]]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'input.cnf'
            for equations in systems:
                nvars, clauses = lab.cnf(equations, 3)
                for bits in itertools.product((0, 1), repeat=3):
                    units = [[i+1 if b else -(i+1)] for i, b in enumerate(bits)]
                    path.write_text(lab.cnf_text(nvars, clauses + units))
                    result = subprocess.run(['cadical', str(path)], capture_output=True, text=True, timeout=5)
                    self.assertIn(result.returncode, [10, 20])
                    self.assertEqual(result.returncode == 10, lab.evaluate(equations, bits))

    def test_wdsat_export_projection_equivalence(self):
        solver = Path(__file__).resolve().parents[1] / 'wdsat_solver'
        systems = [[{1, 2, 3}, {3, 4, 0}], [{0}], [set()], [{1, 2, 4, 0}]]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'input.anf'
            for equations in systems:
                for bits in itertools.product((0, 1), repeat=3):
                    # Unit constraints are ANFs of the form x + value = 0.
                    units = [{1 << i, 0} if b else {1 << i} for i, b in enumerate(bits)]
                    path.write_text(lab.anf_text(equations + units, 3))
                    result = subprocess.run([str(solver), '-i', str(path), '-x'],
                                            capture_output=True, text=True, timeout=5)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertEqual('UNSAT' not in result.stdout, lab.evaluate(equations, bits))

    def test_verified_curve_witness(self):
        F = GF(128, name='a')
        E = EllipticCurve(F, [1, 1, 0, 0, 1])
        basis = [F(1), F.gen(), F.gen()**2]
        r = max(p for p, _ in E.cardinality().factor())
        points = []
        for x in lab.span(basis):
            try:
                points.extend(P for P in E.lift_x(x, all=True) if r*P == E(0))
            except ValueError:
                pass
        R = next(P+Q for P in points for Q in points if P+Q != E(0))
        result = lab.solve_verified(E, basis, 2, R, r, 'cadical', 10)
        self.assertEqual(result['status'], 'verified_relation')
        recovered = [E(F.from_integer(x), F.from_integer(y)) for x, y in result['witness']]
        self.assertEqual(sum(recovered, E(0)), R)
        self.assertTrue(all(r*P == E(0) for P in recovered))

    def test_affine_factor_base(self):
        F = GF(128, name='a')
        E = EllipticCurve(F, [1, 1, 0, 0, 1])
        bases, _, _ = lab.candidates(F, 3, 1)
        basis = next(b for name, b in bases if name.startswith('invariant'))
        values = lab.span(basis, F(1))
        self.assertEqual(len(values), 8)
        self.assertTrue(all(x*x in values for x in values))
        self.assertTrue(all(x.trace() == 1 for x in values))
        r = max(p for p, _ in E.cardinality().factor())
        points = []
        for x in values:
            try:
                points.extend(P for P in E.lift_x(x, all=True) if r*P == E(0))
            except ValueError:
                pass
        self.assertTrue(points)
        R = next(P+Q for P in points for Q in points if P+Q != E(0))
        result = lab.solve_verified(E, basis, 2, R, r, 'cadical', 10, F(1))
        self.assertEqual(result['status'], 'verified_relation')
        recovered = [E(F.from_integer(x), F.from_integer(y)) for x, y in result['witness']]
        self.assertEqual(sum(recovered, E(0)), R)
        self.assertTrue(all(P[0] in values for P in recovered))

    def test_trace_filter_on_general_binary_curves(self):
        F = GF(32, name='a')
        for A, B in [(F(0), F(1)), (F(1), F(1)), (F.gen(), F.gen())]:
            E = EllipticCurve(F, [1, A, 0, 0, B])
            for P in E:
                if P != E(0) and P.order() % 2:
                    self.assertEqual(P[0].trace(), A.trace())
            basis = [F(1), F.gen()]
            equations, names = lab.descend(F, basis, 2, F.gen(), B, F(1), A.trace())
            unfiltered, _ = lab.descend(F, basis, 2, F.gen(), B, F(1))
            for bits in itertools.product((0, 1), repeat=len(names)):
                xs = [F(1) + sum((basis[j]*bits[i*2+j] for j in range(2)), F(0))
                      for i in range(2)]
                expected = lab.evaluate(unfiltered, bits) and all(x.trace() == A.trace() for x in xs)
                self.assertEqual(lab.evaluate(equations, bits), expected)


if __name__ == '__main__':
    unittest.main()
