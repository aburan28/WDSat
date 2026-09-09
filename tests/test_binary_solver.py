"""Independent truth-table oracle for small Boolean systems (not ECDLP evidence)."""
import itertools
from pathlib import Path
import random
import re
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def run(text, flags=()):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / 'input.anf'
        path.write_text(text)
        return subprocess.run([str(ROOT / 'wdsat_solver'), '-i', str(path), *flags],
                              capture_output=True, text=True, timeout=10)


def assignment(output, n):
    # Last line is a decision counter, which can also resemble a bit string.
    lines = output.strip().splitlines()
    if len(lines) < 2 or not re.fullmatch('[01]{%d}' % n, lines[-2]):
        raise AssertionError('Missing complete Boolean model: ' + output)
    return tuple(map(int, lines[-2]))


class BinarySolverTests(unittest.TestCase):
    def test_initial_conflict(self):
        for flags in [(), ('-x',)]:
            result = run('p cnf 1 2\n1 0\n-1 0\n', flags)
            self.assertEqual(result.returncode, 0)
            self.assertIn('UNSAT', result.stdout)

    def test_partial_order_completed(self):
        for flags in [('-g', '1'), ('-x', '-g', '1')]:
            result = run('p cnf 3 1\nx 2 3 0\n', flags)
            self.assertEqual(result.returncode, 0, result.stderr)
            bits = assignment(result.stdout, 3)
            self.assertEqual(bits[1] ^ bits[2], 1)

    def test_mixed_anf_and_general_cnf(self):
        text = 'p cnf 3 3\nx 1 2 T 0\n1 2 3 0\n-1 -2 0\n'
        for flags in [(), ('-x',)]:
            result = run(text, flags)
            self.assertEqual(result.returncode, 0, result.stderr)
            bits = assignment(result.stdout, 3)
            self.assertEqual(bits[0] ^ bits[1], 0)
            self.assertTrue(bits[0] or bits[1] or bits[2])
            self.assertFalse(bits[0] and bits[1])

    def test_symmetry_completes_auxiliary_variables(self):
        # Symmetric in x1,x2, with unresolved auxiliary XOR in x3,x4.
        for flags in [(), ('-x',)]:
            result = run('p cnf 4 1\nx 3 4 0\n',
                         (*flags, '-b', '-l', '1', '-m', '2'))
            self.assertEqual(result.returncode, 0, result.stderr)
            bits = assignment(result.stdout, 4)
            self.assertLessEqual(bits[0], bits[1])
            self.assertEqual(bits[2] ^ bits[3], 1)

    def test_invalid_order(self):
        for order in ['0', '-1', '4', '1,1', '1,', ',1', '1,,2', 'abc']:
            result = run('p cnf 3 1\nx 2 3 0\n', ('-g', order))
            self.assertNotEqual(result.returncode, 0, order)

    def test_incompatible_symmetry_order(self):
        result = run('p cnf 3 1\nx 2 3 0\n',
                     ('-g', '2,1', '-b', '-l', '1', '-m', '2'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('canonical', result.stderr)

    def test_random_anf_against_truth_tables(self):
        rng = random.Random(20260909)
        for trial in range(100):
            n = rng.randrange(2, 7)
            monomials = [p for degree in range(1, min(n, 3) + 1)
                         for p in itertools.combinations(range(1, n + 1), degree)]
            equations = [(rng.sample(monomials, rng.randrange(1, min(8, len(monomials)) + 1)),
                          rng.randrange(2)) for _ in range(rng.randrange(1, 9))]
            def valid(bits):
                return all((sum(all(bits[i - 1] for i in mon) for mon in mons) % 2) == rhs
                           for mons, rhs in equations)
            models = [bits for bits in itertools.product((0, 1), repeat=n) if valid(bits)]
            rows = []
            for mons, rhs in equations:
                terms = [str(mon[0]) if len(mon) == 1 else
                         '.' + str(len(mon)) + ' ' + ' '.join(map(str, mon)) for mon in mons]
                # WDSat XOR clauses evaluate to true, not zero.
                rows.append('x ' + ' '.join(terms + (['T'] if rhs == 0 else [])) + ' 0')
            text = f'p cnf {n} {len(rows)}\n' + '\n'.join(rows) + '\n'
            for flags in [(), ('-x',), ('-g', str(n)), ('-x', '-g', str(n))]:
                result = run(text, flags)
                self.assertEqual(result.returncode, 0, (trial, flags, result.stderr))
                if models:
                    self.assertNotIn('UNSAT', result.stdout, (trial, flags, text))
                    self.assertTrue(valid(assignment(result.stdout, n)), (trial, flags, text))
                else:
                    self.assertIn('UNSAT', result.stdout, (trial, flags, text))


if __name__ == '__main__':
    unittest.main()
