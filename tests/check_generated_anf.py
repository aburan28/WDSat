"""Cross-check generated ANF against its retained CaDiCaL CNF result."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def valid(text, bits):
    for row in text.splitlines()[1:]:
        tokens = row.split()
        if not tokens:
            continue
        if tokens[0] != 'x':
            lit = int(tokens[0])
            if bool(bits[abs(lit)-1]) != (lit > 0):
                return False
            continue
        parity = 0
        pos = 1
        while tokens[pos] != '0':
            token = tokens[pos]
            pos += 1
            if token == 'T':
                parity ^= 1
            elif token.startswith('.'):
                degree = int(token[1:])
                parity ^= int(all(bits[int(v)-1] for v in tokens[pos:pos+degree]))
                pos += degree
            else:
                parity ^= bits[int(token)-1]
        if parity != 1:
            return False
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('artifacts', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Choose a new output path')
    solver = Path(__file__).resolve().parents[1] / 'wdsat_solver'
    checks = []
    for path in sorted(args.artifacts.rglob('input.anf')):
        expected = json.loads(path.with_name('attempt-0000.json').read_text())['returncode']
        text = path.read_text()
        nvars = int(text.splitlines()[0].split()[2])
        for flags in [[], ['-x'], ['-g', '1'], ['-x', '-g', '1']]:
            result = subprocess.run([str(solver), '-i', str(path), *flags],
                                    capture_output=True, text=True, timeout=10)
            assert result.returncode == 0, (path, flags, result.stderr)
            unsat = 'UNSAT' in result.stdout
            assert unsat == (expected == 20), (path, flags, result.stdout)
            if not unsat:
                model = result.stdout.strip().splitlines()[-2]
                assert len(model) == nvars and set(model) <= {'0', '1'}, (path, model)
                assert valid(text, list(map(int, model))), (path, flags, model)
            checks.append(dict(input=str(path), flags=flags, status='UNSAT' if unsat else 'validated_model',
                               stdout=result.stdout, stderr=result.stderr, returncode=result.returncode,
                               input_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    assert checks, 'No input files found'
    args.output.write_text(json.dumps(dict(scope='solver cross-check; not an ECDLP claim',
        solver_sha256=hashlib.sha256(solver.read_bytes()).hexdigest(), checks=checks), indent=2) + '\n')
    print(f'{len(checks)} WDSat checks passed against retained CaDiCaL statuses and direct ANF evaluation')


if __name__ == '__main__':
    main()
