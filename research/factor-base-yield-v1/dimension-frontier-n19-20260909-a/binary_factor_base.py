"""Small-field factor-base/SAT laboratory. Run with ``sage -python``.

No Groebner basis calls. Boolean monomials are bit masks; coefficients are
in GF(2^n). Squaring transforms coefficients only, using b^2=b immediately.
This is a diagnostic implementation, not a full index-calculus attack.
"""
import argparse
from contextlib import nullcontext
import hashlib
import itertools
import json
from pathlib import Path
import random
import platform
import subprocess
import tempfile
import time

from sage.all import GF, EllipticCurve, Matrix, PolynomialRing, set_random_seed
from sage.env import SAGE_VERSION


def span(basis, offset=0):
    values = [basis[0].parent()(offset)]
    for b in basis:
        values += [x + b for x in values]
    return values


def basis_rank(basis):
    return Matrix(GF(2), [b._vector_() for b in basis]).rank()


def candidates(F, dimension, seed):
    """Polynomial/random controls and every invariant subspace of this dimension.

    Odd extension degree only: T^n-1 is square-free, and the normal-basis
    Frobenius module is cyclic. Each invariant subspace is a sum of distinct
    irreducible-factor kernels. No assumption that every dimension exists.
    """
    n = F.degree()
    if n % 2 == 0 or not 1 <= dimension <= n:
        raise ValueError('Use odd n and 1 <= dimension <= n')
    a = F.gen()
    result = [('polynomial', [a**i for i in range(dimension)])]
    rng = random.Random(seed)
    random_basis = []
    while len(random_basis) < dimension:
        b = F.from_integer(rng.randrange(F.order()))
        if basis_rank(random_basis + [b]) > len(random_basis):
            random_basis.append(b)
    result.append(('random', random_basis))
    if dimension < n:
        # For odd n, Tr(1)=1. These vectors span ker(Tr), so 1+V
        # supplies a same-dimension trace-one control for affine candidates.
        trace_basis = []
        for i in range(1, n):
            b = a**i + F((a**i).trace())
            if basis_rank(trace_basis + [b]) > len(trace_basis):
                trace_basis.append(b)
            if len(trace_basis) == dimension:
                break
        result.append(('trace-zero', trace_basis))
    ring = PolynomialRing(GF(2), 'T')
    T = ring.gen()
    factors = [f for f, _ in (T**n - 1).factor()]
    frob = Matrix(GF(2), [((a**i)**2)._vector_() for i in range(n)]).transpose()
    kernels = [f(frob).right_kernel().basis() for f in factors]
    dimensions = {0}
    for f in factors:
        dimensions |= {d + int(f.degree()) for d in dimensions}
    for selection in itertools.product((False, True), repeat=len(factors)):
        if sum(f.degree() for f, yes in zip(factors, selection) if yes) != dimension:
            continue
        basis = [F(list(v)) for kernel, yes in zip(kernels, selection) if yes for v in kernel]
        result.append(('invariant-' + ''.join(map(lambda x: str(int(x)), selection)), basis))
    return result, sorted(dimensions), [str(f) for f in factors]


def add(*polys):
    result = {}
    for p in polys:
        for mon, c in p.items():
            result[mon] = result.get(mon, 0) + c
            if result[mon] == 0:
                del result[mon]
    return result


def mul(p, q):
    result = {}
    for i, a in p.items():
        for j, b in q.items():
            mon = i | j
            result[mon] = result.get(mon, 0) + a*b
            if result[mon] == 0:
                del result[mon]
    return result


def s3(x, y, z, B):
    quadratic = add(mul(x, y), mul(x, z), mul(y, z))
    return add({mon: c*c for mon, c in quadratic.items()}, mul(mul(x, y), z), {0: B})


def descend_multi(F, bases, target_x, B, offsets=None, subgroup_trace=None):
    """Return coordinate ANFs equal to zero and an explicit variable map.

    Factor coordinates first, full-field auxiliary coordinates last. The
    chained model's auxiliaries need search: the symmetric-S4 propagation
    theorem from the WDSat paper does not apply to this representation.
    """
    m = len(bases)
    if m not in (2, 3):
        raise ValueError('This diagnostic supports m=2 or m=3')
    if offsets is None:
        offsets = [F(0)]*m
    if len(offsets) != m:
        raise ValueError('One affine offset is required per factor space')
    for basis in bases:
        if not basis or basis_rank(basis) != len(basis):
            raise ValueError('Factor bases must be nonempty and independent')
    names = []
    starts = []
    def variable(name, support):
        p = {}
        for j, coefficient in enumerate(support):
            p[1 << len(names)] = coefficient
            names.append(f'{name}_{j}')
        return p
    xs = []
    for i, (basis, offset) in enumerate(zip(bases, offsets)):
        starts.append(len(names))
        xs.append(add(variable(f'x{i+1}', basis), {0: F(offset)}))
    if m == 2:
        polys = [s3(xs[0], xs[1], {0: target_x}, B)]
    else:
        u = variable('u1', [F.gen()**i for i in range(F.degree())])
        polys = [s3(xs[0], xs[1], u, B), s3(u, xs[2], {0: target_x}, B)]
    equations = []
    for p in polys:
        rows = [set() for _ in range(F.degree())]
        for mon, coefficient in p.items():
            for j, bit in enumerate(coefficient._vector_()):
                if bit:
                    rows[j].add(mon)
        equations.extend(row for row in rows if row)
    if subgroup_trace is not None:
        # Every odd-order subgroup point is a double. In this curve model,
        # x(2P)=lambda^2+lambda+A, hence Tr(x)=Tr(A). This is necessary,
        # not a replacement for lifting or the final subgroup certificate.
        for i in range(m):
            row = {1 << (starts[i]+j) for j, b in enumerate(bases[i]) if b.trace()}
            if F(offsets[i]).trace() != subgroup_trace:
                row.add(0)
            if row:
                equations.append(row)
    return equations, names


def descend(F, basis, m, target_x, B, offset=0, subgroup_trace=None):
    return descend_multi(F, [basis]*m, target_x, B, [offset]*m, subgroup_trace)


def mon_ids(mon):
    return [i + 1 for i in range(mon.bit_length()) if mon & (1 << i)]


def anf_text(equations, nvars):
    # WDSat equations have odd parity. Constants alone need special handling
    # because the legacy parser assumes each XOR row has at least one variable.
    rows = []
    for eq in equations:
        if eq == {0}:
            rows.extend(['1 0', '-1 0'])
            continue
        if not eq:
            continue
        terms = []
        for mon in sorted(eq - {0}):
            ids = mon_ids(mon)
            terms.append(('' if len(ids) == 1 else f'.{len(ids)} ') + ' '.join(map(str, ids)))
        if 0 not in eq:
            terms.append('T')
        rows.append('x ' + ' '.join(terms) + ' 0')
    return f'p cnf {nvars} {len(rows)}\n' + '\n'.join(rows) + '\n'


def cnf(equations, nvars):
    """Tseitin CNF with shared AND gates and linear-size XOR chains."""
    clauses = []
    products = {}
    def fresh():
        nonlocal nvars
        nvars += 1
        return nvars
    def atom(mon):
        ids = mon_ids(mon)
        if len(ids) == 1:
            return ids[0]
        if mon not in products:
            z = fresh()
            clauses.extend([[-z, x] for x in ids])
            clauses.append([z] + [-x for x in ids])
            products[mon] = z
        return products[mon]
    for eq in equations:
        terms = [atom(mon) for mon in sorted(eq - {0})]
        rhs = int(0 in eq)
        if not terms:
            if rhs:
                clauses.append([])
            continue
        acc = terms[0]
        for b in terms[1:]:
            z = fresh()
            clauses.extend([[acc, b, -z], [-acc, -b, -z],
                            [acc, -b, z], [-acc, b, z]])
            acc = z
        clauses.append([acc if rhs else -acc])
    return nvars, clauses


def cnf_text(nvars, clauses):
    return f'p cnf {nvars} {len(clauses)}\n' + ''.join(
        ' '.join(map(str, c)) + ' 0\n' for c in clauses)


def exact_factor_clauses_multi(E, bases, subgroup_order, offsets=None):
    """Block every factor-coordinate assignment that cannot lift as required.

    This table depends only on the frozen factor space and subgroup, never on
    the decomposition target. Variables are the first m consecutive k-blocks.
    """
    if offsets is None:
        offsets = [E.base_field()(0)]*len(bases)
    clauses = []
    valid_by_block = []
    invalid_by_block = []
    start = 0
    for basis, offset in zip(bases, offsets):
        k = len(basis)
        invalid = []
        valid = []
        for assignment in range(1 << k):
            x = E.base_field()(offset)
            for j, coefficient in enumerate(basis):
                if assignment & (1 << j):
                    x += coefficient
            try:
                lifts = E.lift_x(x, all=True)
            except ValueError:
                lifts = []
            (valid if any(subgroup_order*P == E(0) for P in lifts) else invalid).append(assignment)
        for assignment in invalid:
            clauses.append([-(start+j+1) if assignment & (1 << j) else start+j+1
                            for j in range(k)])
        valid_by_block.append(valid)
        invalid_by_block.append(invalid)
        start += k
    return clauses, valid_by_block, invalid_by_block


def exact_factor_clauses(E, basis, m, subgroup_order, offset=0):
    clauses, valid, invalid = exact_factor_clauses_multi(
        E, [basis]*m, subgroup_order, [offset]*m)
    return clauses, valid[0], invalid[0]


def evaluate(equations, bits):
    return all(sum(all(bits[i-1] for i in mon_ids(mon)) for mon in eq) % 2 == 0
               for eq in equations)


def witness_multi(E, bases, R, bits, subgroup_order, offsets=None):
    """Recover signs, require subgroup membership, and verify the group sum."""
    if offsets is None:
        offsets = [E.base_field()(0)]*len(bases)
    xs = []
    start = 0
    for basis, offset in zip(bases, offsets):
        xs.append(sum((basis[j]*bits[start+j] for j in range(len(basis))),
                      E.base_field()(offset)))
        start += len(basis)
    choices = []
    for x in xs:
        try:
            lifts = E.lift_x(x, all=True)
        except ValueError:
            return None
        choices.append([P for P in lifts if subgroup_order * P == E(0)])
    for points in itertools.product(*choices):
        if sum(points, E(0)) == R:
            return [[int(P[0].to_integer()), int(P[1].to_integer())] for P in points]
    return None


def witness(E, basis, m, R, bits, subgroup_order, offset=0):
    return witness_multi(E, [basis]*m, R, bits, subgroup_order, [offset]*m)


def solve_verified(E, basis, m, R, r, solver, timeout, offset=0, trace_filter=False,
                   artifact_dir=None, exact_factor_filter=False,
                   factor_bases=None, factor_offsets=None):
    encoding_start = time.perf_counter()
    if trace_filter and r % 2 == 0:
        raise ValueError('The subgroup trace filter requires odd subgroup order')
    bases = factor_bases if factor_bases is not None else [basis]*m
    offsets = factor_offsets if factor_offsets is not None else [offset]*m
    if len(bases) != m or len(offsets) != m:
        raise ValueError('m must match factor_bases and factor_offsets')
    equations, names = descend_multi(E.base_field(), bases, R[0], E.a6(), offsets,
                                    E.a2().trace() if trace_filter else None)
    nvars, clauses = cnf(equations, len(names))
    factor_clauses = []
    valid_factor_assignments = []
    invalid_factor_assignments = []
    if exact_factor_filter:
        factor_clauses, valid_factor_assignments, invalid_factor_assignments = \
            exact_factor_clauses_multi(E, bases, r, offsets)
        clauses.extend(factor_clauses)
    encoding_seconds = time.perf_counter() - encoding_start
    input_sha256 = hashlib.sha256(cnf_text(nvars, clauses).encode()).hexdigest()
    original_clauses = len(clauses)
    if artifact_dir is not None:
        artifact_dir = Path(artifact_dir)
        artifact_dir.mkdir(parents=True, exist_ok=False)
        (artifact_dir / 'input.anf').write_text(anf_text(equations, len(names)))
        (artifact_dir / 'variables.json').write_text(json.dumps(names, indent=2) + '\n')
        (artifact_dir / 'factor-filter.json').write_text(json.dumps(dict(
            enabled=exact_factor_filter, subgroup_order=int(r),
            valid_assignments=valid_factor_assignments,
            invalid_assignments=invalid_factor_assignments,
            clauses=factor_clauses), indent=2) + '\n')
    rejects = 0
    start = time.perf_counter()
    solver_seconds = 0.0
    with (tempfile.TemporaryDirectory() if artifact_dir is None else nullcontext(artifact_dir)) as directory:
        while True:
            path = Path(directory) / f'attempt-{rejects:04d}.cnf'
            path.write_text(cnf_text(nvars, clauses))
            remaining = timeout - (time.perf_counter() - start)
            if remaining <= 0:
                status, cert = 'timeout', None
                break
            tick = time.perf_counter()
            try:
                result = subprocess.run([solver, str(path)], capture_output=True, text=True, timeout=remaining)
            except subprocess.TimeoutExpired:
                status, cert = 'timeout', None
                solver_seconds += time.perf_counter() - tick
                if artifact_dir is not None:
                    (artifact_dir / f'attempt-{rejects:04d}.timeout.json').write_text(
                        json.dumps(dict(command=[solver, str(path)], timeout_seconds=remaining)) + '\n')
                break
            solver_seconds += time.perf_counter() - tick
            if artifact_dir is not None:
                prefix = artifact_dir / f'attempt-{rejects:04d}'
                prefix.with_suffix('.stdout').write_text(result.stdout)
                prefix.with_suffix('.stderr').write_text(result.stderr)
                prefix.with_suffix('.json').write_text(json.dumps(dict(
                    command=result.args, returncode=result.returncode,
                    input_sha256=hashlib.sha256(path.read_bytes()).hexdigest())) + '\n')
            if result.returncode == 20 and 's UNSATISFIABLE' in result.stdout:
                status, cert = 'exhausted_polynomial_models', None
                break
            if result.returncode != 10 or 's SATISFIABLE' not in result.stdout:
                raise RuntimeError(f'Solver failed: {result.returncode}: {result.stderr}')
            values = {abs(int(v)): int(v) > 0 for line in result.stdout.splitlines()
                      if line.startswith('v ') for v in line.split()[1:] if int(v)}
            if any(i not in values for i in range(1, len(names)+1)):
                raise RuntimeError('Solver returned an incomplete model')
            bits = [int(values[i]) for i in range(1, len(names)+1)]
            if not evaluate(equations, bits):
                raise RuntimeError('Solver model fails original ANF')
            cert = witness_multi(E, bases, R, bits, r, offsets)
            if cert is not None:
                status = 'verified_relation'
                break
            # Lift failure is independent of the chain auxiliaries. Block the
            # factor-coordinate projection, not merely this full SAT assignment.
            factor_bits = sum(map(len, bases))
            clauses.append([-(i+1) if b else i+1 for i, b in enumerate(bits[:factor_bits])])
            rejects += 1
    return dict(status=status, witness=cert, rejected_projections=rejects,
                input_sha256=input_sha256, encoding_wall_seconds=encoding_seconds,
                solve_wall_seconds=solver_seconds, total_wall_seconds=time.perf_counter()-start,
                boolean_variables=len(names), cnf_variables=nvars, cnf_clauses=original_clauses,
                anf_equations=len(equations), exact_factor_filter=exact_factor_filter,
                factor_filter_valid_assignments=[len(row) for row in valid_factor_assignments],
                factor_filter_invalid_assignments=[len(row) for row in invalid_factor_assignments],
                factor_filter_clauses=len(factor_clauses))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=7)
    parser.add_argument('--dimension', type=int, default=3)
    parser.add_argument('--m', type=int, choices=[2, 3], default=3)
    parser.add_argument('--a', type=lambda x: int(x, 0), default=1)
    parser.add_argument('--b', type=lambda x: int(x, 0), default=1)
    parser.add_argument('--seed', type=int, default=20260909)
    parser.add_argument('--targets', type=int, default=4)
    parser.add_argument('--solver', default='cadical')
    parser.add_argument('--timeout', type=float, default=30)
    parser.add_argument('--no-trace-filter', action='store_true',
                        help='Disable necessary odd-subgroup trace constraints for ablation')
    parser.add_argument('--artifacts', type=Path,
                        help='New directory retaining source, ANF/CNF, variable maps and solver output')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 3 <= args.n <= 13 or args.n % 2 == 0:
        parser.error('This exhaustive diagnostic is limited to odd 3 <= n <= 13')
    if args.targets < 1 or args.timeout <= 0:
        parser.error('targets and timeout must be positive')
    if not 1 <= args.dimension <= min(args.n, 6):
        parser.error('Exhaustive diagnostic requires 1 <= dimension <= min(n, 6)')
    if not 0 <= args.a < 2**args.n or not 1 <= args.b < 2**args.n:
        parser.error('Require 0 <= a < 2^n and 0 < b < 2^n in polynomial-bit representation')
    if args.output.exists():
        parser.error('Output already exists; choose a new path')
    if args.artifacts is not None:
        if args.artifacts.exists():
            parser.error('Artifact directory already exists; choose a new path')
        args.artifacts.mkdir(parents=True)
        (args.artifacts / 'binary_factor_base.py').write_bytes(Path(__file__).read_bytes())
    set_random_seed(args.seed)
    F = GF(2**args.n, name='a')
    E = EllipticCurve(F, [1, F.from_integer(args.a), 0, 0, F.from_integer(args.b)])
    koblitz = E.a2()**2 == E.a2() and E.a6()**2 == E.a6()
    order = E.cardinality()
    r = max(p for p, _ in order.factor())
    if r == 2:
        parser.error('Diagnostic requires an odd prime subgroup')
    if (order // r) % r == 0:
        parser.error('Diagnostic requires a unique order-r subgroup')
    P = next(Q for Q in E if Q != E(0) and r * Q == E(0))
    rng = random.Random(args.seed)
    scalars = rng.sample(range(1, int(r)), min(args.targets, int(r)-1))
    bases, dims, factors = candidates(F, args.dimension, args.seed)
    report = dict(scope='toy implementation diagnostic; not harness-admitted research evidence',
                  environment=dict(sage=SAGE_VERSION, python=platform.python_version(),
                                   platform=platform.platform(),
                                   solver_version=subprocess.run([args.solver, '--version'],
                                       capture_output=True, text=True, check=True).stdout.strip()),
                  parameters=vars(args) | {'output': str(args.output),
                                          'artifacts': str(args.artifacts) if args.artifacts else None},
                  field_modulus=str(F.modulus()), curve_order=int(order), subgroup_order=int(r),
                  curve_type='Koblitz over F2' if koblitz else 'general ordinary binary',
                  subgroup_generator=[int(P[0].to_integer()), int(P[1].to_integer())],
                  target_scalars=scalars, invariant_dimensions=dims, frobenius_factors=factors, candidates=[])
    affine_bases = [(name, basis, F(0)) for name, basis in bases]
    affine_bases += [(name + '-affine-one', basis, F(1)) for name, basis in bases
                     if F(1) not in span(basis)]
    for name, basis, offset in affine_bases:
        values = span(basis, offset)
        fb = []
        for x in values:
            try:
                fb.extend(Q for Q in E.lift_x(x, all=True) if r*Q == E(0))
            except ValueError:
                pass
        # Exhaustive group-law control includes repeated summands and infinity
        # intermediate sums; no Semaev or SAT code is used in this oracle.
        sums = {E(0)}
        for _ in range(args.m):
            sums = {S + Q for S in sums for Q in fb}
        entry = dict(name=name, basis=[int(b.to_integer()) for b in basis],
                     affine_offset=int(offset.to_integer()),
                     subspace_size=len(values), subgroup_factor_points=len(fb),
                     frobenius_invariant=all(x*x in values for x in values),
                     exact_nonzero_target_coverage=len(sums - {E(0)}), runs=[])
        if entry['frobenius_invariant'] and koblitz:
            remaining_points = set(fb)
            orbit_sizes = []
            while remaining_points:
                Q = next(iter(remaining_points))
                orbit = set()
                for _ in range(args.n):
                    orbit.update([Q, -Q])
                    Q = E(Q[0]**2, Q[1]**2)
                if not orbit <= set(fb):
                    raise RuntimeError('Factor base is not closed under signed Frobenius')
                remaining_points -= orbit
                orbit_sizes.append(len(orbit))
            entry['signed_frobenius_orbit_sizes'] = sorted(orbit_sizes)
            entry['orbit_columns'] = len(orbit_sizes)
        for scalar in scalars:
            R = scalar * P
            result = solve_verified(E, basis, args.m, R, r, args.solver, args.timeout, offset,
                                    trace_filter=not args.no_trace_filter,
                                    artifact_dir=args.artifacts / name / str(scalar) if args.artifacts else None)
            result.update(target_scalar=scalar, target=[int(R[0].to_integer()), int(R[1].to_integer())],
                          group_oracle_decomposable=R in sums)
            if result['status'] == 'verified_relation' and R not in sums:
                raise RuntimeError('Relation disagrees with independent group oracle')
            entry['runs'].append(result)
        report['candidates'].append(entry)
    report['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
