"""Validated mixed ANF/CNF adapter for the native WDSat binary."""
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

import binary_factor_base as lab


def formula_valid(text,bits):
    for row in text.splitlines()[1:]:
        tokens=row.split()
        if not tokens:continue
        if tokens[0]=='x':
            parity=0;pos=1
            while tokens[pos]!='0':
                token=tokens[pos];pos+=1
                if token=='T':parity^=1
                elif token.startswith('.'):
                    degree=int(token[1:]);parity^=int(all(bits[int(v)-1] for v in tokens[pos:pos+degree]));pos+=degree
                else:parity^=bits[int(token)-1]
            if parity!=1:return False
        elif not any((int(token)>0)==bool(bits[abs(int(token))-1]) for token in tokens if token!='0'):
            return False
    return True


def mixed_text(equations,nvars,clauses):
    anf=lab.anf_text(equations,nvars).splitlines();header=anf[0].split();header[3]=str(int(header[3])+len(clauses))
    return '\n'.join([' '.join(header),*anf[1:],*[' '.join(map(str,c))+' 0' for c in clauses]])+'\n'


def solve(E,bases,offsets,R,r,solver,timeout,artifact_dir):
    start=time.perf_counter();equations,names=lab.descend_multi(E.base_field(),bases,R[0],E.a6(),offsets,E.a2().trace())
    clauses,valid,invalid=lab.exact_factor_clauses_multi(E,bases,r,offsets);text=mixed_text(equations,len(names),clauses)
    encoding_seconds=time.perf_counter()-start;artifact_dir=Path(artifact_dir);artifact_dir.mkdir(parents=True,exist_ok=False);path=artifact_dir/'input.anf';path.write_text(text)
    (artifact_dir/'variables.json').write_text(json.dumps(names,indent=2)+'\n');(artifact_dir/'filter.json').write_text(json.dumps(dict(valid=valid,invalid=invalid,clauses=clauses),indent=2)+'\n')
    tick=time.perf_counter()
    try:result=subprocess.run([str(Path(solver).resolve()),'-i',str(path),'-x'],capture_output=True,text=True,timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        return dict(status='timeout',encoding_seconds=encoding_seconds,solver_seconds=time.perf_counter()-tick,
                    stdout=(exc.stdout or b'').decode() if isinstance(exc.stdout,bytes) else exc.stdout,
                    stderr=(exc.stderr or b'').decode() if isinstance(exc.stderr,bytes) else exc.stderr,
                    input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),witness=None)
    solver_seconds=time.perf_counter()-tick;(artifact_dir/'stdout').write_text(result.stdout);(artifact_dir/'stderr').write_text(result.stderr)
    if result.returncode!=0:raise RuntimeError(f'WDSat failed {result.returncode}: {result.stderr}')
    if 'UNSAT' in result.stdout:return dict(status='unsat',encoding_seconds=encoding_seconds,solver_seconds=solver_seconds,
        stdout=result.stdout,stderr=result.stderr,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),witness=None)
    models=[line for line in result.stdout.splitlines() if re.fullmatch(r'[01]{%d}'%len(names),line)]
    if not models:raise RuntimeError('WDSat returned no complete model')
    bits=list(map(int,models[-1]))
    if not formula_valid(text,bits):raise RuntimeError('WDSat model fails independent mixed formula evaluator')
    witness=lab.witness_multi(E,bases,R,bits,r,offsets)
    return dict(status='verified_relation' if witness else 'algebraic_nonrelation',encoding_seconds=encoding_seconds,
                solver_seconds=solver_seconds,stdout=result.stdout,stderr=result.stderr,
                input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),witness=witness)
