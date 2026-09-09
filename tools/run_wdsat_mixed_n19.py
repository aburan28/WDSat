"""Run and independently validate WDSat on retained n19 mixed ANF/CNF inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

from sage.all import GF, EllipticCurve
import binary_factor_base as lab
from factor_base_sat_rank import RelationRank, columns, point_key
from factor_base_frobenius_sat import rewritten_coefficients, frobenius_power
from factor_base_yield_screen import subgroup_generator


def mixed_text(anf, filter_data):
    lines=anf.splitlines();header=lines[0].split();clauses=filter_data['clauses']
    header[3]=str(int(header[3])+len(clauses));rows=[' '.join(header),*lines[1:]]
    rows.extend(' '.join(map(str,clause))+' 0' for clause in clauses)
    return '\n'.join(rows)+'\n'


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
        else:
            if not any((int(token)>0)==bool(bits[abs(int(token))-1]) for token in tokens if token!='0'):return False
    return True


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,required=True);parser.add_argument('--solver',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();config_path=args.source/'n19-wdsat-amendment.json';config=json.loads(config_path.read_text())
    admission=json.loads((args.source/'n19-sat-admission-20260909-a'/'result.json').read_text());args.output.mkdir(parents=True,exist_ok=False)
    for path in [Path(__file__),args.solver,config_path]:(args.output/path.name).write_bytes(path.read_bytes())
    F=GF(2**19,name='a');E=EllipticCurve(F,[1,1,0,0,1]);r=int(admission['subgroup_order']);P=subgroup_generator(E,r);Q=E(F.from_integer(admission['Q'][0]),F.from_integer(admission['Q'][1]))
    piP=frobenius_power(E,P,1);walk=E(0);lam=None
    for j in range(r):
        if walk==piP:lam=j;break
        walk+=P
    if lam is None:raise RuntimeError('No Frobenius eigenvalue')
    wanted=['trace-k7-best','sparse-k6-best'];states={}
    for name in wanted:
        arm=next(a for a in admission['arms'] if a['name']==name);arm=dict(arm)
        arm['cofactor_multiplier']=1;states[name]=dict(arm=arm,columns=columns(E,r,arm),ranks={mode:RelationRank(r) for mode in ['plain','gaussian']},runs=[])
    raw=args.output/'raw';raw.mkdir()
    for name,state in states.items():
        arm=state['arm'];base=[F.from_integer(v) for v in arm['basis']];offset=F.from_integer(arm['offset']);bases=[[v**(2**e) for v in base] for e in arm['exponents']];offsets=[offset**(2**e) for e in arm['exponents']]
        source_arm=args.source/'n19-sat-admission-20260909-a'/'raw'/name
        for target_index in range(12):
            original=source_arm/str(target_index);anf=(original/'input.anf').read_text();filter_data=json.loads((original/'factor-filter.json').read_text());mixed=mixed_text(anf,filter_data)
            for mode,flags in [('plain',[]),('gaussian',['-x'])]:
                run_dir=raw/name/str(target_index)/mode;run_dir.mkdir(parents=True);input_path=run_dir/'input.anf';input_path.write_text(mixed)
                start=time.perf_counter()
                try:result=subprocess.run([str(args.solver),'-i',str(input_path),*flags],capture_output=True,text=True,timeout=config['watchdog_seconds'])
                except subprocess.TimeoutExpired as exc:
                    record=dict(arm=name,target_index=target_index,mode=mode,status='timeout',seconds=time.perf_counter()-start,stdout=(exc.stdout or b'').decode() if isinstance(exc.stdout,bytes) else exc.stdout,stderr=(exc.stderr or b'').decode() if isinstance(exc.stderr,bytes) else exc.stderr)
                    state['runs'].append(record);print(json.dumps({k:record[k] for k in ['arm','target_index','mode','status','seconds']}),flush=True);continue
                seconds=time.perf_counter()-start;(run_dir/'stdout').write_text(result.stdout);(run_dir/'stderr').write_text(result.stderr)
                if result.returncode!=0:raise RuntimeError(f'WDSat process failed {result.returncode}: {result.stderr}')
                if 'UNSAT' in result.stdout:status='unsat';bits=None;witness=None
                else:
                    models=[line for line in result.stdout.splitlines() if re.fullmatch(r'[01]{%d}'%len(json.loads((original/'variables.json').read_text())),line)]
                    if not models:raise RuntimeError('WDSat returned no complete model')
                    bits=list(map(int,models[-1]));
                    if not formula_valid(mixed,bits):raise RuntimeError('WDSat model fails independent mixed-formula evaluation')
                    target=next(t for t in arm['runs'] if t['target_index']==target_index);R=target['a']*P+target['b']*Q
                    witness=lab.witness_multi(E,bases,R,bits,r,offsets);status='verified_relation' if witness else 'algebraic_nonrelation'
                    if witness:
                        coefficients=rewritten_coefficients(E,witness,arm['exponents'],lam,r,state['columns']);lhs=E(0)
                        for key,pos in state['columns'].items():lhs+=coefficients[pos]*E(F.from_integer(key[0]),F.from_integer(key[1]))
                        if lhs!=R:raise RuntimeError('Rewritten relation failed')
                        d=state['ranks'][mode].add(coefficients,target['a'],target['b'])
                        if d is not None and d*P!=Q:raise RuntimeError('Dependency log failed')
                record=dict(arm=name,target_index=target_index,mode=mode,status=status,seconds=seconds,returncode=result.returncode,witness=witness,stdout=result.stdout,stderr=result.stderr,rank=len(state['ranks'][mode].rows))
                state['runs'].append(record);print(json.dumps({k:record[k] for k in ['arm','target_index','mode','status','seconds','rank']}),flush=True)
    report=dict(scope='n19 WDSat mixed-formula admission; not a performance breakthrough',solver_sha256=hashlib.sha256(args.solver.read_bytes()).hexdigest(),arms=[])
    for name,state in states.items():
        summaries={}
        for mode in ['plain','gaussian']:
            rows=[r for r in state['runs'] if r['mode']==mode];summaries[mode]=dict(statuses={s:sum(r['status']==s for r in rows) for s in sorted(set(r['status'] for r in rows))},seconds=sum(r['seconds'] for r in rows),rank=len(state['ranks'][mode].rows))
        report['arms'].append(dict(name=name,summaries=summaries,runs=state['runs']))
    (args.output/'result.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
