"""Create a durable compressed archive and concise evidence index."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile


RUNS=['run-20260909-a','sparse-20260909-a','sat-rank-20260909-a',
      'lift-filter-20260909-a','frobenius-schedule-20260909-a',
      'frobenius-schedule-20260909-b','frobenius-sat-20260909-a']


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path)
    args=parser.parse_args();archive=args.root/'run-artifacts-20260909.tar.gz';index=args.root/'evidence-index-20260909.json'
    if archive.exists() or index.exists():parser.error('Archive or index exists')
    files=[]
    with tarfile.open(archive,'w:gz',compresslevel=9) as tar:
        for name in RUNS:
            path=args.root/name
            if not path.exists():raise FileNotFoundError(path)
            tar.add(path,arcname=name)
            files.extend(p for p in path.rglob('*') if p.is_file())
        extra=args.root/'cofactor-control-20260909.json';tar.add(extra,arcname=extra.name);files.append(extra)
    def load(run,n):return json.loads((args.root/run/f'n{n}-result.json').read_text())
    lift={str(n):[{k:a[k] for k in ['name','rank','statuses','measured_seconds','first_extraction']}
                  for a in load('lift-filter-20260909-a',n)['arms']] for n in [11,17]}
    frob={str(n):[{k:a[k] for k in ['name','rank','statuses','measured_seconds','first_extraction']}
                  for a in load('frobenius-sat-20260909-a',n)['arms']] for n in [11,17]}
    sparse={str(n):{k:load('sparse-20260909-a',n)[k] for k in
                    ['actual_subgroup_points','sign_folded_columns','exact_nonzero_target_coverage','encoding']}
            for n in [7,11,17]}
    data=dict(scope='implementation and toy research evidence; no breakthrough claim',runs=RUNS,
              archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,
              archived_file_count=len(files),file_sha256={str(p.relative_to(args.root)):sha(p) for p in files},
              sparse_selection=sparse,exact_lift_sat=lift,frobenius_sat=frob,
              caveats=['n=17 batches have 16 targets','watchdog timeouts are operational observations',
                       'Frobenius image method is prior art','no independent adversarial review yet'])
    index.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(archive=str(archive),bytes=archive.stat().st_size,files=len(files))))


if __name__=='__main__':main()
