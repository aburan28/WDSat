"""Archive completed extraction/frontier follow-ups with a hash index."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile

RUNS=['dimension-frontier-20260909-a','extraction-replication-20260909-a','frontier-extraction-20260909-a']
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path);args=parser.parse_args()
    archive=args.root/'followup-artifacts-20260909.tar.gz';index=args.root/'followup-evidence-index-20260909.json'
    if archive.exists() or index.exists():parser.error('Outputs already exist')
    files=[]
    with tarfile.open(archive,'w:gz',compresslevel=9) as tar:
        for name in RUNS:
            path=args.root/name
            if not path.exists():raise FileNotFoundError(path)
            tar.add(path,arcname=name);files.extend(p for p in path.rglob('*') if p.is_file())
    extraction=json.loads((args.root/'extraction-replication-20260909-a'/'result.json').read_text())
    frontier=json.loads((args.root/'frontier-extraction-20260909-a'/'result.json').read_text())
    def arms(report):
        return [{k:row[k] for k in ['name','column_count','rank','statuses','charged_seconds','extraction'] if k in row}
                for row in report['arms']]
    data=dict(scope='toy extraction evidence; independent replication pending',runs=RUNS,
              archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,file_count=len(files),
              file_sha256={str(p.relative_to(args.root)):sha(p) for p in files},
              extraction_arms=arms(extraction),frontier_arms=arms(frontier))
    index.write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(dict(bytes=archive.stat().st_size,files=len(files))))

if __name__=='__main__':main()
