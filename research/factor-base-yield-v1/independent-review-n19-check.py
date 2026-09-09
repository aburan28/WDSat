import json, hashlib, random, re, statistics
from pathlib import Path
from sage.all import GF, PolynomialRing, EllipticCurve, matrix, vector, ZZ
root=Path(__file__).resolve().parent
ring=PolynomialRing(GF(2),'z');z=ring.gen();F=GF(2**19,'a',modulus=z**19+z**5+z**2+z+1)
E=EllipticCurve(F,[1,1,0,0,1]);r=262543;L=GF(r);O=E(0)
assert E.cardinality()==2*r and ZZ(r).is_prime()
def key(p): return None if p==O else (int(p[0].to_integer()),int(p[1].to_integer()))
def point(pair): return E(F.from_integer(pair[0]),F.from_integer(pair[1]))
def phi(p,e): return O if p==O else E(p[0]**(2**(e%19)),p[1]**(2**(e%19)))
P=point([393246,31114]);Q=point([267253,240749]);lam=84450
assert P!=O and r*P==O and Q==31415*P and phi(P,1)==lam*P
assert (lam*lam-lam+2)%r==0 and pow(lam,19,r)==1

def lift(x):
 if not x:return []
 c=x+1+1/(x*x)
 if c.trace():return []
 h=sum((c**(4**i) for i in range(10)),F(0));assert h*h+h==c
 pts=[E(x,x*h),E(x,x*(h+1))]
 return [p for p in pts if r*p==O]

def eval_formula(text,bits):
 lines=text.splitlines();assert int(lines[0].split()[3])==len(lines)-1
 for line in lines[1:]:
  tok=line.split();assert tok[-1]=='0'
  if tok[0]!='x':assert any((int(t)>0)==bool(bits[abs(int(t))-1]) for t in tok[:-1]);continue
  mon=[];at=1
  while at<len(tok)-1:
   t=tok[at];at+=1
   if t=='T':mon.append(1)
   elif t.startswith('.'):
    deg=int(t[1:]);mon.append(all(bits[int(v)-1] for v in tok[at:at+deg]));at+=deg
   else:mon.append(bits[int(t)-1])
  assert sum(mon)%2==1

configs={}; alltot=[]
for folder,seed in [('n19-wdsat-extraction-20260909-a',202609104),('n19-wdsat-replication-20260909-b',202609105),('n19-wdsat-replication-20260909-c',202609106)]:
 path=root/folder;d=json.loads((path/'result.json').read_text());targets=json.loads((path/'targets.json').read_text())
 assert d['P']==list(key(P)) and d['Q']==list(key(Q)) and d['frobenius_eigenvalue']==lam
 rng=random.Random(seed+19);gen=[]
 while len(gen)<96:
  a,b=rng.randrange(r),rng.randrange(r);R=a*P+b*Q
  if R!=O:gen.append((a,b,R))
 assert len(targets)==96
 for i,(a,b,R) in enumerate(gen):assert targets[i]==dict(index=i,a=a,b=b,target=list(key(R)))
 for arm in d['arms']:
  name=arm['name'];k=len(arm['basis']);exps=arm['exponents']
  if name not in configs:
   base=[F.from_integer(v) for v in arm['basis']];off=F.from_integer(arm['offset']);assign=[[off+sum((base[j] for j in range(k) if mask>>j&1),F(0)) for mask in range(2**k)]]
   assign=[ [x**(2**e) for x in assign[0]] for e in exps]
   pts=[[p for x in block for p in lift(x)] for block in assign]
   for block in assign:assert all(x.trace()==1 for x in block)
   columns=sorted({min(key(p),key(-p)) for p in pts[0]});assert len(columns)==arm['column_count']
   valid=[[i for i,x in enumerate(block) if lift(x)] for block in assign]
   invalid=[[i for i in range(2**k) if i not in v] for v in valid]
   clauses=[[(-1 if mask>>j&1 else 1)*(t*k+j+1) for j in range(k)] for t,row in enumerate(invalid) for mask in row]
   pairs={key(p1+p2) for p1 in pts[0] for p2 in pts[1]}
   configs[name]=(assign,pts,columns,valid,invalid,clauses,pairs)
   print('factor',name,'valid',[len(x) for x in valid],'points',[len(x) for x in pts],'columns',len(columns),'pair_count',len(pairs),flush=True)
  assign,pts,columns,valid,invalid,clauses,pairs=configs[name]
  rows=[];aa=[];bb=[];extracted=None;rank=0;nsat=nunsat=0
  assert [x['target_index'] for x in arm['runs']]==list(range(len(arm['runs'])))
  for rec in arm['runs']:
   idx=rec['target_index'];a,b,R=gen[idx];reused=folder.endswith('extraction-20260909-a') and idx<12
   if not reused or rec['status']=='verified_relation':assert rec['a']==a and rec['b']==b
   raw=root/'n19-wdsat-20260909-c'/'raw'/name/str(idx)/'gaussian' if reused else path/'raw'/name/str(idx)
   txt=(raw/'input.anf').read_text();out=(raw/'stdout').read_text();assert out==rec['stdout']
   assert (raw/'stderr').read_text()==rec['stderr']
   if not reused:
    assert hashlib.sha256((raw/'input.anf').read_bytes()).hexdigest()==rec['input_sha256']
    assert json.loads((raw/'filter.json').read_text())==dict(valid=valid,invalid=invalid,clauses=clauses)
   cnf=[list(map(int,x.split()[:-1])) for x in txt.splitlines()[1:] if not x.startswith('x')]
   assert cnf==clauses
   reachable=any(key(R-p3) in pairs for p3 in pts[2])
   if 'reachable' in rec:assert reachable==rec['reachable']
   if rec['status']=='unsat':assert not reachable and 'UNSAT' in out;nunsat+=1;continue
   assert rec['status']=='verified_relation' and reachable
   nsat+=1
   bits=list(map(int,next(x for x in reversed(out.splitlines()) if re.fullmatch('[01]{%d}'%(3*k+19),x))))
   eval_formula(txt,bits)
   xs=[assign[t][sum(bits[t*k+j]<<j for j in range(k))] for t in range(3)]
   u=F.from_integer(sum(bits[3*k+j]<<j for j in range(19)))
   def s3(x,y,t):return (x*y+x*t+y*t)**2+x*y*t+1
   assert s3(xs[0],xs[1],u)==0 and s3(u,xs[2],R[0])==0
   w=[point(v) for v in rec['witness']]
   assert [p[0] for p in w]==xs and all(r*p==O for p in w) and sum(w,O)==R
   coeff=[0]*len(columns)
   for p,e in zip(w,exps):
    pre=phi(p,-e);canon=min(key(pre),key(-pre));assert phi(pre,e)==p
    coeff[columns.index(canon)]+=(1 if key(pre)==canon else -1)*pow(lam,e,r)
   assert coeff==rec['coefficients']
   assert sum((c*point(col) for c,col in zip(coeff,columns)),O)==R
   rows.append(coeff);aa.append(a);bb.append(b);M=matrix(L,rows);rank=int(M.rank())
   assert rec.get('rank_after',rec.get('rank'))==rank
   for dep in M.left_kernel().basis():
    aval=dep*vector(L,aa);bval=dep*vector(L,bb)
    if bval:
     recovered=int(-aval/bval);assert recovered==31415 and recovered*P==Q
     assert extracted is None;extracted=(idx,recovered)
    else:assert aval==0
  assert extracted==(arm['extraction']['target_index'],arm['extraction']['log'])
  assert extracted[0]==len(arm['runs'])-1 and nsat==arm['extraction']['relations'] and rank==arm['rank']
  print('verified',folder,name,'SAT',nsat,'UNSAT',nunsat,'rank',rank,'first_dependency',extracted,flush=True)
  alltot.append(dict(folder=folder,arm=name,sat=nsat,unsat=nunsat,rank=rank,extraction=list(extracted)))
print('ALL PASS',sum(x['sat'] for x in alltot),'models',sum(x['unsat'] for x in alltot),'semantic UNSAT checks',flush=True)
# Independent full reachability counts using a subgroup walk and Python sets,
# separate from the producer's cyclic bitset implementation.
wanted={key(p) for cfg in configs.values() for block in cfg[1] for p in block}
logs={};walk=O
for scalar in range(1,r):
 walk+=P;tag=key(walk)
 if tag in wanted:logs[tag]=scalar
assert set(logs)==wanted
frontier=json.loads((root/'dimension-frontier-n19-20260909-a/result.json').read_text())
coverage=[]
for name,cfg in configs.items():
 point_sets=cfg[1];sets=[{logs[key(p)] for p in block} for block in point_sets]
 sums={(a+b)%r for a in sets[0] for b in sets[1]}
 sums={(a+b)%r for a in sums for b in sets[2]}
 count=len(sums-{0})
 candidate=next(c for c in frontier['candidates'] if c['name']+'-best'==name)
 expected=next(c['coverage'] for c in candidate['controls'] if c['name']=='best')
 assert count==expected
 record=dict(arm=name,coverage=count,nonzero_fraction=count/(r-1))
 coverage.append(record);print('coverage',record,flush=True)
(root/'independent-review-n19-check-result.json').write_text(json.dumps(dict(checks=alltot,coverage=coverage),indent=2)+'\n')
