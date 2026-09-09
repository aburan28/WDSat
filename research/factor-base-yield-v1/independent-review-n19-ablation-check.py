"""Independent n19 ablation certificate audit; no producer-module imports or log oracle.
Run with DOT_SAGE in a writable location and sage -python.
"""
import hashlib,itertools,json,random,re,statistics
from pathlib import Path
from sage.all import GF,PolynomialRing,EllipticCurve,matrix,vector,ZZ
root=Path(__file__).resolve().parent
rawroot=root/'n19-same-dimension-ablation-20260909-a'
data=json.loads((rawroot/'result.json').read_text())
protocol=json.loads((rawroot/'n19-same-dimension-ablation.json').read_text())
pr=PolynomialRing(GF(2),'z');z=pr.gen()
F=GF(2**19,'a',modulus=z**19+z**5+z**2+z+1)
E=EllipticCurve(F,[1,1,0,0,1]);O=E(0);r=262543;L=GF(r)
assert ZZ(r).is_prime() and E.cardinality()==2*r
P=E(F.from_integer(393246),F.from_integer(31114));assert r*P==O and P!=O
assert data['P']==[393246,31114] and data['subgroup_order']==r
pol=PolynomialRing(L,'t');t=pol.gen()
roots=[int(x) for x,m in (t*t-t+2).roots()]
def phi(p,e): return O if p==O else E(p[0]**(2**(e%19)),p[1]**(2**(e%19)))
lam,=[x for x in roots if x*P==phi(P,1)]
assert lam==84450==data['frobenius_eigenvalue']
def key(p): return None if p==O else (int(p[0].to_integer()),int(p[1].to_integer()))
def point(pair):return E(F.from_integer(pair[0]),F.from_integer(pair[1]))
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


def parse_equations(text):
 out=[]
 for line in text.splitlines()[1:]:
  tok=line.split()
  if tok[0]!='x':continue
  equation={0};pos=1
  while tok[pos]!='0':
   val=tok[pos];pos+=1
   if val=='T':mon=0
   elif val.startswith('.'):
    deg=int(val[1:]);mon=sum(1<<(int(v)-1) for v in tok[pos:pos+deg]);pos+=deg
   else:mon=1<<(int(val)-1)
   if mon in equation:equation.remove(mon)
   else:equation.add(mon)
  out.append(equation)
 return out

# Independently expand over an ordinary multivariate polynomial ring,
# then reduce Boolean exponents. This differs from the producer's bit-mask algebra.
poly=PolynomialRing(F,37,'b');variables=poly.gens()
def reduced_coordinate_equations(expr):
 coeff={}
 for powers,c in expr.dict().items():
  mon=sum(1<<i for i,e in enumerate(powers) if e)
  coeff[mon]=coeff.get(mon,F(0))+c
 rows=[set() for i in range(19)]
 for mon,c in coeff.items():
  v=int(c.to_integer())
  for bit in range(19):
   if v>>bit&1:rows[bit].add(mon)
 return [row for row in rows if row]

def s3_expr(x,y,z):return (x*y+x*z+y*z)**2+x*y*z+1

configs={};result_checks=[]
for arm in data['streams'][0]['arms']:
 name=arm['name'];basis=arm['basis'];k=len(basis)
 assert k==6 and arm['offset']==1 and arm['exponents']==[0,1,2]
 basic=[F.from_integer(x) for x in basis]
 assigns=[[F(1)+sum((b for j,b in enumerate(basic) if mask>>j&1),F(0)) for mask in range(64)]]
 assigns=[[x**(2**e) for x in assigns[0]] for e in [0,1,2]]
 pts=[[p for x in block for p in lift(x)] for block in assigns]
 cols=sorted({min(key(p),key(-p)) for p in pts[0]})
 valid=[[i for i,x in enumerate(block) if lift(x)] for block in assigns]
 invalid=[[i for i in range(64) if i not in v] for v in valid]
 clauses=[[(-1 if mask>>j&1 else 1)*(block*6+j+1) for j in range(6)] for block,inv in enumerate(invalid) for mask in inv]
 # Different grouping from producer: blocks zero and two first, then block one.
 pair={p0+p2 for p0 in pts[0] for p2 in pts[2]}
 covered={key(s+p1) for s in pair for p1 in pts[1]};covered.discard(None)
 assert len(covered)==arm['reachable_targets']
 cover_digest=hashlib.sha256(json.dumps(sorted(covered),separators=(',',':')).encode()).hexdigest()
 xs=[poly(1)+sum((variables[6*i+j]*b**(2**i) for j,b in enumerate(basic)),poly(0)) for i in range(3)]
 u=sum((variables[18+j]*F.gen()**j for j in range(19)),poly(0))
 first_equations=reduced_coordinate_equations(s3_expr(xs[0],xs[1],u))
 configs[name]=dict(assigns=assigns,pts=pts,cols=cols,valid=valid,invalid=invalid,clauses=clauses,covered=covered,
  coverage_digest=cover_digest,xs=xs,u=u,first_equations=first_equations,basis=basis)
 print('factor',name,'columns',len(cols),'coverage',len(covered),'digest',cover_digest,flush=True)

for stream,spec in zip(data['streams'],protocol['streams']):
 seed=stream['target_seed'];secret=stream['secret'];Q=point(stream['Q'])
 assert seed==spec['target_seed'] and secret==spec['secret'] and Q==secret*P
 rng=random.Random(seed+19);targets=[]
 while len(targets)<protocol['maximum_targets']:
  a,b=rng.randrange(r),rng.randrange(r);R=a*P+b*Q
  if R!=O:targets.append((a,b,R))
 for arm in stream['arms']:
  name=arm['name'];cfg=configs[name];cols=cfg['cols'];rows=[];aa=[];bb=[];extracted=None;sat=unsat=0
  assert arm['basis']==cfg['basis'] and arm['offset']==1 and arm['exponents']==[0,1,2]
  assert len(cols)==arm['column_count'] and arm['reachable_targets']==len(cfg['covered'])
  assert [rec['target_index'] for rec in arm['runs']]==list(range(len(arm['runs'])))
  for rec in arm['runs']:
   idx=rec['target_index'];a,b,R=targets[idx];assert (rec['a'],rec['b'])==(a,b)
   folder=rawroot/'raw'/str(seed)/name/str(idx);txt=(folder/'input.anf').read_text();out=(folder/'stdout').read_text()
   assert hashlib.sha256((folder/'input.anf').read_bytes()).hexdigest()==rec['input_sha256']
   assert out==rec['stdout'] and (folder/'stderr').read_text()==rec['stderr']
   assert json.loads((folder/'filter.json').read_text())==dict(valid=cfg['valid'],invalid=cfg['invalid'],clauses=cfg['clauses'])
   assert json.loads((folder/'variables.json').read_text())==[f'x{i}_{j}' for i in range(1,4) for j in range(6)]+[f'u1_{j}' for j in range(19)]
   assert list(map(int,txt.splitlines()[0].split()[2:]))==[37,len(txt.splitlines())-1]
   assert [list(map(int,line.split()[:-1])) for line in txt.splitlines()[1:] if not line.startswith('x')]==cfg['clauses']
   expected=cfg['first_equations']+reduced_coordinate_equations(s3_expr(cfg['u'],cfg['xs'][2],R[0]))
   assert parse_equations(txt)==expected
   reachable=key(R) in cfg['covered'];assert rec['reachable']==reachable
   assert rec['pipeline_seconds']>=rec['encoding_seconds']+rec['solver_seconds']
   if rec['status']=='unsat':
    assert not reachable and 'UNSAT' in out and rec['witness'] is None
    assert rec['rank_after']==(matrix(L,rows).rank() if rows else 0)
    unsat+=1;continue
   assert rec['status']=='verified_relation' and reachable;sat+=1
   bits=list(map(int,next(v for v in reversed(out.splitlines()) if re.fullmatch('[01]{37}',v))))
   eval_formula(txt,bits)
   xs=[cfg['assigns'][i][sum(bits[6*i+j]<<j for j in range(6))] for i in range(3)]
   u=F.from_integer(sum(bits[18+j]<<j for j in range(19)))
   assert s3_expr(xs[0],xs[1],u)==0 and s3_expr(u,xs[2],R[0])==0
   points=[point(p) for p in rec['witness']]
   assert [p[0] for p in points]==xs and all(r*p==O for p in points) and sum(points,O)==R
   coeff=[0]*len(cols)
   for e,p in enumerate(points):
    pre=phi(p,-e);col=min(key(pre),key(-pre));assert phi(pre,e)==p
    coeff[cols.index(col)]+=(1 if col==key(pre) else -1)*pow(lam,e,r)
   assert coeff==rec['coefficients']
   assert sum((c*point(col) for c,col in zip(coeff,cols)),O)==R
   rows.append(coeff);aa.append(a);bb.append(b);M=matrix(L,rows);rank=int(M.rank());assert rank==rec['rank_after']
   for dep in M.left_kernel().basis():
    ca=dep*vector(L,aa);cb=dep*vector(L,bb)
    if cb:
     assert extracted is None
     d=int(-ca/cb);assert d*P==Q and d==rec['verified_extracted_log']==secret
     extracted=(idx,d)
    else:assert ca==0
  assert extracted==(arm['extraction']['target_index'],arm['extraction']['log'])
  assert extracted[0]==len(arm['runs'])-1 and sat==arm['extraction']['relations'] and rank==arm['rank']
  solver=sum(rec['solver_seconds'] for rec in arm['runs']);pipeline=sum(rec['pipeline_seconds'] for rec in arm['runs'])
  assert solver==arm['solver_seconds']==arm['extraction']['solver_seconds']
  assert pipeline==arm['pipeline_seconds']==arm['extraction']['pipeline_seconds']
  record=dict(seed=seed,secret=secret,arm=name,targets=len(arm['runs']),sat=sat,unsat=unsat,rank=rank,
   first_dependency=list(extracted),solver_seconds=solver,pipeline_seconds=pipeline)
  result_checks.append(record);print('verified',record,flush=True)

# Reproduce the sparse candidate using rational/subgroup membership predicates only.
# Compute actual half-trace lifts and test subgroup membership, without a log table.
directions=[F.gen()**i+F((F.gen()**i).trace()) for i in range(1,19)]
membership={};best=None;candidates=0
for subset in itertools.combinations(range(18),6):
 vectors=[int(directions[i].to_integer()) for i in subset];values=[1]
 for b in vectors:values += [x^b for x in values]
 for value in values:
  if value not in membership:
   x=F.from_integer(value);assert x and x.trace()==1
   membership[value]=bool(lift(x))
   assert membership[value]==(not (x+1+1/(x*x)).trace())
 score=sum(membership[x] for x in values);exponents=tuple(i+1 for i in subset)
 candidate=(score,-sum(exponents),tuple(-i for i in exponents))
 if best is None or candidate>best[0]:best=(candidate,exponents,vectors)
 candidates+=1
assert list(best[1])==[3,4,7,10,12,15] and best[0][0]==45
assert best[2]==configs['density-selected-sparse-k6']['basis']
selection=dict(candidates=candidates,membership_inputs=len(membership),winner_exponents=list(best[1]),winner_basis=best[2],valid_assignments=int(best[0][0]),uses_logarithms=False)
print('membership-only selection',selection,flush=True)
metrics={name:{metric:statistics.median(row[metric] for row in result_checks if row['arm']==name) for metric in ['targets','solver_seconds','pipeline_seconds']} for name in configs}
ratios={metric:metrics['fixed-trace-k6'][metric]/metrics['density-selected-sparse-k6'][metric] for metric in metrics['fixed-trace-k6']}
report=dict(snapshot='3c5e675b0ccd3685076fa5e58b4262ed2007d24c',eigenvalue_roots=roots,eigenvalue=lam,checks=result_checks,
 coverage=[dict(arm=name,count=len(cfg['covered']),fraction=len(cfg['covered'])/(r-1),digest=cfg['coverage_digest'],columns=len(cfg['cols']),valid_assignments=[len(row) for row in cfg['valid']],membership_clauses=len(cfg['clauses'])) for name,cfg in configs.items()],selection=selection,medians=metrics,ratios_of_medians=ratios,
 counts=dict(formulas=sum(x['targets'] for x in result_checks),sat_models=sum(x['sat'] for x in result_checks),semantic_unsat=sum(x['unsat'] for x in result_checks)))
(root/'independent-review-n19-ablation-check-result.json').write_text(json.dumps(report,indent=2)+'\n')
print('ALL PASS',report['counts'],metrics,ratios,flush=True)
