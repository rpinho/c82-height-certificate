"""Independent root coverage audit using roots in X and the quadratic inverse of T.
No imports from the submitted verifier or inherited checkers.
"""
import json, sys, math, hashlib, time
from pathlib import Path
from fractions import Fraction
from flint import fmpz_poly, acb, arb, ctx
ctx.prec=256;ctx.threads=1
start=time.time();data=json.loads(Path(sys.argv[1]).read_bytes());records=[];D=0;mah=arb(0)
for role in ['denominator','numerator']:
 for rec in data[role]:
  a=rec['coeffs'];d=len(a)-1
  # Binomial expansion of sum a_k z^k(1-z)^k, rather than Horner composition.
  coeff=[0]*(2*d+1)
  for k,c in enumerate(a[::-1]):
   for j in range(k+1):coeff[k+j]+=c*math.comb(k,j)*(-1)**j
  composed=fmpz_poly(coeff); original=fmpz_poly(a[::-1])
  assert composed.degree()==2*d and int(composed.leading_coefficient())==a[0]*(-1)**d
  # Verify symmetry exactly by substitution z -> 1-z.
  reflected=[sum(coeff[k]*math.comb(k,j)*(-1)**j for k in range(j,len(coeff))) for j in range(len(coeff))]
  assert reflected==coeff
  rootsX=original.complex_roots();assert sum(m for r,m in rootsX)==d
  rootsZ=[]
  for x,m in rootsX:
   sq=(1-4*x).sqrt()
   rootsZ.extend([((1+sq)/2,m),((1-sq)/2,m)])
  assert sum(m for r,m in rootsZ)==composed.degree()
  # An independently composed polynomial must vanish on every mapped root ball.
  for z,m in rootsZ:
   value=acb(0)
   for c in coeff[::-1]:value=value*z+c
   assert value.contains(0)
  # Pairwise disjoint balls, combined with the certified X multiplicities, ensure
  # no quadratic preimage was counted twice. No X=1/4 occurs in this certificate.
  atquarter=sum(Fraction(c,4**k) for k,c in enumerate(a[::-1]));assert atquarter!=0
  for j,(z,m) in enumerate(rootsZ):
   for other,_ in rootsZ[:j]:assert not z.overlaps(other),(rec['name'],z,other)
  count=sum(m for _,m in rootsZ)
  entry=dict(role=role,name=rec['name'],degree=d,composed_degree=composed.degree(),enclosed_count=count,maximum_multiplicity=max(m for _,m in rootsZ),min_distance_to_circle=str(min((abs(abs(z)-1).lower() for z,m in rootsZ))),exact_symmetry=True,independent_binomial_composition=True)
  if role=='denominator':
   m=rec['multiplicity'];D+=d*m
   v=arb(abs(a[0])).log()+sum((multiplicity*abs(z).max(1).log() for z,multiplicity in rootsZ),arb(0))
   mah+=m*v;entry['mahler']=str(v)
  records.append(entry)
result=dict(sha256=hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest(),passed=True,precision=ctx.prec,D=D,weighted_denominator_root_count=2*D,mahler_normalized=str(mah/D),largest_original_degree=max(x['degree'] for x in records),largest_composed_degree=max(x['composed_degree'] for x in records),records=records,seconds=time.time()-start)
Path(sys.argv[2]).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2),flush=True)
