"""Fresh C82 integrator: 4-node Gaussian quadrature and 8th derivative bounds.

Written for this review; imports no claimant/search/inherited verifier code.
Roots are isolated in X and mapped through (1 +/- sqrt(1-4X))/2.
The complete circle [0,1] is integrated, without symmetry doubling.

D_t^k log|exp(i*tau*t)-r| has (for k>=2) an Eulerian numerator:
 its modulus <= tau^k |r| E_{k-1}(|r|) / distance^k.
For k=8, E_7 = 1+120x+1191x^2+2416x^3+1191x^4+120x^5+x^6.
4-point Gaussian quadrature has error <= width^9 (4!)^4 / (9 (8!)^3) sup |g^(8)|.
Crossing cells use the whole-cell range enclosure, without Taylor quadrature.
"""
import argparse, hashlib, json, math, time, multiprocessing
from fractions import Fraction as Q
from pathlib import Path
from flint import arb,acb,fmpz_poly,ctx

def qball(x):
 x=Q(x);return arb(x.numerator)/arb(x.denominator)

def initialize(path,precision,tol):
 global GROUPS,TAU,TOL,NODES,WEIGHTS,CONSTANT,MAHLER,GAUSS_ERROR
 ctx.prec=precision;ctx.threads=1
 cert=json.loads(Path(path).read_bytes())
 D=sum((len(r['coeffs'])-1)*r['multiplicity'] for r in cert['denominator'])
 TAU=2*arb.pi();TOL=qball(tol);GROUPS=[];CONSTANT=arb(0);MAHLER=arb(0)
 for kind in ['denominator','numerator']:
  for rec in cert[kind]:
   a=rec['coeffs'];exponent=Q(-rec['multiplicity'],D) if kind=='denominator' else Q(rec['q_rational'])/D
   if not exponent:continue
   w=qball(exponent);CONSTANT+=w*arb(abs(a[0])).log()
   xr=fmpz_poly(a[::-1]).complex_roots();assert sum(m for r,m in xr)==len(a)-1
   roots=[]
   for x,m in xr:
    s=(1-4*x).sqrt()
    for r in [(1+s)/2,(1-s)/2]:
     ar=abs(r)
     # Eulerian polynomial evaluated by Horner with positive terms.
     E=arb(1)
     for c in [120,1191,2416,1191,120,1]:E=E*ar+c
     roots.append((r,w*m,abs(w)*m*TAU**8*ar*E,m))
     if kind=='denominator':MAHLER+=(-w)*m*ar.max(1).log()
   assert sum(m for _,m in xr)*2==2*(len(a)-1)
   GROUPS.append((w,roots))
   if kind=='denominator':MAHLER-=w*arb(abs(a[0])).log()
 s30=arb(30).sqrt();inner=((15-2*s30)/35).sqrt();outer=((15+2*s30)/35).sqrt()
 NODES=[-outer,-inner,inner,outer]
 WEIGHTS=[(18-s30)/36,(18+s30)/36,(18+s30)/36,(18-s30)/36]
 GAUSS_ERROR=qball(Q(math.factorial(4)**4,9*math.factorial(8)**3))

def value_at(t):
 u=acb(0,TAU*t).exp();v=CONSTANT
 for w,rs in GROUPS:
  prod=arb(1)
  for r,wm,c8,multiplicity in rs:
   # All actual factors are squarefree; retain multiplicity support explicitly.
   prod*=abs(u-r) if multiplicity==1 else abs(u-r)**multiplicity
  v+=w*prod.log()
 return v

def integrate_tile(index):
 # 64 disjoint exact dyadic cells cover the full circle.
 started=time.time();stack=[(index,6)];answer=arb(0);counts=[0,0,0];splits=0;maximum=6
 while stack:
  k,depth=stack.pop();maximum=max(maximum,depth)
  width=arb(1)/(1<<depth);mid=arb(2*k+1)/(1<<(depth+1));half=width/2
  center=acb(0,TAU*mid).exp();radius=TAU*half
  hi=CONSTANT;lo=CONSTANT;valid_hi=True;valid_lo=True;separated=True;M8=arb(0)
  for weight,roots in GROUPS:
   for r,wm,c8,multiplicity in roots:
    dist=abs(center-r);near=dist-radius;far=dist+radius
    if wm>0:
     hi+=wm*far.log()
     if near>0:lo+=wm*near.log()
     else:valid_lo=False
    else:
     lo+=wm*far.log()
     if near>0:hi+=wm*near.log()
     else:valid_hi=False
    if near>0:M8+=c8/near**8
    else:separated=False
  if valid_hi and hi<0:
   counts[0]+=1;continue
  positive=valid_lo and valid_hi and lo>0
  if not positive and separated:
   # Direct interval enclosure of g' on the arc, followed by the mean value
   # theorem. This retains cancellation in signed derivatives and uses no
   # midpoint Taylor sign test or fourth derivative from verify_exact.py.
   arc=acb(0,TAU*arb(mid,half)).exp();derivative=acb(0)
   for weight,roots in GROUPS:
    for r,wm,c8,multiplicity in roots:derivative+=wm*acb(0,TAU)*arc/(arc-r)
   center_value=value_at(mid)
   reach=abs(derivative.real)*half
   if center_value+reach<0:
    counts[0]+=1;continue
   positive=bool(center_value-reach>0)
  if positive:
   remainder=GAUSS_ERROR*width**9*M8
   if remainder<TOL*width/4:
    v=sum((weight*value_at(mid+half*node) for node,weight in zip(NODES,WEIGHTS)),arb(0))*half
    assert v.is_finite()
    answer+=(v+remainder).upper();counts[1]+=1;continue
  elif valid_hi and hi.max(0)<TOL/4:
   answer+=(width*hi.max(0)).upper();counts[2]+=1;continue
  assert depth<80,('unresolved',k,depth,str(hi))
  stack.extend([(2*k+1,depth+1),(2*k,depth+1)]);splits+=1
 m,e=answer.upper().man_exp()
 return dict(tile=index,m=int(m),e=int(e),counts=counts,splits=splits,maxdepth=maximum,seconds=time.time()-started)

def main():
 a=argparse.ArgumentParser();a.add_argument('certificate');a.add_argument('--precision',type=int,default=256);a.add_argument('--tol',default='1e-14');a.add_argument('--workers',type=int,default=2);a.add_argument('--target',default='0.253123274621');a.add_argument('--out',default='evidence/gauss.json');o=a.parse_args()
 begin=time.time();initialize(o.certificate,o.precision,o.tol)
 print('independent Gauss start; Mahler',MAHLER,'groups',len(GROUPS),'full circle',flush=True)
 acc=arb(0);tiles=[]
 with multiprocessing.Pool(o.workers,initializer=initialize,initargs=(o.certificate,o.precision,o.tol)) as pool:
  for result in pool.imap_unordered(integrate_tile,range(64),chunksize=1):
   acc+=arb(result['m'])*arb(2)**result['e'];tiles.append(result)
   print('tile',result['tile'],'seconds',round(result['seconds'],2),'counts',result['counts'],'depth',result['maxdepth'],flush=True)
 assert sorted(x['tile'] for x in tiles)==list(range(64))
 bound=MAHLER+acc
 result=dict(sha256=hashlib.sha256(Path(o.certificate).read_bytes()).hexdigest(),method='4-point Gauss / 8th derivative; quadratic root mapping; whole circle',precision=o.precision,tolerance=o.tol,target=o.target,passed=bool(bound.upper()<qball(o.target)),mahler=str(MAHLER),integral_upper=str(acc.upper()),bound_upper=str(bound.upper()),tiles=tiles,seconds=time.time()-begin)
 Path(o.out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='tiles'},indent=2),flush=True)
 raise SystemExit(0 if result['passed'] else 1)
if __name__=='__main__':main()
