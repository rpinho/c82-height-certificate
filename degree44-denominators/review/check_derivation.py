"""Symbolic checks of the fresh Gauss derivation and submitted Taylor formulas."""
from pathlib import Path
import json, math
from fractions import Fraction
import sympy as S
u,r=S.symbols('u r');v=S.log(u-r);forms={}
for k in range(1,9):
 v=S.cancel(u*S.diff(v,u))
 if k in [2,4,8]:forms[k]=v
assert S.cancel(forms[2]+r*u/(u-r)**2)==0
assert S.cancel(forms[4]+r*u*(u*u+4*u*r+r*r)/(u-r)**4)==0
E=u**6+120*r*u**5+1191*r*r*u**4+2416*r**3*u**3+1191*r**4*u*u+120*r**5*u+r**6
assert S.cancel(forms[8]+r*u*E/(u-r)**8)==0
sqrt30=S.sqrt(30);inner=S.sqrt((15-2*sqrt30)/35);outer=S.sqrt((15+2*sqrt30)/35)
nodes=[-outer,-inner,inner,outer];weights=[(18-sqrt30)/36,(18+sqrt30)/36,(18+sqrt30)/36,(18-sqrt30)/36]
for k in range(8):assert S.simplify(sum(w*x**k for w,x in zip(weights,nodes))-(S.Rational(2,k+1) if k%2==0 else 0))==0
error8=S.simplify(S.Rational(2,9)-sum(w*x**8 for w,x in zip(weights,nodes)))
c=Fraction(math.factorial(4)**4,9*math.factorial(8)**3)
assert S.simplify(error8-S.Rational(c.numerator,c.denominator)*2**9*math.factorial(8))==0
# Midpoint degree-4 term: integral x^4/4! on [-w/2,w/2] = w^5/1920.
x,w=S.symbols('x w',positive=True)
assert S.integrate(x**4/24,(x,-w/2,w/2))==w**5/1920
assert S.integrate(x,(x,0,w/2))*2==w*w/4
assert S.integrate(x*x/2,(x,-w/2,w/2))==w**3/24
out=dict(passed=True,derivatives_checked=[2,4,8],gauss_exact_through_degree=7,gauss_width9_constant=str(c),midpoint_width5_constant='1/1920',crossing_slope_constant='1/4',crossing_curvature_constant='1/24',formulas={k:str(v) for k,v in forms.items()})
Path('evidence/derivation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
