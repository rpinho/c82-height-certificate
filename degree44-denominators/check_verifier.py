"""Sanity checks of mathematical derivative bounds and structural rejection."""
import json,copy,tempfile,pathlib
import sympy as sp
from verify_exact import check
from flint import arb

x,r=sp.symbols('x r');v=sp.log(x-r)
der=[]
for k in range(4):v=sp.simplify(x*sp.diff(v,x));der.append(v)
assert sp.simplify(der[1]+r*x/(x-r)**2)==0
assert sp.simplify(der[3]+r*x*(x*x+4*r*x+r*r)/(x-r)**4)==0
# Multiplication by (2*pi*i)^k gives the circle derivatives. In particular
# f''=+(2*pi)^2 r*x/(x-r)^2, |f''''| <= (2*pi)^4 |r|(1+4|r|+|r|^2)/d^4.
cert=dict(denominator=[dict(coeffs=[1,-10],multiplicity=1)],numerator=[])
path=pathlib.Path('evidence/verifier-sanity.json');path.write_text(json.dumps(cert))
a=check(str(path),'2.303','1e-9');b=check(str(path),'2.302','1e-9')
assert a['passed'] and not b['passed']
assert arb(a['bound_upper']).contains(arb(10).log())
# A nonzero, everywhere-positive integrand exercises midpoint curvature and
# its remainder. Here B = (1/2) log M(100*(T-100)) = log(100).
positive=dict(denominator=[dict(coeffs=[1,-10],multiplicity=1)],
 numerator=[dict(coeffs=[100,-10000],q_rational='1/2')])
positive_path=pathlib.Path('evidence/verifier-positive-sanity.json')
positive_path.write_text(json.dumps(positive))
positive_result=check(str(positive_path),'4.606','1e-9')
up=arb(positive_result['bound_upper']);exact=arb(100).log()
assert positive_result['passed'] and positive_result['counts'][1]>0
assert up>=exact and up<exact+arb('1e-8')
rejects=[]
mutations=[lambda c:c['denominator'][0].update(coeffs=[2,-10]),lambda c:c['denominator'][0].update(multiplicity=0),lambda c:c.update(numerator=[dict(coeffs=[1,-10],q_rational='1/2')]),lambda c:c.update(numerator=[dict(coeffs=[1,0],q_rational='-1')]),lambda c:c.update(numerator=[dict(coeffs=[1,0],q_rational='1')])]
for i,mutate in enumerate(mutations):
 c=copy.deepcopy(cert);mutate(c);path.write_text(json.dumps(c))
 try:check(str(path),'3','1e-9')
 except AssertionError:rejects.append(i)
assert len(rejects)==len(mutations)
path.write_text(json.dumps(cert))
json.dump(dict(derivative_identities=True,closed_form_log10=a,positive_integrand_log100=positive_result,below_log10_rejected=not b['passed'],structural_corruptions_rejected=rejects),open('evidence/verifier-sanity-results.json','w'),indent=2)
print('PASS: derivative identities, log(10), positive-integrand log(100), 5 structural corruptions')
