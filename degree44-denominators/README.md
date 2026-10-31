# C82 < 0.253123275 with degree-44 denominators

Unrefereed. This certifies an upper bound on the essential minimum of the Zhang-Zagier height, constant 82a in [Tao's optimization-constants table](https://teorth.github.io/optimizationproblems/constants/82a.html). No human expert has reviewed it line by line.

The rigorous upper endpoint is 0.25312327462072 from `verify_exact.py` at 256 bits. The separate integrator in `review/` gives 0.25312327462066. Strongest earlier bound I found: 0.2534350402, [PR #280](https://github.com/teorth/optimizationproblems/pull/280) by vidigoat, still open.

## What's here

- `certificate.json`: 10 monic denominator factors, D = 999,999,999,980, and 29 numerator factors with exact rational exponents. SHA-256 `6a24d45c18012bf7cc3f0ea05f7bcec685a43fd2f7e942f1ca9ce69fd558f056`.
- `PROOF.md`: why the certificate bounds the essential minimum, and how the integration is made rigorous.
- `METHOD.md`: how the factors were found.
- `verify_exact.py`, `check_verifier.py`: the checker and its controls. It imports no search code.
- `review/`: a separate AI review's exact algebra audit, root audit, derivative identities and whole-circle Gaussian integrator. None of it shares code with `verify_exact.py`.
- `receipts/`: output of every run below.

## Reproduce

Needs python-flint 0.9.0, and sympy for the controls. Run from this folder. The released copy of `verify_exact.py` passes with the same endpoint.

    mkdir -p evidence
    python verify_exact.py certificate.json --target 0.253123274621 --precision 256 --tol 1e-11 --out out.json
    python check_verifier.py
    python review/audit_algebra.py certificate.json algebra.json
    python review/audit_roots.py certificate.json roots.json
    python review/check_derivation.py
    python review/integrate_gauss.py certificate.json --workers 8 --precision 256 --tol 1e-14 --target 0.253123274621 --out gauss.json

The two earlier public checkers also pass on this certificate. [vidigoat's `vcheck.py`](https://github.com/vidigoat/optimization-records/tree/32898ac1dc0d75bf7f6d8a4e0a779eea629467a4/tao/certificates/c82) gives 0.25312327463 at tolerance 1e-9, and this repository's `independent_check.py` gives 0.2531251 at tolerance 1e-4. Receipts: `receipts/vcheck.json`, `receipts/independent-coarse.json`.

The search is not included and isn't needed to check the bound. It stopped on a time cap, so this says nothing about optimality.

## Credit

Doche's limit-point formula and factors, Grinsztajn's factors, and Épi's Ed5 and unequal multiplicities. The certificate format, joint LP and dual lattices come from the earlier certificates in this repository. The iterated IRLS/LLL numerator search and public baseline are vidigoat's. New here: the denominator search by complementary-measure monic IRLS/CVP, its eight degree-44 factors, and 23 numerator factors.

## AI use

Analysis and proof drafting by AI agents under my direction; I reviewed the results at summary level. AI agents also prepared the certificate, checkers and this text. This does not claim a line-by-line human verification. The separate review is also AI, so it isn't a second party.

Code MIT, as in the repository licence.
