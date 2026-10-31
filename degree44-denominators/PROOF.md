# Why the certificate gives an upper bound on the essential minimum

All logarithms are natural. Integer coefficient arrays in the JSON are in descending order.
Write T(z)=z(1-z), D=deg Q, and

    Q(X) = product_j R_j(X)^m_j,
    A(X) = product_i S_i(X)^q_i,       q_i >= 0 rational.

The verifier checks that each R_j is monic up to sign, m_j is a positive
integer, deg A < D, and gcd(S_i,R_j)=1 whenever q_i>0. Here A denotes a
formal product with rational exponents; no branch of a complex power is used.
Only log|A|=sum q_i log|S_i| appears.

The upper-bound formula is

    B = [ log M(Q(T(z)))
          + integral_0^1 max(0, sum_i q_i log|S_i(T(exp(2 pi i t)))|
                                   - log|Q(T(exp(2 pi i t)))|) dt ] / D.

This is the logarithm of Doche 2001b, equation (9). The following direct
argument also verifies the infinite-sequence step without requiring a
numerical sequence of huge irreducible polynomials.

Take a common denominator b of the q_i and put

    U(z)=product_i S_i(T(z))^(b q_i),   W(z)=Q(T(z))^b.

These are coprime integer polynomials, deg U < deg W=2bD, and W is monic
up to sign. In particular U/W is not constant. For n>=1, set

    F_n(z)=U(z)^n-W(z)^n.

Its leading coefficient has modulus 1 and its degree is 2bDn. Factorization
using the n-th roots of unity and multiplicativity of Mahler measure give

    log M(F_n)/n = (1/n) sum_{k=1}^n log M(U-exp(2 pi i k/n) W).

The function of the phase on the right is continuous: polynomial degree
and leading-coefficient modulus stay fixed, roots vary continuously as a
multiset, and log M is the sum of log max(1,|root|). Thus this Riemann sum
converges. Jensen's formula, first in z and then in the phase, gives

    lim_n log M(F_n)/n
      = integral_0^1 log max(|U(exp(2 pi i t))|,|W(exp(2 pi i t))|) dt
      = b D B.

The exchange of integrals is valid. The positive part of log|U-yW| is
uniformly bounded on the two circles, while Jensen's phase integral is
log max(|U|,|W|), bounded below because U and W have no common zero.
Consequently the negative part is integrable as well.

For each prime p, use the integer polynomial

    G_p(z) = F_p(z)/(U(z)-W(z)).

It is monic up to sign, symmetric under z -> 1-z, and has degree
2bD(p-1). For distinct primes p and r, G_p and G_r are coprime. Indeed,
a common root cannot be a zero of W by coprimality of U,W. Its ratio U/W
would have to be both a nontrivial p-th and a nontrivial r-th root of
unity, which is impossible. (At U/W=1 the quotient evaluates to
p W^(p-1), so this is not a root.) This remains valid in the presence of
ramification or repeated roots of individual G_p.

As G_p is symmetric and monic up to sign, the average of h(alpha)+h(1-alpha)
over its roots, counted with multiplicities, equals

    2 log M(G_p)/deg G_p.

To justify this even for reducible G_p, factor into primitive irreducible
integer polynomials; the leading coefficients of all factors have modulus
1. Each factor contributes its degree times h(alpha), and the same applies
to h(1-alpha). Symmetry makes the two aggregate sums equal.

The displayed averages tend to B, since subtracting the fixed
log M(U-W) does not affect the normalized limit. Each G_p therefore has
an irreducible factor whose roots have Zhang-Zagier height no greater than
that average. Coprimality for distinct primes makes these algebraic numbers
distinct. Hence C82 <= B. No assumption that the search LP is optimal,
no measure-approximation theorem, and no irreducibility of the R_j is needed.

# Rigorous integration in verify_exact.py

The code composes the polynomials exactly in Z[z], encloses every composed
root (with multiplicity) using Arb, and encloses their Mahler measures by
Jensen's formula. No search or upstream verification code is imported.

For u(t)=exp(2*pi*i*t), a panel of width w and midpoint c satisfies
|u(t)-u(c)| <= pi*w. Triangle inequalities therefore enclose all root
distances on the panel. Signed weighted log-distance bounds establish
whether the normalized log-ratio g is negative throughout. Such panels
contribute zero, including any numerator zeros. Other panels require
positive lower distances from every root before Taylor estimates are used.

On crossing panels let M2 bound |g''| and h=w/2. Taylor's theorem and the
1-Lipschitz positive-part map give the upper estimate

    integral max(0,g) <= w*max(0,g(c)) + |g'(c)|*w^2/4 + M2*w^3/24.

The sign can also be established by g(c) +/- (|g'(c)|*h+M2*h^2/2).
Positive panels use the fourth-order remainder described below.

Dyadic recursion tiles [0,1/2] exactly; real coefficients give
g(1-t)=g(t), so doubling covers the entire circle. All evaluations,
comparisons and accumulated upper endpoints use Arb. Tolerance controls
subdivision and is not an assumed discretization error. The checker aborts
if an unresolved panel reaches depth 55; it does not assume denominator
roots on the unit circle can be handled. Every panel of our certificate
resolved successfully. The fields integral_upper and bound_upper enclose
computed upper estimates, rather than giving two-sided enclosures of the
true integral or exact construction value.

# Fourth-order bound and independent checks

`verify_exact.py` is written separately from the search and imports neither search code nor the older verifier. It checks the same exact polynomial/rational conditions, independently composes with T, and encloses complete composed root multisets in Arb. Dyadic panel endpoints are exact. Whole-panel root-distance triangle inequalities enclose g and establish sign.

On a separated panel, with u = exp(2*pi*i*t), each log-distance term has

    f''(t) = Re [(2*pi)^2 r*u/(u-r)^2],
    |f''''(t)| <= (2*pi)^4 |r| (1 + 4|r| + |r|^2) / d^4,

where d is the whole-panel lower bound on |u-r|. These follow by applying (u*d/du)^k to log(u-r); symbolic identities are independently checked in check_verifier.py. Absolute signed weights give a valid fourth derivative bound M4 for g. On positive panels, midpoint Taylor quadrature gives

    integral g <= width*g(mid) + width^3*g''(mid)/24 + M4*width^5/1920.

Odd terms integrate to zero. The previous second derivative / positive-part Lipschitz enclosure is used for crossing panels. When the integral contribution is accepted, its Arb upper endpoint is accumulated; doubling covers the whole circle by conjugation symmetry. The tolerance controls subdivision, not a guessed error bar. The final upper endpoint is compared to an exact rational target.

The first improved certificate replaces R2 with a new degree-44 denominator and reweights the remaining library. All proof conditions are checked on the frozen JSON, so the same infinite-sequence argument above applies directly.
