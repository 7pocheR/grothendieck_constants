# A real Grothendieck upper bound from a five-variable sign function

**Unpublished computer-assisted research.** With the numerical enclosures
specified below, our work gives
\[
 K_G^{\mathbb R}<1.779754412112.
 \tag{1}
\]
The argument uses one primary Gaussian coordinate, four auxiliary Gaussian
coordinates, signed tensor preprocessing, and a complete coefficient
estimate. The mathematical reduction is stated here in full. The finite
numerical enclosures and the exact consequences that can be checked from
this package have distinct scopes, described in the final section.

## The constant and the sign function

For a finite real matrix \(M=(M_{ij})\), define
\[
 S(M)=\max_{\sigma_i,\tau_j\in\{-1,1\}}
        \left|\sum_{ij}M_{ij}\sigma_i\tau_j\right|,
\qquad
 V(M)=\sup_{d\ge1}\max_{x_i,y_j\in S^{d-1}}
        \left|\sum_{ij}M_{ij}\langle x_i,y_j\rangle\right|.
\]
Then \(K_G^{\mathbb R}=\sup_{M\ne0}V(M)/S(M)\). The vector
dimension can be restricted to the number of row and column vectors,
since they span a space of at most that dimension.

All Gaussian variables below are independent standard real Gaussians
unless correlations are specified. Let
\[
 B(x)=\ell x+\sum_{l=1}^{6}b_l\sin(3lx/4),\qquad
 S(x)=\sum_{i=1}^{4}B(x_i),
\]
\[
 A(s)=\sum_{l=1}^{4}a_l\sin(2ls),\qquad U(s)=s+A(s),\qquad
 f(w,x)=\operatorname{sign}(w+U(S(x))).
 \tag{2}
\]
Every coefficient in (2), and every coefficient of the polynomial \(P\)
used below, is specified as an exact rational number in
`certificate/physical.json`. This file is part of the definition of the
candidate. Its SHA-256 is
`3077b5a8011868ff8a9200590bebb73cb41ad3ba10b09eb37860d5074c3b1653`.
The collective-amplitude field is zero. In particular, there is no
additional term in (2). The auxiliary polynomial is \(Q(z)=z\), and
\[
 P(z)=\sum_{k=1,3,\ldots,81}p_k z^k,
 \qquad q:=\sum_k|p_k|=\frac{999999999999}{10^{12}}<1.
 \tag{3}
\]
The signed value \(P(1)\) is different from \(q\). All estimates refer
to scalar radius one; (3) is already part of the preprocessing, not a
factor to apply again to the final bound.

Write \(\psi_j=\mathrm{He}_j/\sqrt{j!}\) for the orthonormal
probabilists' Hermite polynomials, and let
\[
 A_{j,\beta}=\mathbb E\left[f(W,X)\psi_j(W)
                         \prod_{i=1}^{4}\psi_{\beta_i}(X_i)\right],
 \qquad \beta\in\mathbb N^4.
\]
Conditioning on \(X\) shows that the zero set of the sign has Gaussian
measure zero. Parseval therefore gives \(\sum_{j,\beta}A_{j,\beta}^2=1\).
The scalar kernel is
\[
 H(z)=\frac\pi2\sum_{j,\beta}A_{j,\beta}^2P(z)^j(-z)^{|\beta|}
     =\sum_{n\ge1}c_n z^n.
 \tag{4}
\]
The minus sign reflects all four auxiliary coordinates, while the primary
coordinate is not reflected. Global oddness of \(f\) implies that
\(A_{j,\beta}=0\) when \(j+|\beta|\) is even. Thus (4) has real
coefficients and only positive odd powers. Moreover,
\[
 \sum_n|c_n|\le\frac\pi2\sum_{j,\beta}A_{j,\beta}^2q^j\le\frac\pi2.
 \tag{5}
\]
Absolute convergence justifies coefficient rearrangements and boundary
evaluation, including \(z=\pm1\).

## The matrix implication

For each unit vector \(x\), form the unit vectors
\[
 L_P(x)=\bigoplus_k\sqrt{|p_k|}\,x^{\otimes k}
                  \oplus\sqrt{1-q}\,e_L,
\quad
 R_P(y)=\bigoplus_k\operatorname{sign}(p_k)\sqrt{|p_k|}\,y^{\otimes k}
                  \oplus\sqrt{1-q}\,e_R,
\]
where the two padding directions are orthogonal to one another and to all
tensor summands. Their cross inner product is
\(P(\langle x,y\rangle)\). Use one common Gaussian field on the
primary tensor space and four independent common Gaussian fields on the
original vector space. Negate the right projections of the latter four
fields. Applying (2) on each side produces jointly defined signs for the
entire matrix, with
\[
 \mathbb E[\sigma_i\tau_j]=\frac2\pi H(\langle x_i,y_j\rangle).
 \tag{6}
\]
The Hermite covariance identity proves (6) for interior correlations;
absolute convergence and Gaussian continuity extend it to the endpoints.
This construction supplies the joint distributions, rather than relying
on an entrywise covariance test.

Choose vectors attaining \(V(M)\), with positive global orientation.
Compactness of their finite Gram matrices gives attainment, or one can
use a sequence tending to the supremum. Every tensor family is again a
family of unit vectors, so
\(\left|\sum M_{ij}\langle x_i,y_j\rangle^n\right|\le V(M)\).
Equations (5) and (6) give
\[
 S(M)\ge\mathbb E\sum_{ij}M_{ij}\sigma_i\tau_j
 \ge \frac2\pi\left(c_1-\sum_{n\ge3}|c_n|\right)V(M).
 \tag{7}
\]
Consequently, if
\(\Gamma\le c_1-\sum_{n\ge3}|c_n|\) and \(\Gamma>0\), then
\(K_G^{\mathbb R}\le\pi/(2\Gamma)\). This applies to every finite
real matrix and every real vector dimension. It uses neither an
inverse-series theorem nor a sign assumption on the nonlinear coefficients.

## The finite coefficient calculation

Set \(N=403\), \(M=31\), \(L=40\), \(\Delta=2\pi/L\), and use
the positive Fourier nodes \(u_n=n\Delta\), \(1\le n\le F=768\).
Only indices with \(j+|\beta|\le403\) can contribute to (4) through
degree 403, since \(P(0)=0\) and its least degree is one. Work throughout
in the quotient by \(z^{404}\), retaining all 404 coordinates.

Let \(\phi(u)=e^{-u^2/2}/\sqrt{2\pi}\). The conditional primary
coefficients have the Fourier representation
\[
 q_j(s):=\mathbb E_W[\operatorname{sign}(W+s)\psi_j(W)]
 =\sqrt{2/\pi}\,i^{j-1}\int_{\mathbb R}
        V_j(u)e^{ius}\,du,
 \quad V_j(u)=\frac{\phi(u)u^{j-1}}{\sqrt{j!}}.
 \tag{8}
\]
For \(j=0\), interpret \(V_0=\phi(u)/u\) as a symmetric principal
value; \(q_0(s)=\operatorname{erf}(s/\sqrt2)\) fixes its constant.
Expand the outer phase through degree 31. With
\[
 R(\zeta)=\sum_{l=1}^{4}\frac{a_l}{2}(\zeta^l-\zeta^{-l}),
 \qquad R(\zeta)^m/m!=\sum_k r_{m,k}\zeta^k,
\]
the resulting weight is
\[
 \widetilde V_j(u)=\sum_{m=0}^{31}\sum_k
        r_{m,k}(u-2k)^m V_j(u-2k).
 \tag{9}
\]
All shifts and their signed coefficients are included. Only the term
\(m=j=0\) has a principal-value singularity. For \(m\ge1\), cancel
the corresponding power of \(u-2k\) before numerical evaluation.

Define the real functions \(p_k\) by
\(\mathbb E[\psi_k(X)e^{iuB(X)}]=i^k p_k(u)\). For the finite
calculation they are approximated by the direct signed spatial lattice
with spacing \(1/256\) and indices \(-5120,\ldots,5120\), using one
normalization factor \(1/\sqrt{k!}\). The Gaussian integral and omitted
lattice allowances are accounted for separately. Set
\[
 K_{uv}(z)=\sum_{k=0}^{403}\widetilde p_k(u)\widetilde p_k(v)z^k,
 \quad W_{uv}(z)=\sum_{j=0}^{403}\widetilde V_j(u)\widetilde V_j(v)z^j.
\]
The positive-frequency contribution is
\[
 4\Delta^2\operatorname{Odd}_{\le403}
     \sum_{u,v\in\{u_1,\ldots,u_{768}\}}
       W_{uv}(P(z))K_{uv}(-z)^4.
 \tag{10}
\]
The fourth power counts all ordered auxiliary multiindices, so no additional
multinomial or factorial is inserted. The scalar normalization follows
from \((\pi/2)(\sqrt{2/\pi})^2=1\). A triangular summation has
295,296 pairs, 768 diagonal pairs, and 589,824 terms with multiplicity.
All sums remain signed until the final coefficient objective.

There are two zero-frequency corrections. For odd \(j,k\ge1\), let
\[
 a_j^{(0)}=2\Delta\sum_{n=1}^{768}
      \widetilde V_j(u_n)\widetilde p_0(u_n)^4,
 \qquad z_j=\Delta\widetilde V_j(0),
\]
\[
 b_k^{(0)}=2\Delta\sum_{n=1}^{768}
       \widetilde V_0(u_n)\widetilde p_k(u_n)\widetilde p_0(u_n)^3,
\]
\[
 t_k=\Delta\phi(0)\left[\ell\mathbf1_{k=1}
       +\sum_{l=1}^{6}b_l(3l/4)^k e^{-(3l/4)^2/2}/\sqrt{k!}\right].
\]
Add to (10)
\[
 \sum_{j\text{ odd}}(2a_j^{(0)}z_j+z_j^2)P(z)^j
    -4\sum_{k\text{ odd}}(2b_k^{(0)}t_k+t_k^2)z^k.
 \tag{11}
\]
The second term comes from the principal-value product limit. Exactly one
auxiliary index can be positive at zero frequency; there are four choices,
and reflection gives the minus sign. Assigning a finite value to
\(\widetilde V_0(0)\) would be incorrect. Both cross terms and squares
in (11) are required, with no additional \(4\Delta^2\) multiplier.

The analytic error of (8)--(11) includes the phase Taylor remainder,
every nonzero Fourier alias, all frequencies above 768, the complete
spatial strip error, and the omitted spatial lattice. An additional
phase-bin allowance is retained conservatively even though direct spatial
summation needs no bin approximation. The [head estimates](head_estimates.md)
give every formula and its hypotheses; `certificate/head_error_input.json`
gives all parameters and the 404 individual frequency choices.

If \(\delta\) bounds the Euclidean error in the Hermite coefficient
array, Bessel's inequality, the fourfold tensor estimate, and
Cauchy--Schwarz give
\[
 \sum_{j,\beta}|A_{j,\beta}^2-\widetilde A_{j,\beta}^2|
 \le \delta(2+\delta).
\]
Composition by \(P\) and \(-z\), followed by restriction of degrees,
contracts coefficient absolute norm. Thus
\[
 E=\frac\pi2\delta(2+\delta)
   \le E_{\rm cert}\simeq1.5978022940933675\cdot10^{-10}
 \tag{12}
\]
is a complete scalar error allowance. The exact upper bound is stored in
`certificate/scalar_bounds.json`. The number of auxiliary multiindices
does not produce another multiplicity: Bessel bounds the whole array.

For the 404 finite coefficient intervals \([l_n,u_n]\) in
`certificate/head_coefficients.json`, form
\[
 g=l_1-\sum_{n=3,5,\ldots,403}\max(|l_n|,|u_n|)-E_{\rm cert}.
 \tag{13}
\]
This yields \(g\simeq0.8826334911000137\). All arithmetic interval
widths and outward enlargements are retained in the endpoints. The
coefficient objective is 1-Lipschitz for the coefficient absolute norm,
which proves that (12) is subtracted exactly once in (13).

## Complete localization of the infinite tail

The following argument concerns the true kernel (4), not its finite
approximation. First, its boundary derivatives exist to every fixed order.
Put \(L(x)=U(S(x))\). Every positive derivative of \(L\) is bounded.
For \(0<r<1\), smoothing only the primary coordinate gives
\[
 G_r(w,x)=\operatorname{erf}
       \left(\frac{\sqrt r\,w+L(x)}{\sqrt{2(1-r)}}\right).
\]
Every derivative of fixed positive order is a finite sum of bounded
boundary derivatives times a polynomial times a real Gaussian exponential.
It has finite Gaussian Sobolev norm. Its Hermite coefficients are
\(r^{j/2}A_{j,\beta}\). Summing the squared norms of its mixed
derivatives gives
\[
 r^p S_{p,q}(r),\qquad
 S_{p,q}(r)=\sum_{j,\beta}(j)_p(|\beta|)_q
                           r^{j-p}A_{j,\beta}^2<\infty.
 \tag{14}
\]
Here the sums over auxiliary derivatives are ordered, producing the
falling factorial \((|\beta|)_q\). Divide by \(r^p\) for \(r>0\);
the nonnegative series gives the limit at zero. Apply (14) at any radius
strictly between \(q\) in (3) and one. It dominates the differentiated
kernel uniformly for \(|P(z)|\le q\) and \(|z|\le1\), and proves
the required boundary regularity. No differentiation of the unsmoothed
sign is assumed.

Near real auxiliary correlations, (14) and the finite chain rule bound
the derivatives through order four. The computation retains the signed
conditional Fourier expansions before taking their Gaussian squared
norms. The supplied source files specify the conditional terms, analytic
remainders, nonnegative radial interpolation and complete circle
evaluation. The underlying mixed-energy calculation comprises 126
enclosures and yields 14 mixed derivative bounds on each of 2,048 panels.
These numerical enclosures are an explicit component of the result.

Away from real auxiliary correlations, higher derivatives can be bounded
directly from the full Gaussian density. For one correlated real Gaussian
pair with complex correlation \(a\), the absolute density has mass
\[
 d_0(a)=\sqrt{\frac{|1-a^2|}{1-(\operatorname{Re}a)^2}}
\]
and, after normalization in the plus/minus coordinates, independent real
Gaussian variances \(\sigma_{a,\pm}=|1\pm a|^2/(1\pm\operatorname{Re}a)\).
Use one primary pair with \(a=P(z)\), and four auxiliary pairs with
\(b=-z\). Denominators are required to be strictly positive on each
panel where this estimate is used. Also \(1-a^2\) avoids the nonpositive
real axis when \(|\operatorname{Re}a|<1\), giving the consistent
Gaussian square-root branch. Absolute integrability permits
differentiation under the integral.

For \(a(u)=P(ze^u)\), \(b(u)=-ze^u\), let \(L_j\) be the \(j\)-th
derivative at zero of the log density
\[
 -\tfrac12\log(1-a(u)^2)-2\log(1-b(u)^2)
 -\sum_{\epsilon=\pm}\left(
 \frac{S_{a,\epsilon}}{2(1+\epsilon a(u))}
 +\frac{S_{b,\epsilon}}{2(1+\epsilon b(u))}\right)+\mathrm{constant}.
\]
Normalize the four squared-coordinate sums by their variances. They are
independent chi-square variables of degrees \((1,1,4,4)\). Each \(L_j\)
is affine in these variables. Define the complete derivative polynomial
\[
 R_0=1,\qquad
 R_n=\sum_{k=1}^{n}\binom{n-1}{k-1}R_{n-k}L_k.
\]
Since \(|f|\le1\), Cauchy--Schwarz gives
\[
 |(z\partial_z)^nH(z)|\le
 \frac\pi2 d_0(a)d_0(b)^4\sqrt{\mathbb E|R_n|^2}.
 \tag{15}
\]
There is no missing auxiliary or determinant factor. Angular derivatives
multiply (15) by \(i^n\), which has modulus one.

The expectation in (15) is a finite polynomial norm. For a chi-square
variable \(Y=2T\), where \(T\) has Gamma shape \(\alpha\),
\[
 Y^k=2^k k!\sum_{n=0}^{k}(-1)^n
       \frac{(\alpha+n)_{k-n}}{(k-n)!}L_n^{(\alpha-1)}(T),
 \qquad \mathbb E|L_n^{(\alpha-1)}(T)|^2=(\alpha)_n/n!.
\]
Expansion in the product Laguerre basis with shapes
\((1/2,1/2,2,2)\) therefore evaluates the squared norm with positive
rational weights. Interval coefficients bound it on a complete angular
panel. This includes the varying Gaussian variances in the coefficients;
the normalized chi-square variables have a fixed distribution.

Choose a real, even, \(\pi\)-periodic cutoff \(\chi\). In the first
quadrant put
\[
 a=592/2048=37/128,\quad b=1248/2048=39/64,\quad
 t=(2\theta/\pi-a)/(b-a).
\]
Use \(\chi=1\) below the transition, \(\chi=0\) above it, and
\[
 \chi(t)=1-12012\int_0^t u^6(1-u)^6\,du
 \tag{16}
\]
inside. Derivatives one through six vanish at both joins. The even,
periodic extension is \(C^6\). The degree-13 Bernstein coefficients of
(16) are seven ones followed by seven zeros; derivatives and rational
subinterval restriction give exact bounds for each cutoff derivative.

Set \(F_0(\theta)=\chi(\theta)H(e^{i\theta})\) and
\(F_1(\theta)=(1-\chi(\theta))H(e^{i\theta})\). Both have only odd
Fourier indices, possibly negative. With normalized circle measure,
Parseval and Cauchy--Schwarz imply
\[
 \sum_{n>403,\ n\text{ odd}}|c_n|
 \le\frac{\|F_0^{(4)}\|_2}{\sqrt{14}\,403^{7/2}}
     +\frac{\|F_1^{(6)}\|_2}{\sqrt{22}\,403^{11/2}}.
 \tag{17}
\]
Indeed \(\sum_{n>N,\ n\text{ odd}}n^{-2m}\le
[2(2m-1)N^{2m-1}]^{-1}\) for odd \(N\): integrate each decreasing
summand over the length-two interval immediately before its index.
Parseval already includes negative Fourier indices; they introduce no
extra factor of two in (17).

Let \(M_j\) bound \(|\partial_\theta^j H|\) uniformly on a panel.
Use \(M_0=\pi/2\), the complete mixed estimates for \(j\le4\), and
(15) wherever applicable. The full Leibniz bounds are
\[
 L_4=\sum_{k=0}^{4}\binom4k|\chi^{(k)}|M_{4-k},\qquad
 L_6=|1-\chi|M_6+
          \sum_{k=1}^{6}\binom6k|\chi^{(k)}|M_{6-k}.
 \tag{18}
\]
No fifth or sixth derivative bound is needed where its multiplier is
identically zero. In particular the density estimate is not used at
the real auxiliary singularities. The first-quadrant panels are
\([i\pi/4096,(i+1)\pi/4096]\), \(0\le i<2048\). Fourfold symmetry
and normalized circle measure give
\[
 \|F_0^{(4)}\|_2^2\le Q_4:=\frac1{2048}\sum_i L_{4,i}^2,
 \qquad
 \|F_1^{(6)}\|_2^2\le Q_6:=\frac1{2048}\sum_i L_{6,i}^2.
 \tag{19}
\]
The exact upper values of \(Q_4,Q_6\) and the resulting upper bound
\[
 T\ge\sqrt{Q_4/(14\,403^7)}+\sqrt{Q_6/(22\,403^{11})}
   \simeq0.00004187278135447257
 \tag{20}
\]
are in `certificate/scalar_bounds.json`. The rational square-root test
in `verify.py` checks (20) directly. Equations (17)--(20) cover every
omitted coefficient starting at degree 405. The complete-circle
fourth-derivative estimate is an ingredient of this construction; it is
not another coefficient tail to subtract. The errors in derivative
enclosures are already included in \(Q_4,Q_6\).

## Exact endpoint and numerical premises

Using the exact values in (12), (13), and (20), integer/rational arithmetic
gives
\[
 g-T>\Gamma_0:=\frac{17651832366373}{20000000000000}
              =0.88259161831865.
 \tag{21}
\]
The remaining rational margin is about \(9.1763\cdot10^{-15}\).
Alternating-series bounds in Machin's identity
\(\pi=16\arctan(1/5)-4\arctan(1/239)\) give
\[
 \pi<2\Gamma_0(1.779754412112).
\]
Together with (7), this proves (1) from the numerical enclosures.

The finite coefficient calculation uses outward real-ball arithmetic
with 512-bit primary and characteristic primitives and complete signed
pair sums. Its saved consequences and both zero corrections were
recomputed with exact rational aggregation and 640-bit arithmetic. The
last computation does not independently rerun all transient pair products.
The localized tail uses the complete mixed-energy enclosures and 512-bit
directional-density panel calculations. The software versions and native
library identities are recorded in `certificate/software.json`.

The standard-library verifier included here checks all supplied coefficient
intervals, the final error accounting, the last square-root inequality and
the endpoint. It does not establish that the supplied intervals enclose the
primitive Gaussian integrals, independently repeat every pair product, or
recompute (19) from all physical panel data. Those remain computer-assisted
execution premises of the result. `certificate/external_dependencies.json`
identifies the retained bulk data; hashes establish identity, not validity.
The [reproduction specification](reproduction.md) explains the additional
calculations and data needed to examine these premises. The small bundle
is self-contained for the mathematical reduction and endpoint arithmetic,
but not for full numerical reproduction.
