# A computer-assisted lower bound for the real Grothendieck constant

For a nonzero finite real matrix $A=(A_{ij})$, define
$$
 S_{\mathbb R}(A)=\max_{\epsilon_i,\delta_j\in\{-1,1\}}
                 \left|\sum_{i,j}A_{ij}\epsilon_i\delta_j\right|,
 \qquad
 V_{\mathbb R}(A)=\sup_{\|u_i\|=\|v_j\|=1}
                 \left|\sum_{i,j}A_{ij}\langle u_i,v_j\rangle\right|.
$$
The vector supremum allows any common finite-dimensional real Euclidean
space. The real Grothendieck constant is
$K_G^{\mathbb R}=\sup_{A\ne0}V_{\mathbb R}(A)/S_{\mathbb R}(A)$,
where all finite matrix sizes are allowed.

**Theorem.**
$$
 K_G^{\mathbb R}\ge\frac{1625\pi}{2917}>\frac74.
$$
A particular finite matrix of order $(40\cdot2^{40})^{100000}+1$,
defined by (25)--(26), has ratio strictly greater than $7/4$.

Put $a=13/8$, $\nu=\sqrt{2/\pi}$, and $C=2917/2000$.
The proof first establishes
$$
 |a\langle P_1f,P_1g\rangle-\langle P_3f,P_3g\rangle|
 \le \nu^2 C
 \tag{1}
$$
for every finite Gaussian dimension and every pair of measurable real
sign functions. Sections 1--5 prove (1) using an analytic circle
inequality and two finite certificates for ternary functions. Section 6
constructs finite real matrices from (1), with explicit approximation
errors. The additional argument in
[Joint and individual Gaussian noise constraints](NOISE_MATHEMATICS.md)
gives a second proof of the small-first-moment estimate. The proof here
uses only the circle inequality and the two ternary certificates.

## 1. Normalization and midpoints

Use standard Gaussian probability $\gamma_n$. Define the probabilists'
Hermite polynomials by
$e^{ts-t^2/2}=\sum_{j\ge0}H_j(s)t^j/j!$, and write
$\psi_j=H_j/\sqrt{j!}$. For multi-indices, use the product basis
$\psi_\alpha$, and let $P_j$ be the orthogonal projection onto
total degree $j$. The generating identity
$$
 \mathbb E e^{u\cdot X-|u|^2/2}e^{v\cdot X-|v|^2/2}=e^{u\cdot v}
$$
proves orthonormality. Comparison of homogeneous degrees also proves
orthogonal invariance of each total-degree space. Completeness follows,
for example, because a function orthogonal to all polynomials has zero
Gaussian-weighted Laplace transform, and hence zero Fourier transform:
Cauchy--Schwarz makes that transform entire, and uniqueness of the
Fourier transform gives the conclusion.

For arbitrary measurable signs $f,g$, define the original midpoints
$$
 h_0=(f+g)/2,\qquad k_0=(f-g)/2.
$$
They take values in $\{-1,0,1\}$, satisfy $h_0k_0=0$ and
$h_0^2+k_0^2=1$, and obey the exact identity
$$
 L_a=a(\|P_1h_0\|_2^2-\|P_1k_0\|_2^2)
       -\|P_3h_0\|_2^2+\|P_3k_0\|_2^2.                 \tag{2}
$$
For $n\ge1$, choose a unit coordinate $S$ in the direction of
the first coefficient vector of $h_0$, choosing any direction if
that vector is zero. Set
$$
 x=\|P_1h_0\|_2/\nu,\qquad \mathbb E Sh_0=\nu x.
$$
Since $|h_0|\le1$, $\mathbb E Sh_0\le\mathbb E|S|=\nu$,
so $0\le x\le1$. Dimension zero has $L_a=0$ and needs no
coordinate choice.

For the noise argument, define $q_o(X)=(q(X)-q(-X))/2$ and
$h=(h_0)_o$, $k=(k_0)_o$. Odd-degree Hermite polynomials change
sign under global inversion, so
$$
 P_jq_o=P_jq\quad(j\text{ odd}).                         \tag{3}
$$
Consequently the first coefficient vector, its orientation, and $x$
are exactly the same for $h$ and $h_0$. The cubic projections
and objective are also unchanged. However, $h,k$ need not be ternary.
They satisfy
$$
 |h|+|k|=\max(|f_o|,|g_o|)\le1.
$$
The ternary inequalities below are applied to $k_0$, and the noise
inequalities in the accompanying note are applied to $h,k$.
These are different functions with the same required odd projections.

## 2. An analytic bound for the first and third Gaussian projections

All circle integrals in this section use probability measure
$d\theta/(2\pi)$. Let $E_j$ denote real circle frequency $j$.
Every real trigonometric polynomial $p$ with frequencies one and
three satisfies
$$
 \left(\int|p|\right)^2
 \le \frac9{\pi^2}\|E_1p\|_2^2+
       \frac8{\pi^2}\|E_3p\|_2^2.                       \tag{4}
$$

Here is a proof covering every coefficient direction. A nonzero such
polynomial is antiperiodic with period $\pi$. Multiplication of its
Laurent polynomial by $z^3$ gives a nonzero polynomial of degree at
most six. There are at most six distinct circle zeros. Its sign
$\sigma$ therefore has either one or three sign changes on a
half-circle whose endpoints avoid zeros. Zeros of even multiplicity do
not affect the sign. Antiperiodicity makes the number of changes odd.

For one change, $\sigma$ is a translate, possibly negated, of
$\operatorname{sgn}\sin\theta$, and
$\|E_j\sigma\|_2^2=8/(\pi^2j^2)$ for odd $j$.
For three changes at $r_1<r_2<r_3<r_1+\pi$, the jumps alternate.
Distributional integration by parts gives
$$
 \|E_j\sigma\|_2^2=
 \frac8{\pi^2j^2}|e^{ijr_1}-e^{ijr_2}+e^{ijr_3}|^2,
 \qquad j\text{ odd}.                                  \tag{5}
$$
In this normalization the complex Fourier coefficient has modulus
$2|e^{ijr_1}-e^{ijr_2}+e^{ijr_3}|/(\pi j)$, and the real
projection energy is twice its square.

Set $a_0=r_2-r_1$, $b_0=r_3-r_2$,
$u=\cos((a_0+b_0)/2)$, $v=\cos((a_0-b_0)/2)$.
Then $0<u\le v\le1$. With $T_3(t)=4t^3-3t$, the squared
moduli at frequencies one and three are respectively
$$
 Z_1=1+4u^2-4uv,\qquad
 Z_3=1+4T_3(u)^2-4T_3(u)T_3(v).
$$
Direct multiplication gives
$$
 9-8Z_1-Z_3
 =4u(v-u)\left[8+(4u^2-3)\{4(u^2+uv+v^2)-3\}\right]. \tag{6}
$$
The bracket is nonnegative. For $u\ge\sqrt3/2$, both factors of
its product are nonnegative. For $u\le\sqrt3/2$, it equals
$$
 (2u-1)^2(4u^2+8u+5)
 +4(1-v)(1+v+u)(3-4u^2),                                \tag{7}
$$
whose terms are nonnegative. Thus, in both sign-change cases,
$$
 \frac89\|E_1\sigma\|_2^2+\|E_3\sigma\|_2^2
 \le\frac8{\pi^2}.
$$
Weighted Cauchy--Schwarz applied to $\langle\sigma,p\rangle$
proves (4). The zero polynomial is immediate. The polynomial identities
(6)--(7) also cover limiting, coincident-zero configurations by
continuity.

Let $F_1,F_3$ be arbitrary real polynomials in Gaussian chaoses one
and three in any fixed finite dimension. For independent standard
Gaussian vectors $X,Y$, put
$$
 p_{X,Y}(\theta)=(F_1+F_3)(X\cos\theta+Y\sin\theta).
$$
This polynomial has only the circle frequencies one and three. For any
angles its covariance is
$$
 \mathbb E p(\theta)p(\varphi)=
 \|F_1\|_2^2\cos(\theta-\varphi)
 +\|F_3\|_2^2\cos^3(\theta-\varphi).                    \tag{8}
$$
To verify (8), the two Gaussian arguments have marginal covariance
$I$ and cross covariance $\cos(\theta-\varphi)I$.
The generating identity gives
$\mathbb E\psi_\alpha(U)\psi_\beta(V)
=\mathbf1_{\alpha=\beta}\rho^{|\alpha|}$, including negative
$\rho$. The two polynomials have finite expansions, so this proves
(8) without an infinite-series interchange. The joint covariance is
positive semidefinite for every angle, including correlations $\pm1$.

Since $\cos^3 t=(3\cos t+\cos3t)/4$, circle orthogonality gives
$$
 \mathbb E\|E_1p\|_2^2=\|F_1\|_2^2+\tfrac34\|F_3\|_2^2,
 \qquad \mathbb E\|E_3p\|_2^2=\tfrac14\|F_3\|_2^2.
$$
Fubini and Cauchy--Schwarz are justified by finite Gaussian polynomial
second moments. Averaging (4), then applying Jensen to the square root,
proves
$$
 \mathbb E|F_1+F_3|
 \le \left(\frac9{\pi^2}\|F_1\|_2^2+
             \frac{35}{4\pi^2}\|F_3\|_2^2\right)^{1/2}. \tag{9}
$$
For any measurable $|f|\le1$, put
$A=9/\pi^2$, $B=35/(4\pi^2)$,
$F_1=P_1f/A$, $F_3=P_3f/B$. If
$Z=\|P_1f\|_2^2/A+\|P_3f\|_2^2/B$, then
$Z=\mathbb E f(F_1+F_3)\le\mathbb E|F_1+F_3|\le\sqrt Z$.
This includes $Z=0$ and proves
$$
 \frac{\pi^2}{9}\|P_1f\|_2^2+
 \frac{4\pi^2}{35}\|P_3f\|_2^2\le1.                    \tag{10}
$$
The parallelogram identity applied separately to $f,g$ yields
$$
 \|P_3h_0\|_2^2+\|P_3k_0\|_2^2
 +\frac{35}{36}(\|P_1h_0\|_2^2+\|P_1k_0\|_2^2)
 \le\frac{35}{4\pi^2}.
$$
Substitution in (2), retaining the signs of every coefficient, gives
$$
 \frac{L_a}{\nu^2}
 \le\frac{35}{8\pi}+\frac{47}{72}x^2
 -\frac{2\|P_3h_0\|_2^2}{\nu^2}
 -\left(a+\frac{35}{36}\right)\frac{\|P_1k_0\|_2^2}{\nu^2}.
                                                               \tag{11}
$$
The last two terms are nonpositive. Consequently (11) proves
$L_a/\nu^2<C$ throughout $0\le x\le11/50$. For example,
using $\pi>3141592653589793/10^{15}$, the gap at the right endpoint
is greater than
$$
 \frac{19396081925963203499}{565486677646162740000}>0.
$$
The coefficient of $x^2$ is positive, so the same endpoint check
covers the entire interval.

## 3. The two ternary inequalities and their transverse implication

For a measurable ternary function $u$ on the Gaussian line, set
$m_j=\mathbb E u\psi_j$ and
$\beta=\mathbb E|S|u(S)^2$. Define
$$
 Q_p(m)=(m_0,m_2)M_e(m_0,m_2)^T+(m_1,m_3)M_o(m_1,m_3)^T,
$$
where
$$
 M_e=\begin{pmatrix}p&r\\r&q\end{pmatrix},
 \quad q=\frac{p+a}{a+1},\quad r=\sqrt{(p-1)q},
 \qquad
 M_o=\begin{pmatrix}139/118&21\sqrt6/118\\
                     21\sqrt6/118&63/59\end{pmatrix}.
$$
The two targets are the following fully quantified inequalities:
$$
\begin{array}{ll}
 p=21/20:& Q_p(m)\le\tfrac32\nu\beta+\tfrac{26}{125}\nu^2,
                                                        \quad(T_A)\\
 p=13/10:& Q_p(m)\le\tfrac{11}{5}\nu\beta+\tfrac1{6000}\nu^2.
                                                        \quad(T_B)
\end{array}
$$
Their $(q,r^2)$ values are $(107/105,107/2100)$ and
$(39/35,117/350)$, respectively.

For each parity block $M$, all three differences
$$
 M,\quad M-\operatorname{diag}(1,0),\quad
 M-\operatorname{diag}(-a,1)                             \tag{12}
$$
are positive semidefinite, while $M$ and $2I-M$ are positive
definite. These are rational determinant tests: for entries
$P,R,S$, use the diagonal tests and determinants
$(P-d_0)(S-d_1)-R^2$, with
$(d_0,d_1)=(0,0),(1,0),(-13/8,1)$, and use
$(2-P)(2-S)-R^2>0$. All entries other than the sign of $R$
are rational; the specified $R$ is positive.
In fact the odd block has the same parameterization, with
$p=139/118$, $q=63/59$. In all three blocks, the determinants
in (12) are respectively $q,0,0$. The determinants of $2I-M$
are $37/42$, $2/7$, and $34/59$ for the A-even, B-even,
and common odd blocks. The required diagonal inequalities are strict
for $M$ and $2I-M$, and nonnegative for the two differences.

Write $X=(S,T)$, and let
$K_j(T)=\mathbb E_S k_0(S,T)\psi_j(S)$.
Fubini and conditional Cauchy--Schwarz make these $L^2$ functions.
If $\Pi_d$ denotes transverse degree $d$, then
$$
\begin{split}
 \|P_3k_0\|_2^2-a\|P_1k_0\|_2^2
 ={}&\|\Pi_3K_0\|_2^2+\|\Pi_2K_1\|_2^2
   +\|\Pi_1K_2\|_2^2+\|\Pi_0K_3\|_2^2\\
 &-a\|\Pi_1K_0\|_2^2-a\|\Pi_0K_1\|_2^2
 \le{}&\mathbb E_TQ_p(K(T)).                            \tag{13}
\end{split}
$$
For the even pair $(K_0,K_2)$, the left quadratic form is
$\operatorname{diag}(1,0)$ on transverse degree three and
$\operatorname{diag}(-a,1)$ on degree one. For the odd pair
$(K_1,K_3)$, these are the forms on degrees two and zero,
respectively. On each orthogonal complement it is zero. Thus (12)
proves (13) on the entire space. This argument may be made with those
finite orthogonal decompositions alone; no uncomputed transverse tail
is discarded without a positive semidefinite bound. Dimension one has
a singleton transverse space and satisfies the same argument.

Each section $k_0(\cdot,T)$ is ternary almost everywhere. Let
$\beta(T)=\mathbb E_S|S|k_0(S,T)^2$. The budget is
$$
 \mathbb E_T\beta(T)\le\nu(1-x),                        \tag{14}
$$
because $\nu x=\mathbb E Sh_0\le\mathbb E|S||h_0|
=\nu-\mathbb E|S|k_0^2$. Therefore either target, with its own
positive $t$ and $\kappa$, implies
$$
 L_a/\nu^2\le ax^2-\|P_3h_0\|_2^2/\nu^2
                    +t(1-x)+\kappa.                   \tag{15}
$$

## 4. Finite certificates for both ternary inequalities

Order the coordinates as $0,1,2,3$, and factor the block matrix of
$Q_p$ as $LL^T$, with the positive Cholesky diagonal and positive
off-diagonal entries in the two parity blocks. For $v\in\mathbb R^4$,
put
$$
 F(v)=\nu^{-1}\mathbb E\max\{t|S|,2|\psi(S)^TLv|\}-|v|^2.
$$
Completing the square gives the exact identity
$$
 \sup_{u\text{ ternary}}(Q_p(m)/\nu^2-t\beta/\nu)
 =\sup_{v\in\mathbb R^4}F(v)-t.                         \tag{16}
$$
Indeed, $Q_p(m)/\nu^2=\sup_v\{2m^TLv/\nu-|v|^2\}$.
The two operations are suprema, so their interchange is valid. At
fixed $v$, the pointwise maximum over $u\in\{-1,0,1\}$ is
$(2|\psi^TLv|-t|S|)_+$, attained by a measurable threshold rule.
All functions involved are Gaussian-integrable polynomials. No minimax
theorem or exchange of an infimum and a supremum occurs here.

Because $0<Q_p<2I$, for each unit $w$ the expectation part of
$F(rw)$ is Lipschitz in $r$ with constant at most
$2\sqrt2/\nu=2\sqrt\pi$. Hence
$$
 F(r_2w)-F(r_1w)
 \le(r_2-r_1)(2\sqrt\pi-r_1-r_2)\le0
 \quad(r_2\ge r_1\ge\sqrt\pi).
$$
The supremum occurs in the ball of radius $\sqrt\pi<2$.
Changing the sign of the even parity block or of the odd parity block
does not change $F$: use reflection of $S$ and global negation
of the polynomial. Thus every needed value has a representative in
$$
 D=[0,2]\times[0,2]\times[-2,2]\times[-2,2].             \tag{17}
$$
If either leading parity coordinate is zero the same conclusion holds
without choosing a sign for it. There is no omitted zero-coordinate
direction or unbounded parameter range.

The function $F(v)+|v|^2$ is convex. Iterated interpolation in
the four coordinates of a box $B=\prod_i[l_i,u_i]$ therefore gives
$$
 F(v)\le\max_{w\in\operatorname{Vert}(B)}F(w)
          +\frac14\sum_i(u_i-l_i)^2,
 \qquad v\in B.                                       \tag{18}
$$
The exact quadratic interpolation error before maximizing is
$\sum_i(v_i-l_i)(u_i-v_i)$. This proves the assertion on every
point and face of the box.

A finite sufficient certificate consists of binary subdivisions of
(17), bisecting a longest side with the smallest coordinate index
breaking ties, together with a rigorous rational upper bound for each
vertex value. Distinct leaves must be prefix-free and must have exact
weights $2^{-\text{path length}}$ summing to one. Those conditions
prove coverage: every missing node of the finite binary subdivision
would leave positive dyadic weight. Shared boundary faces are covered
by the closed boxes. For every leaf require
$$
 \max_{w\in\operatorname{Vert}(B)}U(w)
 +\frac14\sum_i(u_i-l_i)^2\le t+\kappa,
 \qquad F(w)\le U(w).                                  \tag{19}
$$
Equations (16)--(19) reduce each target, with all its function
quantifiers, to finite subdivision coverage and finite vertex integrals.

For completeness, the vertex integral has a fully bounded elementary
verification. For the cubic $P(s)=2\psi(s)^TLv$,
$$
 \max(t|s|,|P(s)|)=\tfrac12(|P(s)+ts|+|P(s)-ts|).
$$
On an interval $[l,u]$, an exact polynomial Taylor expansion at
the midpoint $m$, with $r=(u-l)/2$, gives
$$
 |P(s)-P(m)|\le\sum_{j=1}^{\deg P}|P^{(j)}(m)|r^j/j!.
$$
This proves a sign when the resulting interval avoids zero. On such
cells, integrate $P$ with its proved sign. The needed Gaussian
monomial moments are
$$
 I_0(l,u)=\Phi(u)-\Phi(l),\quad I_1(l,u)=\phi(l)-\phi(u),
 \quad I_j(l,u)=l^{j-1}\phi(l)-u^{j-1}\phi(u)+(j-1)I_{j-2}(l,u).
$$
On every remaining cell, use lower bound zero and upper bound
$\sup_{[l,u]}|P|\,I_0(l,u)$. Subdivision improves this bound;
stopping at any finite depth preserves validity because the full
remainder is included. Starting with the unit intervals in $[-12,12]$
gives a finite covering without any root assumptions.

For the two-sided tails at $R=12$, put
$$
 J_0=\operatorname{erfc}(R/\sqrt2),\quad J_1=\nu e^{-R^2/2},
 \quad J_j=\nu R^{j-1}e^{-R^2/2}+(j-1)J_{j-2}.
$$
If $P(s)=\sum c_js^j$, then the complete omitted Gaussian integral
is at most $\sum_j|c_j|J_j$. Outward arithmetic for every coefficient,
range, moment, sum, and final rational rounding proves (19). Floating
priorities for choosing the next cell have no bearing on the inequality.

The complete rational data are in
[certificates_02/fibers.json](certificates_02/fibers.json). The domain and
both matrices are specified in that file. For each target, `leaf_paths`
lists the binary paths in the subdivision just defined, and every row
of `vertices` is $(v_0,v_1,v_2,v_3,U(v))$, with all five entries
given as exact rational strings. The vertex list is exactly the union
of the vertices of the listed leaves. The lists have the following sizes:
target A has 926 leaves and 2,527 vertices, with maximum path length 22;
target B has 1,040 leaves and 3,111 vertices, with maximum path length 32.
For each target the paths are prefix-free and
$\sum_{\text{leaves}}2^{-\text{path length}}=1$.

The procedure above encloses every vertex integral at 192-bit precision.
It starts from the 24 unit intervals in $[-12,12]$, and subdivides
until the sum of unresolved upper contributions is at most $2^{-52}$,
or until 20,000 subdivisions have been performed. At either stopping
condition the complete unresolved contribution and both infinite tails
remain in the upper bound. The final value is rounded upward to a
multiple of $2^{-60}$. Rational assembly of these bounds gives
$$
\begin{split}
 \min_{B\text{ in A}}\left\{\tfrac32+\tfrac{26}{125}
 -\max_{v\in\operatorname{Vert}(B)}U(v)
 -\tfrac14\sum_i(u_i-l_i)^2\right\}
 &=\frac{20720380967041001}{144115188075855872000}>0,\\
 \min_{B\text{ in B}}\left\{\tfrac{11}{5}+\tfrac1{6000}
 -\max_{v\in\operatorname{Vert}(B)}U(v)
 -\tfrac14\sum_i(u_i-l_i)^2\right\}
 &=\frac{311040589043071}{27021597764222976000}>0.
\end{split}
$$
Thus every leaf satisfies (19), proving $(T_A)$ and $(T_B)$.
Section 7 specifies the arithmetic implementation of this finite part.

## 5. Coverage of the remaining first-moment intervals

Using $(T_A)$, (15) gives
$$
 L_a/\nu^2\le ax^2+\tfrac32(1-x)+\tfrac{26}{125}.
$$
This is convex in $x$. At $x=11/50$ and $x=2/3$, its gaps
below $C$ are $37/20000$ and $509/18000$. Thus it covers
the entire closed interval $[11/50,2/3]$.

For $x\ge2/3$,
$\mathbb E S^3h_0\le\mathbb E|S|^3=2\nu$. Therefore
$$
 m_3=\mathbb E h_0\psi_3(S)
 \le\nu(2-3x)/\sqrt6\le0,
 \quad
 \|P_3h_0\|_2^2\ge m_3^2\ge\nu^2(3x-2)^2/6.          \tag{20}
$$
The sign condition is necessary for the last squared inequality and
is why (20) is used only on this interval. Using $(T_B)$, (15)
becomes, with $y=1-x\in[0,1/3]$,
$$
 L_a/\nu^2\le
 ax^2-(3x-2)^2/6+\tfrac{11}{5}(1-x)+\tfrac1{6000}
 =C-y/20+y^2/8\le C-y/120\le C.                       \tag{21}
$$
Equations (11), (15), and (21) cover all of $[0,1]$, including
both junctions and the endpoints. They first prove the upper bound
for $L_a$. Apply that same global result to $(-f,g)$ to prove
the lower bound and hence (1). The value of $x$ may change when a
sign is negated; global coverage makes this immaterial.

## 6. A finite-matrix implication with all errors

For each finite $n$, define the real finite-rank kernel
$$
 T_n(X,Y)=a\sum_{i=1}^n X_iY_i
           -\sum_{|\alpha|=3}\psi_\alpha(X)\psi_\alpha(Y).
                                                               \tag{22}
$$
It represents $aP_1-P_3$, has $L^2$ operator norm $a$, and
lies in $L^2(\gamma_n\otimes\gamma_n)$. Its squared kernel norm
is $a^2n+\binom{n+2}{3}$. It has negative cubic eigenvalues.
The definition of the real Grothendieck constant permits general real
matrices and independent row and column signs; it does not require a
positive semidefinite matrix or equal signs on the two sides.

Write $D=\nu^2C=2917/(1000\pi)$. Equation (1) gives
$$
 \left|\iint f(X)T_n(X,Y)g(Y)\,d\gamma_n(X)d\gamma_n(Y)\right|
 \le D                                                     \tag{23}
$$
for every measurable sign pair, uniformly in $n$.

Use the unit-vector map $U(X)=X/|X|$, defined arbitrarily at zero.
Let $w_j=\sum_i\|P_jU_i\|_2^2$. Rotational symmetry gives
$w_1=(\mathbb E|X|)^2/n$, and orthogonality gives
$w_1+w_3\le1$. Put $Z=|X|/\sqrt n$. Since
$\mathbb EZ^2=1$, $\mathbb E(Z^2-1)^2=2/n$, and
$(Z-1)^2\le(Z^2-1)^2$,
$2-2\mathbb EZ\le2/n$. Hence
$$
 w_1\ge(1-1/n)^2,\qquad
 \iint T_n(X,Y)\langle U(X),U(Y)\rangle\,d\gamma_n^2
 =aw_1-w_3\ge(a+1)(1-1/n)^2-1.                        \tag{24}
$$
This tends to $a$. It is enough to retain the full orthogonal
complement in the inequality $w_3\le1-w_1$; no Hermite tail
limit is assumed.

For integers $R\ge2$, $q\ge0$, set $\eta=2^{-q}$, partition
$[-R,R)^n$ into cubes of side $\eta$, and add its complementary
cell. For each cell $B$ set
$\mu_\alpha(B)=\int_B\psi_\alpha\,d\gamma_n$. Define the finite
matrix
$$
 A_{B,B'}=a\sum_i\mu_{e_i}(B)\mu_{e_i}(B')
              -\sum_{|\alpha|=3}\mu_\alpha(B)\mu_\alpha(B').
                                                               \tag{25}
$$
The matrix has order $(2R2^q)^n+1$. Cube moments factor into
$$
 M_0(l,u)=\Phi(u)-\Phi(l),\qquad
 M_j(l,u)=\frac{\phi(l)H_{j-1}(l)-\phi(u)H_{j-1}(u)}{\sqrt{j!}},
 \quad j=1,2,3.                                       \tag{26}
$$
These identities follow from $(\phi H_{j-1})'=-\phi H_j$.
Complementary moments are the full moment minus the finite sum of
cube moments. For the odd total degrees in (25) they are zero by
global inversion symmetry, up to Gaussian-null boundary faces.
Equations (25)--(26) specify exact,
computable real entries; there is no matrix-entry rounding error.

Any row signs and independent column signs define two measurable step
functions. Expanding their finite bilinear form in (25) gives exactly
the left side of (23). Thus $S_{\mathbb R}(A)\le D$ for all
assignments simultaneously. No maximizer is passed through a limit.

For vector values, assign to each cube the unit direction of its
center and to the complementary cell the vector $e_1$. Every
center is nonzero, since each coordinate is a half-integer multiple
of $\eta$. Call this unit step function $G$. Its matrix value
equals its kernel vector value by the same finite expansion. Inside
the box and outside the unit ball,
$$
 |U(X)-G(X)|\le\frac{2|X-\operatorname{center}(B)|}{|X|}
 \le\sqrt n\eta.
$$
On the exceptional events the distance is at most two. Markov's
inequality and the Gaussian Laplace integral give
$$
 \Pr(|X|<1)\le e\,3^{-n/2}<3\,3^{-n/2}.
$$
Integration of $s\phi(s)$ gives
$\Pr(|S|\ge R)\le(\nu/R)e^{-R^2/2}\le e^{-R^2/2}$.
A union bound on coordinates consequently proves
$$
 \|U-G\|_{L^2(\gamma_n;\mathbb R^n)}^2
 \le n2^{-2q}+12\,3^{-n/2}+4ne^{-R^2/2}=:\epsilon^2.     \tag{27}
$$
All Gaussian tails are present in (27). Overlap of exceptional events
only enlarges this upper bound. The vector-valued extension of
$aP_1-P_3$ has the same norm $a$, by applying its scalar norm
bound coordinatewise and summing squares. Therefore
$$
 V_{\mathbb R}(A)\ge(a+1)(1-1/n)^2-1-2a\epsilon.         \tag{28}
$$
Here the error estimate follows from
$\langle U,TU\rangle-\langle G,TG\rangle
=\langle U-G,TU\rangle+\langle G,T(U-G)\rangle$, with both
$L^2$ vector norms equal to one.

Taking $n\ge2$ and $R=q=n$ produces a sequence of specified finite matrices
with $\epsilon\to0$. Polynomial factors times $4^{-n}$,
$3^{-n/2}$, and $e^{-n^2/2}$ tend to zero by the ratio test.
For sufficiently large $n$ the vector value is positive, so the
matrix is nonzero. Every nonzero finite matrix has positive sign norm:
the average squared bilinear form over independent uniform signs is
$\sum_{B,B'}A_{B,B'}^2>0$. Equations (23) and (28) then imply
$$
 K_G^{\mathbb R}\ge\frac aD=\frac{1625\pi}{2917}.        \tag{29}
$$

There is also a single fully specified finite example. Choose
$$
 n=100000,\qquad R=20,\qquad q=40.
$$
Since $e>2$, $3^{-50000}<2^{-100}$, and $e<3$, (27) gives
$$
 \epsilon^2<100000\,2^{-80}+12\,2^{-100}+400000\,2^{-200}
 <10^{-18}.
$$
Its vector value is therefore greater than
$129995799761/80000000000$. With the same rational lower bound
for $\pi$, (23) proves the strict finite-matrix ratio bound
$$
 \frac{V_{\mathbb R}(A)}{S_{\mathbb R}(A)}
 >\frac{408393849526687368661439473}{233360000000000000000000000}
 =\frac74+
   \frac{13849526687368661439473}{233360000000000000000000000}
 >\frac74.                                             \tag{30}
$$
The order of this matrix is $(40\cdot2^{40})^{100000}+1$.
The statement specifies its entries and vectors exactly; it does not
require storing an enumeration of its entries.

A rational matrix can also be obtained with a specified total entry
error. Write $N=(40\cdot2^{40})^{100000}+1$,
$W=129995799761/80000000000$, and
$\overline D=2917000000000000/3141592653589793>D$.
For every entry of (25), compute a rational interval of radius less than
$10^{-6}/N^2$ and take its midpoint, obtaining a rational matrix
$\widetilde A$. These intervals can be computed to arbitrary accuracy
from (26): the exponential and error-function series converge at every
finite argument, and beyond their decreasing-term threshold their
alternating tails are bounded by the next term. Rational square-root
bisection and arbitrarily long versions of the $\pi$ series below
complete these computations. Every requested accuracy is positive and
there are finitely many entries, so the construction terminates.

Then $\sum_{i,j}|\widetilde A_{ij}-A_{ij}|<10^{-6}$.
For every sign assignment and every unit-vector assignment, changing
the matrix changes its bilinear value by at most this sum. Therefore
$$
 \frac{V_{\mathbb R}(\widetilde A)}{S_{\mathbb R}(\widetilde A)}
 >\frac{W-10^{-6}}{\overline D+10^{-6}}
 =\frac{136131199399758360492666491}{77786750442470762394480000}
 =\frac74+
   \frac{4386125434526302326491}{77786750442470762394480000}
 >\frac74.
$$
The positive vector lower bound also proves that this rational matrix
is nonzero, so its sign norm is positive.

The rational bound for $\pi$ used above follows from
$\pi=16\arctan(1/5)-4\arctan(1/239)$ and alternating series
bounds. The tangent subtraction identity verifies this formula in the
first quadrant. Thirty-two terms for $\arctan(1/5)$ and eight
terms for $\arctan(1/239)$, retaining the next term as the error,
give
$3141592653589793/10^{15}<\pi<3141592653589794/10^{15}$.

## 7. Finite arithmetic and implementation dependencies

The finite vertex procedure is implemented in
[check_fibers.py](check_fibers.py), in `make_verifier`. Its inputs are
the displayed rational matrices, a rational vertex, and the positive
numbers $t$ and $\nu$. The coefficient construction is exactly
$$
\begin{split}
 d_0&=2\sqrt p\,v_0,&
 d_2&=2\{\sqrt{r^2/p}\,v_0+\sqrt{q-r^2/p}\,v_2\},\\
 d_1&=2\sqrt{p_o}\,v_1,&
 d_3&=2\{\sqrt{r_o^2/p_o}\,v_1+
                   \sqrt{q_o-r_o^2/p_o}\,v_3\}.
\end{split}
$$
Here $p_o=139/118$, $q_o=63/59$,
$r_o^2=1323/6962$, and
$P(s)=d_0+d_1s+d_2(s^2-1)/\sqrt2+d_3(s^3-3s)/\sqrt6$.
Applying the absolute-integral procedure to $P+ts$ and $P-ts$
computes an upper bound for $F(v)$. The polynomial coefficients,
Taylor ranges, Gaussian moments, all accumulated remainders, and the
division by $\nu$ are enclosed outward.

The implementation uses CPython integer arithmetic and its standard
`fractions.Fraction` type for exact rational data. Nonrational
calculations use `python-flint`, with FLINT/Arb and its GMP/MPFR
arithmetic dependencies. Required real-ball operations are addition,
subtraction, multiplication, division by a positive enclosure, integer
powers, square root, exponential, error function and complementary
error function, $\pi$, outward bounds, and integer ceiling.
The additional noise proof also uses logarithm, sine, cosine, and inverse
sine. [FLINT's Arb specification](https://flintlib.org/doc/arb.html)
states the enclosure contract; the
[python-flint implementation](https://raw.githubusercontent.com/flintlib/python-flint/0.8.0/src/flint/types/arb.pyx)
supplies the corresponding outward-bound and comparison operations.

An interval comparison succeeds only when it proves the comparison for
all enclosed values. Upper and absolute upper bounds are extracted
outward. A final ceiling must yield one exact integer; the resulting
rational bound is compared back to the preceding upper bound. Thus a
precision failure cannot silently replace an uncertain inequality by a
positive decision. A floating-point priority changes only the order of
subdivision; every unresolved cell is included at termination.

The exact coverage and final rational sums are also implemented in
[check_certificate_arithmetic_02.py](check_certificate_arithmetic_02.py).
That calculation checks the supplied finite bounds and their assembly;
the vertex evaluator establishes the bounds themselves.
[check_final_exact_02.py](check_final_exact_02.py), using the elementary
polynomial operations in [check_exact.py](check_exact.py), checks the
circle identities, matrix determinants, rational $\pi$ bounds, and
the finite-matrix margins. The mathematical certificates consist of
the rational data and these explicitly bounded evaluations, together
with the analytic proofs above.
