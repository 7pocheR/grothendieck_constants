# Quintic Gaussian operator and finite real matrices

Put
$$
 a=\frac{33}{20},\qquad \lambda=\frac3{20},\qquad
 \nu=\sqrt{\frac2\pi},\qquad C=\frac{3681}{2500}.
$$
This document specifies a finite certificate for the bound
$K_{\mathbb R}\ge1375\pi/2454$, and proves its analytic implications.
The finite hypothesis is the collection of six-dimensional vertex inequalities
in Section 6. The agreement integrals and scalar interval inequalities are
given as rational certificates in `SCALAR_MATHEMATICAL_DATA_01.json`.

## 1. Statement and normalization

For a real finite matrix $A=(A_{ij})$, write
$$
 B(A)=\max_{\epsilon_i,\delta_j\in\{-1,1\}}
             \left|\sum_{ij}A_{ij}\epsilon_i\delta_j\right|.
$$
The real Grothendieck constant is the supremum, over all nonzero real finite
matrices and all finite-dimensional real unit vectors, of
$\big|\sum A_{ij}\langle u_i,v_j\rangle\big|/B(A)$.
The row and column signs are chosen independently.

Let $\gamma_n$ be standard Gaussian measure, with covariance $I_n$.
Let $P_j$ be the orthogonal projection onto homogeneous Hermite degree
$j$ in $L^2(\gamma_n)$. The desired universal inequality is
$$
 \left|a\langle P_1f,P_1g\rangle-\langle P_3f,P_3g\rangle
                 -\lambda\langle P_5f,P_5g\rangle\right|
 \le \nu^2 C                                                   \tag{1}
$$
for every finite $n\ge1$ and every pair of measurable real sign functions.

**Conditional theorem.** If the explicit finite vertex inequalities in
Section 6 hold for the three matrices in `targets_01.json`, then (1) holds.
Consequently
$$
 K_{\mathbb R}\ge\frac{a}{\nu^2C}=\frac{1375\pi}{2454}.
                                                               \tag{2}
$$
Under the same finite hypothesis, Section 9 specifies one actual finite
matrix whose ratio is strictly greater than $44/25$, hence than $7/4$.

The Hermite normalization is
$\psi_j(s)=H_j(s)/\sqrt{j!}$, with
$H_0=1,H_1=s,H_{j+1}=sH_j-jH_{j-1}$. Products
$\psi_\alpha(X)=\prod_i\psi_{\alpha_i}(X_i)$ form an orthonormal
basis. Their generating function
$e^{t\cdot X-|t|^2/2}$ proves both orthogonality and the identity
$$
 \mathbb E\psi_\alpha(X)\psi_\beta(Y)
      =\mathbf1_{\alpha=\beta}\,r^{|\alpha|},                 \tag{3}
$$
when $X,Y$ are jointly standard Gaussian with cross covariance $rI_n$,
including $-1\le r\le1$. For completeness, a function orthogonal to
all polynomials has identically zero entire Gaussian Laplace transform:
Cauchy--Schwarz gives its convergence on every compact complex set, and all
derivatives at zero vanish. Its Fourier transform therefore vanishes, so
the function is zero. This also justifies the complete Parseval expansions
used below; no Hermite series is numerically truncated.

## 2. Original midpoints and the parameter $x$

For the original signs define
$$
 h=\frac{f+g}{2},\qquad k=\frac{f-g}{2},\qquad
 x=\frac{\|P_1h\|_2}{\nu}.
$$
Pointwise, $(h,k)$ is one of $(1,0),(-1,0),(0,1),(0,-1)$.
In particular $|h|=h^2$, $hk=0$, and $h^2+k^2=1$.
No oddification of the original functions is made.

Choose a unit Gaussian coordinate $S$ in the direction of $P_1h$
when it is nonzero, and any coordinate otherwise. Then
$\mathbb E Sh=\nu x$, and $|h|\le1$ gives $0\le x\le1$.
All occurrences of $x$ below refer to this same quantity.
Put $H_j=\|P_jh\|_2^2$ and $K_j=\|P_jk\|_2^2$ in formulas
where $H_j$ denotes an energy rather than a one-variable polynomial.
The unnormalized bilinear expression in (1) is exactly
$$
 L=a(H_1-K_1)-H_3+K_3-\lambda H_5+\lambda K_5.                \tag{4}
$$
This follows by expanding both midpoints and using symmetry of each
orthogonal projection. Proving the upper bound for every pair proves the
absolute bound as well: replace $g$ by $-g$ to negate $L$.
That replacement is allowed to change $x$.

## 3. A circle estimate, with the quintic term retained

We first prove, for any real trigonometric polynomial $p=p_1+p_3$
with frequencies one and three and normalized circle measure,
$$
 \left(\int|p|\right)^2
       \le\frac9{\pi^2}\int p_1^2+\frac8{\pi^2}\int p_3^2.  \tag{5}
$$
The case $p=0$ is immediate. Otherwise, $p(\theta+\pi)=-p(\theta)$,
and there are at most six zeros on the full circle, counted with
multiplicity: multiplying its Laurent polynomial by $z^3$ gives a
nonzero polynomial of degree at most six. Thus $\sigma=\operatorname{sgn}p$
has either one or three sign changes in a half circle; even-multiplicity
zeros do not contribute. All zero sets here have measure zero.

In the one-change case its frequency energies are
$\|\sigma_j\|_2^2=8/(\pi^2j^2)$, for $j=1,3$.
In the three-change case rotate the changes to $-\alpha,0,\beta$,
where $\alpha,\beta>0$ and $\alpha+\beta<\pi$. Set
$$
 u=\cos\frac{\alpha+\beta}{2},\qquad
 v=\cos\frac{\alpha-\beta}{2},\qquad 0\le u\le v\le1.
$$
Integration of the three jumps and their opposite half-circle jumps gives
$$
 \|\sigma_j\|_2^2=\frac8{\pi^2j^2}Z_j,\quad
 Z_1=1+4u^2-4uv,\quad
 Z_3=1+4T_3(u)^2-4T_3(u)T_3(v),
$$
where $T_3(u)=4u^3-3u$. Equivalently the complex jump sum is
$e^{ij\alpha}-1+e^{-ij\beta}$, whose squared absolute value is
$Z_j$. Direct polynomial identities give
$$
 9-8Z_1-Z_3=4u(v-u)B(u,v),
$$
$$
 B=8+(4u^2-3)\{4(u^2+uv+v^2)-3\}.
$$
If $u\ge\sqrt3/2$, both factors in braces and before it are
nonnegative. If $u\le\sqrt3/2$, use instead the identical expression
$$
 B=(2u-1)^2(4u^2+8u+5)
          +4(1-v)(1+v+u)(3-4u^2)\ge0.
$$
Thus $8Z_1+Z_3\le9$, with boundary and multiple-zero cases included.
Weighted Cauchy--Schwarz applied to $\int p\sigma=\int|p|$ now proves
(5), since
$\|\sigma_1\|_2^2/(9/\pi^2)+\|\sigma_3\|_2^2/(8/\pi^2)\le1$.

For Gaussian polynomials $F=F_1+F_3$ in Hermite degrees one and three,
apply (5) to $F(X\cos\theta+Y\sin\theta)$, where $X,Y$ are
independent standard Gaussian vectors. Formula (3) gives covariance
$\|F_1\|_2^2\cos t+\|F_3\|_2^2\cos^3t$. Since
$\cos^3t=(3\cos t+\cos3t)/4$, the expected two circle energies are
$\|F_1\|_2^2+3\|F_3\|_2^2/4$ and $\|F_3\|_2^2/4$.
Jensen's inequality therefore gives
$$
 (\mathbb E|F|)^2\le\frac9{\pi^2}\|F_1\|_2^2
                       +\frac{35}{4\pi^2}\|F_3\|_2^2.
$$
For any $|w|\le1$, insert
$F=(\pi^2/9)P_1w+(4\pi^2/35)P_3w$.
Writing $D=(\pi^2/9)\|P_1w\|_2^2+(4\pi^2/35)\|P_3w\|_2^2$,
we obtain $D=\langle w,F\rangle\le\mathbb E|F|\le\sqrt D$.
Consequently
$$
 \|P_3w\|_2^2+\frac{35}{36}\|P_1w\|_2^2
                      \le\frac{35}{4\pi^2}.                 \tag{6}
$$

Apply (6) to $f,g$ and average. In addition, full Parseval gives
$H_1+K_1+H_3+K_3+H_5+K_5\le1$. Since $0<\lambda<1$,
$$
 K_3+\lambda K_5
 \le(1-\lambda)\left[\frac{35}{4\pi^2}-H_3
                   -\frac{35}{36}(H_1+K_1)\right]
       +\lambda(1-H_1-K_1-H_3-H_5).
$$
Substitute this inequality in (4), retaining all remaining negative terms:
$$
 \frac L{\nu^2}\le
 \frac{3\pi}{40}+\frac{119}{32\pi}+\frac{97}{144}x^2
       -\frac{1891}{720}\frac{K_1}{\nu^2}
       -2\frac{H_3}{\nu^2}-\frac3{10}\frac{H_5}{\nu^2}.       \tag{7}
$$
Dropping those three nonpositive terms proves the asserted small-$x$
estimate. It is increasing for $x\ge0$. The rational bounds
$$
 \frac{3141592653589793}{10^{15}}<\pi
       <\frac{1570796326794897}{5\cdot10^{14}}                \tag{8}
$$
give a strictly positive gap greater than $2545/10^7$ below $C$
at $x=7/25$. Thus (1)'s upper bound holds on $[0,7/25]$.
The certificates for (8) use the alternating series in
$\pi=16\arctan(1/5)-4\arctan(1/239)$, with 32 and 8 terms.
The identity follows from the tangent addition formula and its branch
$0<4\arctan(1/5)-\arctan(1/239)<\pi/2$.

## 4. Complete transverse reduction

Write the remaining Gaussian coordinates as $T$, independent of $S$,
and define
$k_j(T)=\mathbb E_S k(S,T)\psi_j(S)$, for $0\le j\le5$.
In dimension one, $T$ is a singleton. Put
$$
 D=\operatorname{diag}(1,\lambda,0),\quad
 E=\operatorname{diag}(-a,1,\lambda),\quad
 J=\operatorname{diag}(\lambda,0,0).
$$
Take any real symmetric $M_e,M_o$ each dominating $D$ and $E$,
and set
$$
 Q(k)=(k_0,k_2,k_4)M_e(k_0,k_2,k_4)^T
              +(k_1,k_3,k_5)M_o(k_1,k_3,k_5)^T.
$$
These matrices also dominate $J$ and zero, because $D\succeq J,0$.
Expand each $k_j$ in all Hermite degrees of $T$. At transverse
degree $d$, the desired coefficient of $|k_{j,\beta}|^2$ is
$-a\mathbf1_{j+d=1}+\mathbf1_{j+d=3}+\lambda\mathbf1_{j+d=5}$.
The nonzero even-coordinate blocks are $E,D,J$ at $d=1,3,5$,
respectively. The nonzero odd-coordinate blocks are $E,D,J$ at
$d=0,2,4$, respectively. Every other block is zero.
Coefficients with longitudinal degree above five are also zero.
Applying matrix domination separately for every transverse multi-index gives
$$
 K_3+\lambda K_5-aK_1\le\mathbb E_T Q(k(T)).                 \tag{9}
$$
This statement includes arbitrarily high transverse degrees. Its right
side is a convergent nonnegative quadratic sum by Parseval and boundedness
of the finite matrices; the left involves bounded orthogonal projections.
Taking finite partial sums and their limits proves (9) for all measurable
original signs, without a degree restriction.

For a one-variable ternary $u$, put
$m_j=\mathbb E u(S)\psi_j(S)$ and
$\beta=\mathbb E|S|u(S)^2$. Suppose
$$
 Q(m)/\nu^2\le t\beta/\nu+\kappa\quad\hbox{for every ternary }u,
 \qquad t>0.                                                \tag{10}
$$
For almost every $T$, the actual $k(S,T)$ is ternary. Moreover
$$
 \mathbb E_T\beta(T)=\mathbb E|S|(1-|h|)
       \le\nu-|\mathbb E Sh|=\nu(1-x).
$$
Fubini is applicable because all moments are integrable. Combining (4),
(9), and (10) yields
$$
 L/\nu^2\le ax^2+t(1-x)+\kappa-H_3/\nu^2-\lambda H_5/\nu^2.\tag{11}
$$
The sign of $t$ is essential to this use of the budget.

## 5. Duality and a complete six-dimensional domain

Interleave the two parity matrices to write $Q=LL^T$, using the
lower triangular factor with positive diagonal in each parity block.
Write $c_j(s)=\sum_iL_{ij}\psi_i(s)$ and define
$$
 F(v)=\nu^{-1}\mathbb E\max(t|S|,2|\sum_jv_jc_j(S)|)-\|v\|^2.
                                                               \tag{12}
$$
Completing the square gives
$Q(m)/\nu^2=\sup_v[2\langle L^Tm/\nu,v\rangle-\|v\|^2]$.
Taking the supremum over ternary $u$ and exchanging two suprema, then
choosing $u$ pointwise, gives exactly
$$
 \sup_u\{Q(m)/\nu^2-t\beta/\nu\}=\sup_v F(v)-t.             \tag{13}
$$
The maximizing pointwise choice is
$u=\operatorname{sgn}(\sum v_jc_j)$ when
$2|\sum v_jc_j|>t|S|$, and zero otherwise. A fixed tie rule
is measurable. All relevant polynomial integrals are finite. No minimum
and maximum are interchanged.

The function $F$ is continuous and
$F(v)\le t+2\|L\|\|v\|/\nu-\|v\|^2$, so it has a global
maximizer. Select a pointwise maximizing $u$ at such a maximizer
$v_*$. For this fixed $u$, the completed-square expression has
the unique maximizer $L^Tm/\nu$. If it differed from $v_*$, its
strictly larger value would contradict global maximality in (13). Therefore
$$
 |(v_*)_j|\le b_j:=\nu^{-1}\mathbb E|c_j(S)|.                \tag{14}
$$
This argument covers every global maximizer and all boundary cases.

Reflection $S\mapsto-S$ flips the odd part of $\sum v_jc_j$
without changing (12). Negating all coordinates also preserves (12).
Together these give independent sign reversal of each parity block.
Hence a global maximizer can be chosen with $v_0,v_1\ge0$.
If $B_j\ge b_j$, the full required domain is consequently
$$
 [0,B_0]\times[0,B_1]\times\prod_{j=2}^5[-B_j,B_j].         \tag{15}
$$
Zeros in either first coordinate cause no exception. This is a symmetry
of the integral, not a restriction on the original signs.

## 6. The finite certificate

The six rational matrices are exactly those in `targets_01.json`, with
the even matrix first. Every principal minor of each matrix and of its
differences from $D,E$ is checked in exact rational arithmetic.
The matrices themselves are positive definite. Checking all principal
minors suffices for positive semidefiniteness of a real symmetric matrix:
add $\epsilon I$, expand its principal minors in $\epsilon>0$,
apply the positive-definite criterion, and take $\epsilon\downarrow0$.
Exact rational $N\Delta N^T$ factorizations with positive $\Delta$
give the required Cholesky factors $N\sqrt\Delta$.

The three choices $(t,\kappa; B_0,\ldots,B_5)$ are
$$
\begin{aligned}
 A:&\quad (3/2,6/25;\ 11/8,19/16,13/16,9/8,13/16,5/16),\\
 B:&\quad (2,1/50;\ 3/2,17/16,15/16,7/8,13/16,3/8),\\
 C:&\quad (4371/2000,1/4000;\ 25/16,1,7/8,13/16,3/8,3/8).
\end{aligned}
$$
Their required upper bounds for $F$ are respectively
$87/50,101/50,8743/4000$. All 18 inequalities (14) with these
rational $B_j$ have direct one-dimensional integral certificates.

For each rectangle $R=\prod_j[l_j,u_j]$, convexity of
$F(v)+\|v\|^2$ gives
$$
 \sup_R F\le\max_{z\in\operatorname{vertices}(R)}F(z)
                            +\tfrac14\sum_j(u_j-l_j)^2.     \tag{16}
$$
Indeed the product distribution on vertices with mean $v$ has variance
$\sum_j(v_j-l_j)(u_j-v_j)$, at most the displayed correction.

The finite hypothesis of the theorem is this: the terminal rectangles
form a complete binary subdivision of each domain (15), and every one
of them satisfies
$$
 \max_{z\in\operatorname{vertices}(R)}U(z)
       +\tfrac14\sum_j(u_j-l_j)^2\le t+\kappa,               \tag{17}
$$
where each $U(z)$ is a verified upper bound for the exact function
(12). Every terminal rectangle has 64 vertices. Binary paths specify
successive bisections of the longest side, with the smallest coordinate
index breaking ties. In each initial path of length four, distinct
prefix-free terminal paths with sum $\sum_p2^{-(|p|-4)}=1$
give complete coverage: expand the finite collection to the maximum
depth and count every binary descendant. This proves coverage of all
closed boxes and their shared boundaries, not merely almost everywhere
in the six-dimensional parameter space.

For evaluating a vertex, use
$$
 F(v)=\frac{\mathbb E|2\sum_jv_jc_j(S)+tS|
                 +\mathbb E|2\sum_jv_jc_j(S)-tS|}{2\nu}-\|v\|^2.
$$
Each integral is of an absolute polynomial of degree at most five.
`REPLAY_MATHEMATICS_01.md` specifies outward interval Horner range bounds,
signed Gaussian antiderivatives, exact dyadic rounding, and complete
two-sided tails. An unresolved sign interval contributes its full
polynomial-range upper bound times its Gaussian mass. All such contributions
remain in the result. At radius $R$ the exact absolute monomial tails are
$$
 T_0=\operatorname{erfc}(R/\sqrt2),\quad T_1=\nu e^{-R^2/2},\quad
 T_j=\nu R^{j-1}e^{-R^2/2}+(j-1)T_{j-2}.
$$
The estimate $\sum_{j=0}^5|p_j|T_j$ includes every sign change in both
infinite tails. Thus (17) is a finite collection of explicitly computable
inequalities, with no residual continuum or tail hypothesis.

## 7. Agreement supports and their complete integrals

For any real $\ell,z_3,z_5$, completing two squares in the
longitudinal agreement coefficients gives
$$
 -H_3/\nu^2-\lambda H_5/\nu^2
 \le\ell x+z_3^2+z_5^2/\lambda
        +\nu^{-1}\mathbb E|\ell S+2z_3\psi_3(S)+2z_5\psi_5(S)|.
                                                               \tag{18}
$$
For example, $-y^2\le z^2-2zy$ is the identity
$(y-z)^2\ge0$, and $-\lambda y^2\le z^2/\lambda-2zy$
is analogous. The complete projection energy dominates the square of
its longitudinal coefficient. Add
$\ell(x-\mathbb E Sh/\nu)=0$, and then use $|h|\le1$.
This proves (18) globally for every parameter, not only near a tangency
point or for interval or halfspace signs.

For rational $r\ge0,u,v$, take
$$
 \ell=-u(r-3)-v(r^2-10r+15),\qquad
 z_3=u\sqrt6/2,\quad z_5=v\sqrt{120}/2.
$$
The polynomial under the absolute value is $S Q(S^2)$, where
$$
 Q(y)=A_0+B_0y+vy^2=(y-r)\{vy+u+v(r-10)\},
 \quad A_0=\ell-3u+15v,\quad B_0=u-10v.
$$
The measure $|S|d\gamma_1(S)/\nu$, pushed forward by $S\mapsto S^2$,
has density $e^{-y/2}/2$ on $[0,\infty)$. Thus the constant in
(18) is bounded by
$$
 c\ge\frac32u^2+\frac{30v^2}{\lambda}
                +\frac12\int_0^\infty |Q(y)|e^{-y/2}\,dy.   \tag{19}
$$
The only possible positive roots are $r$ and $10-r-u/v$ if
$v\ne0$; if $v=0,u\ne0$, the sole possible root is $r$.
The zero polynomial is treated separately. Repeated roots, zero roots,
and negative second roots cause no missing interval.

A signed primitive is
$$
 J(y)=-e^{-y/2}[A_0+B_0(y+2)+v(y^2+4y+8)],\qquad J(\infty)=0.
$$
Differentiation gives $J'=Qe^{-y/2}/2$. Split at all distinct positive
roots and determine the sign at any rational interior point. The full
integral is a rational linear combination of $e^{-q}$ with rational
$q\ge0$. This includes the final interval to infinity exactly.

For each of the 1,007 supports in `SCALAR_MATHEMATICAL_DATA_01.json`, (19)
is verified using rational bounds alone. To enclose $e^{-q}$, divide
$q$ by a power of two until $z\le1/2$. The alternating sums
$S_{39}(z)\le e^{-z}\le S_{40}(z)$ are exact rational bounds;
successive squaring, with outward rational dyadic rounding after each
step, recovers an enclosure at the original argument. This proof applies
also to very large arguments and never omits a small exponential.
The zero support gives zero. Every other supplied constant has a strictly
positive rational margin in (19).

## 8. Exact assembly on the entire remaining interval

Let $A_*(x)=\min_i(\ell_i x+c_i)$ be the minimum of the 1,007
proved agreement lines. From (11), (18), and the three bounds (10),
$$
 L/\nu^2\le ax^2+A_*(x)+F_*(x),                            \tag{20}
$$
where
$$
 F_*(x)=\min\{87/50-3x/2,\ 101/50-2x,\
                         8743/4000-4371x/2000\}.
$$
Each fiber inequality holds for the same actual $h,k,x$, and each
agreement line bounds the same negative agreement energy. Therefore
the two finite minima in (20) may be taken independently, without any
common optimizer or minimax assertion.

On $[7/25,1]$, the three fiber intervals end at
$7/25,14/25,663/742,1$. The agreement envelope has 197 intervals.
For every such interval the exact active line is compared with all
1,007 lines at both endpoints; affine differences then prove dominance
on the entire closed interval. The union of the two sets of breakpoints
has 200 distinct endpoints, all listed as rationals in the data file.

Between consecutive endpoints, the right side of (20) is a quadratic
with positive leading coefficient $a$, so its maximum is at an
endpoint. Exact rational evaluation of all 200 endpoints gives
$$
 \max=\frac{1618847619965169}{1099511627776000}<C,
$$
with gap
$$
 C-\max=\frac{366503861067}{5497558138880000}>0.              \tag{21}
$$
All tie points and both ends are included. Together with (7), this
proves the universal upper bound on $[0,1]$ under (17), and hence
(1). No noise inequality, Gaussian isoperimetric comparison, or Borell
theorem enters this route.

## 9. An actual finite matrix and all discretization errors

In dimension $n$, define the finite-rank real kernel
$$
 \mathcal K_n(X,Y)=a\sum_{|\alpha|=1}\psi_\alpha(X)\psi_\alpha(Y)
       -\sum_{|\alpha|=3}\psi_\alpha(X)\psi_\alpha(Y)
       -\lambda\sum_{|\alpha|=5}\psi_\alpha(X)\psi_\alpha(Y). \tag{22}
$$
It is the kernel of $T=aP_1-P_3-\lambda P_5$, whose $L^2$
operator norm is $a$. This kernel and the resulting finite matrix
are not assumed positive semidefinite; the negative coefficients are
retained throughout. The real Grothendieck definition permits arbitrary
real matrices.

Let $U(X)=X/\|X\|$, assigning any unit value at zero, and put
$w_j=\sum_i\|P_jU_i\|_2^2$. Full Parseval gives $\sum_jw_j=1$.
Rotational symmetry gives
$\mathbb E U_iX_j=\delta_{ij}\mathbb E\|X\|/n$, so
$w_1=(\mathbb E\|X\|)^2/n$. For $Z=\|X\|/\sqrt n$,
$\mathbb EZ^2=1$ and $\mathbb E(Z^2-1)^2=2/n$.
Because $(Z-1)^2\le(Z^2-1)^2$, we get $\mathbb EZ\ge1-1/n$.
As $\lambda\le1$, the vector value therefore satisfies
$$
 I_n=a w_1-w_3-\lambda w_5
      \ge(a+1)(1-1/n)^2-1.                                \tag{23}
$$
Every unused Hermite degree is covered by $\sum_jw_j=1$.

Partition $[-R,R)^n$, for positive integer $R$, into cubes of side
$\eta=2^{-q}$, and include its complement as one further cell. For
each cell $B$, set $\mu_{B,\alpha}=\int_B\psi_\alpha\,d\gamma_n$.
Define the actual finite matrix
$$
 A_{BB'}=a\sum_{|\alpha|=1}\mu_{B,\alpha}\mu_{B',\alpha}
       -\sum_{|\alpha|=3}\mu_{B,\alpha}\mu_{B',\alpha}
       -\lambda\sum_{|\alpha|=5}\mu_{B,\alpha}\mu_{B',\alpha}. \tag{24}
$$
Its entries are explicit real numbers. For a cube they are products
of one-dimensional moments
$$
 M_0(l,u)=\Phi(u)-\Phi(l),\quad
 M_j(l,u)=\frac{\phi(l)H_{j-1}(l)-\phi(u)H_{j-1}(u)}{\sqrt{j!}}
 \quad (1\le j\le5).
$$
The formula follows by differentiating $\phi H_{j-1}$.
The complement moments for odd total degree are zero by central
symmetry, with the half-open boundary changing only a Gaussian null set.

Independent signs on cells define measurable sign functions $f,g$.
Finite summation and integration identify their matrix value exactly
with $\langle f,Tg\rangle$. Thus (1), already universal before this
partition is selected, proves $B(A)\le\nu^2C$ simultaneously for
all finite row and column sign choices. There is no discretization error
in this scalar upper bound and no exchange of an optimizer with a limit.

Choose the unit vector of a cube to be its center divided by its norm,
and choose $e_1$ for the complement; denote the resulting step field
by $G$. All cube centers are nonzero because each coordinate lies at
a half mesh point in this even-sized grid. For $X$ in a cube centered
at $c$, the triangle inequality gives
$\|X/\|X\|-c/\|c\|\|\le2\|X-c\|/\|X\|$.
For $n>2$, integration by parts in the radial Gaussian integral gives
$\mathbb E\|X\|^{-2}=1/(n-2)$. Furthermore
$\Pr(|S|\ge R)\le e^{-R^2/2}$: substitute $s=R+t$ in
the two-sided tail and discard the factor $e^{-Rt}\le1$.
The union bound and $\|U-G\|\le2$ outside the cube give
$$
 \|U-G\|_{L^2(\gamma_n;\mathbb R^n)}^2
 \le\epsilon^2:=\frac n{n-2}2^{-2q}+4n e^{-R^2/2}.         \tag{25}
$$
The vector-valued operator $T\otimes I_n$ also has norm $a$.
Since $U,G$ have norm one in this Hilbert space,
$$
 \left|\langle G,TG\rangle-\langle U,TU\rangle\right|
       \le2a\epsilon.
$$
Finite expansion of (24) identifies $\langle G,TG\rangle$ with
its finite vector objective. Consequently the single finite matrix has
ratio at least
$$
 \frac{(a+1)(1-1/n)^2-1-2a\epsilon}{\nu^2C}.              \tag{26}
$$
The positive lower bound on its vector objective also proves that the
matrix is nonzero and $B(A)>0$.

Specifically take $n=100000,R=12,q=40$. The matrix order is
$(3\cdot2^{43})^{100000}+1$; (24) specifies every entry without
requiring its enumeration. The elementary inequality
$e>1+1+1/2+1/6=8/3$ gives
$$
 \epsilon^2<\frac{100000}{99998}2^{-80}+400000(3/8)^{72}
                       <\left(\frac1{5\cdot10^{11}}\right)^2.
$$
Using this rational upper bound for $\epsilon$ and the lower bound
for $\pi$ in (8), the exact rational evaluation of (26) is greater
than $1.76020829$, in particular greater than $44/25>7/4$.
The full rational numerator, denominator, and positive gaps are in the
data file. This is an actual finite example, not only a limiting kernel.

Finally take, for example, $R=q=n$ and $n\to\infty$ in (26).
Both terms in (25) tend to zero and (23) tends to $a$. For every
positive error tolerance some finite matrix therefore has ratio at
least $a/(\nu^2C)$ minus that tolerance. The definition as a
supremum over finite matrices proves (2).

## 10. Arithmetic dependencies

The agreement supports, affine envelopes, constants in (7), and the
finite-matrix inequalities are verified by integer and rational arithmetic
and the explicitly bounded alternating series above. Their checker uses
only the Python standard library. Matrix positivity and rational
factorization also use only exact rational arithmetic.

The six-dimensional integral checker additionally uses python-flint's Arb
real intervals at 256 bits for square roots, $\pi$, exponentials,
$\operatorname{erf}$, and $\operatorname{erfc}$. Its required numerical
contract is outward containment of every exact result in the returned
interval, as specified in the [FLINT real-ball documentation](https://flintlib.org/doc/arb.html).
The conversions `upper()` and `lower()` produce outward exact dyadic
endpoints; `fmpq()` then transfers those endpoints to rational arithmetic,
as documented by [python-flint](https://python-flint.readthedocs.io/en/latest/arb.html).
Gaussian moments, range subdivision, tails, vertex corrections, coverage,
and final comparisons are separately specified formulas, not an appeal to
a numerical optimizer or an opaque integration routine. Correct execution
of these arithmetic operations is the implementation dependency of the
computer-assisted certificate.
