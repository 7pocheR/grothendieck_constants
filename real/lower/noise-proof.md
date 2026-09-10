# Joint and individual Gaussian noise constraints

Use the normalization and original and oddified midpoints in Section 1
of [MATHEMATICS.md](MATHEMATICS.md). In particular,
$h=(h_0)_o$, $k=(k_0)_o$, $|h|+|k|\le1$, and
$x=\|P_1h\|_2/\nu=\|P_1h_0\|_2/\nu$. The choice of $x$
does not change when switching between the two arguments.

This note proves $L_a/\nu^2<C=2917/2000$ on
$[0,11/50]$ from joint and individual Gaussian noise inequalities.
Together with the two ternary inequalities and the finite-matrix
construction in the main proof, it gives an additional complete proof
of $K_G^{\mathbb R}\ge1625\pi/2917$.

## 1. Noise normalization and a Gaussian comparison inequality

For $-1<\rho<1$, let $X,Z$ be independent standard Gaussian
vectors and put $Y=\rho X+\sqrt{1-\rho^2}Z$. Thus
$\mathbb EX_iY_j=\rho\delta_{ij}$. The covariance matrix has
eigenvalues $1+\rho$ and $1-\rho$, so is positive definite.
Define
$T_\rho q(X)=\mathbb E_Zq(\rho X+\sqrt{1-\rho^2}Z)$ and
$\mathcal S_\rho(q)=\langle q,T_\rho q\rangle$.
The Hermite generating identity proves
$T_\rho\psi_\alpha=\rho^{|\alpha|}\psi_\alpha$. Therefore
$$
 \mathcal S_\rho(q)=\sum_{j\ge0}\rho^j\|P_jq\|_2^2.
                                                               \tag{1}
$$
For $|\rho|\le1$, this series is absolutely convergent by
Parseval. It contains only odd terms when $q$ is odd. For positive
$\rho$ all its coefficients are nonnegative.

Let $J_\rho(p,q)$ be the probability that two standard scalar
Gaussians of correlation $\rho$ lie below
$\Phi^{-1}(p)$ and $\Phi^{-1}(q)$, respectively. For arbitrary
measurable $b,c:\mathbb R^n\to[0,1]$,
$$
 J_{-\rho}(\mathbb Eb,\mathbb Ec)
 \le\mathbb E b(X)c(Y)
 \le J_\rho(\mathbb Eb,\mathbb Ec),\qquad 0<\rho<1.     \tag{2}
$$
These are the bounded-function consequences of Borell's Gaussian
comparison theorem. The set and functional forms with this precise
covariance convention are Theorems 1.1 and 1.2 of
[Mossel and Neeman, Robust Optimality of Gaussian Noise Stability](https://arxiv.org/pdf/1210.4126).

A derivation also verifies the signs and the bounded-function extension.
For $u=\Phi^{-1}(p)$, $v=\Phi^{-1}(q)$, direct boundary
differentiation gives $J_{pq}>0$ and
$$
 J_{pp}=-\rho J_{pq}\frac{\phi(v)}{\phi(u)},\qquad
 J_{qq}=-\rho J_{pq}\frac{\phi(u)}{\phi(v)}.
$$
Consequently
$$
 \begin{pmatrix}J_{pp}&\rho J_{pq}\\\rho J_{pq}&J_{qq}\end{pmatrix}
 =\rho J_{pq}
 \begin{pmatrix}-\phi(v)/\phi(u)&1\\1&-\phi(u)/\phi(v)\end{pmatrix}
 \preceq0.                                             \tag{3}
$$
For the Ornstein--Uhlenbeck semigroup $P_t=T_{e^{-t}}$, consider
$R(t)=\mathbb E J_\rho(P_tb(X),P_tc(Y))$, where the outer
Gaussian pair retains correlation $\rho$. The joint stationary
generator is
$$
 \mathcal L=\Delta_X+\Delta_Y+2\rho\sum_i\partial_{X_i}\partial_{Y_i}
             -X\cdot\nabla_X-Y\cdot\nabla_Y.
$$
The chain rule and $\mathbb E\mathcal L=0$ show
$$
 R'(t)=-\mathbb E\{J_{pp}|\nabla P_tb|^2+J_{qq}|\nabla P_tc|^2
                  +2\rho J_{pq}\langle\nabla P_tb,\nabla P_tc\rangle\}
 \ge0                                                     \tag{4}
$$
by (3). For full justification, initially replace $b,c$ by
$\delta+(1-2\delta)b$, $\delta+(1-2\delta)c$. Their values
stay in $[\delta,1-\delta]$, so the derivatives of $J$ are
bounded on their value range. On any time interval bounded away from
zero, the Mehler formula gives bounded first and second spatial
derivatives; multiplication by $X,Y$ is Gaussian-integrable.
Gaussian integration by parts and differentiation under the expectation
in (4) are then justified. Strong $L^2$ continuity of $P_t$ at
zero and convergence to the mean at infinity follow from its Hermite
expansion. Uniform continuity of $J$ on $[0,1]^2$, followed
by $\delta\downarrow0$, gives
$$
 \mathbb E J_\rho(b(X),c(Y))\le J_\rho(\mathbb Eb,\mathbb Ec).
$$
The correlation derivative of $J$ is the positive bivariate
Gaussian density at its two finite thresholds. Since $J_0(p,q)=pq$,
$J_\rho(p,q)\ge pq$ for $\rho>0$, including boundary values
by continuity. This proves the upper inequality of (2). Applying that
upper inequality to $b,1-c$, and using
$p-J_\rho(p,1-q)=J_{-\rho}(p,q)$, proves the lower inequality.

## 2. The absolute-moment and joint inequalities

For $0<\tau<1$, define
$$
 F_\tau(m)=2\{J_\tau(m/2,m/2)-J_{-\tau}(m/2,m/2)\},
 \qquad 0\le m\le1.
$$
If $y$ is odd and $|y|\le B$, with $B>0$, its positive
and negative parts divided by $B$ both have mean
$m/2=\mathbb E|y|/(2B)$. Expanding its noise quadratic form and
applying (2) to the two positive terms and the two cross terms proves
$$
 \mathcal S_\tau(y)\le B^2F_\tau(\mathbb E|y|/B).        \tag{5}
$$
Equality is attained by the odd tail function
$B\operatorname{sgn}(X_1)\mathbf1_{|X_1|\ge r}$, with
$2\Phi(-r)=m$; its positive and negative tails attain the two
comparisons simultaneously. Endpoint masses have their continuous
interpretations. Dimension zero has only the zero odd function.

In particular, for any odd $|q|\le1$, applying (2) to
$(q+1)/2$, which has mean $1/2$, gives
$$
 \mathcal S_\rho(q)\le4J_\rho(1/2,1/2)-1
                         =\nu^2\arcsin\rho.
$$
The last identity follows by integrating the correlation derivative
$1/(2\pi\sqrt{1-\rho^2})$ from zero, where the probability is
$1/4$. Apply this to $f_o$ and $g_o$, and use the
parallelogram identity for the quadratic form $\mathcal S_\rho$:
$$
 \mathcal S_\rho(h)+\mathcal S_\rho(k)
 =\tfrac12(\mathcal S_\rho(f_o)+\mathcal S_\rho(g_o))
 \le\nu^2\arcsin\rho.                                  \tag{6}
$$
This does not require $f_o,g_o$ to be signs.

## 3. Global affine majorants and the exact inflection point

Put $c=\sqrt{(1-\tau)/(1+\tau)}$. If
$u=\Phi^{-1}(m/2)<0$, differentiating the two Gaussian boundaries
gives
$$
 F_\tau'(m)=2\{\Phi(cu)-\Phi(u/c)\},\qquad
 F_\tau''(m)=\frac{c\phi(cu)-c^{-1}\phi(u/c)}{\phi(u)}.  \tag{7}
$$
These derivatives are taken only at interior mass values. In particular
$F_\tau$ is nondecreasing. The second derivative changes sign
exactly when $u^2=r_*^2$, where
$$
 r_*^2=\frac{4\log(1/c)}{c^{-2}-c^2}
       =\frac{2\log(1/c^2)}{c^{-2}-c^2}.                 \tag{8}
$$
Thus $F_\tau$ is convex on $[0,2\Phi(-r_*)]$ and concave
on $[2\Phi(-r_*),1]$. These endpoint statements follow by continuity
from the interior derivative signs; endpoint differentiation is not used.

Fix $0\le r\le r_*$, let $m_r=2\Phi(-r)$, and set
$$
 d_1=2\{\Phi(r/c)-\Phi(cr)\},\qquad
 d_0=F_\tau(m_r)-m_rd_1.
$$
If $d_0\ge0$, then
$$
 F_\tau(m)\le d_0+d_1m\quad(0\le m\le1).              \tag{9}
$$
Concavity proves the tangent bound on the concave interval. On the
remaining convex interval, $F_\tau$ minus the tangent is convex,
so lies below the larger endpoint value. At its right endpoint the
tangent bound applies; at zero the difference is $-d_0\le0$.
This proves globality rather than just a local supporting property.

The exact integral for the tangent height is
$$
 F_\tau(2\Phi(-r))=
 \frac1\pi\int_{-\arcsin\tau}^{\arcsin\tau}
            \exp\left(-\frac{r^2}{1+\sin\theta}\right)d\theta.
                                                               \tag{10}
$$
To derive it, integrate the correlation derivative
$\exp(-r^2/(1+s))/(2\pi\sqrt{1-s^2})$ of the equal-threshold
Gaussian probability over $[-\tau,\tau]$, multiply by two, and
substitute $s=\sin\theta$. The density identity follows by integrating
its mixed spatial derivative over a lower quadrant. Gaussian boundary
terms vanish, and the compact correlation interval stays away from
$\pm1$, so a common integrable Gaussian bound justifies the
derivative and integral interchange.

For rational certificates it is sufficient to use rational
$D_0\ge d_0$, $D_1\ge d_1$, both nonnegative. More generally,
it suffices that
$D_0\ge d_0+\max(0,d_1-D_1)$, $D_1\ge0$.
The latter condition bounds the difference of the two affine functions
over the whole interval $[0,1]$. Rounding the slope upward and then
defining the intercept by subtracting that rounded slope times $m_r$
does not by itself majorize the original tangent to the left of $m_r$.

## 4. Combining distinct correlations

Let $0\le\alpha\le1$, $0<\rho,\tau<1$, and suppose
$F_\tau(m)\le D_0+D_1m$ globally, with $D_0,D_1\ge0$.
By (5), monotonicity, and $\mathbb E|k|\le1-\mathbb E|h|$,
$$
 \mathcal S_\tau(k)\le D_0+D_1(1-\mathbb E|h|).        \tag{11}
$$
Write $A_j=\|P_jh\|_2^2$, $B_j=\|P_jk\|_2^2$, and
$\lambda=\alpha/\rho^3$, $\mu=(1-\alpha)/\tau^3$.
Subtract $\lambda(\mathcal S_\rho(h)+\mathcal S_\rho(k))
+\mu\mathcal S_\tau(k)$ from the exact objective, then use
(6) and (11). The exact remainder is
$$
\begin{split}
 L_a\le{}&\lambda\nu^2\arcsin\rho+\mu(D_0+D_1)
            -\mu D_1\mathbb E|h|\\
 &+(a-\alpha\rho^{-2})A_1
   -(a+\alpha\rho^{-2}+(1-\alpha)\tau^{-2})B_1
   -(1+\alpha)A_3\\
 &-\sum_{j\ge5,\ j\text{ odd}}
    \{\alpha\rho^{j-3}A_j+
       [\alpha\rho^{j-3}+(1-\alpha)\tau^{j-3}]B_j\}.
                                                               \tag{12}
\end{split}
$$
The degree-three coefficient of $k$ is exactly zero. Every other
$k$ coefficient is nonpositive, including degree one. Every
$h$ coefficient above degree one is nonpositive. All these series
converge absolutely: for fixed $\rho,\tau>0$, their coefficient
magnitudes are bounded and $\sum_j(A_j+B_j)<\infty$. Thus the
sign argument includes the entire infinite Hermite remainder.

Let $J$ be a finite set of odd degrees containing three. For
$\alpha>0$, put
$c_3=1+\alpha$, $c_j=\alpha\rho^{j-3}$ for $j\ge5$.
For $\alpha=0$, take $J=\{3\}$. Retain the corresponding
negative $h$ terms in (12). With $S$ along $P_1h$,
$m_j=\mathbb E h\psi_j(S)$ obeys $m_j^2\le A_j$.
For arbitrary real $z_j$,
$$
 -c_jA_j\le-c_jm_j^2
 \le\nu^2z_j^2/c_j-2\nu z_jm_j.
$$
Add the exact zero $\nu\ell(\nu x-\mathbb E Sh)$. At each
point maximize the remaining expression over $|h|\le1$, using
$\sup_{|t|\le1}(At-B|t|)=(|A|-B)_+$ for $B\ge0$.
This proves, universally in every finite dimension and every sign pair,
$$
\begin{split}
 L_a/\nu^2\le{}&
 \frac{\alpha\arcsin\rho}{\rho^3}
 +\frac{(1-\alpha)(D_0+D_1)}{\nu^2\tau^3}
 +(a-\alpha/\rho^2)x^2+\ell x+\sum_{j\in J}z_j^2/c_j\\
 &+\frac1\nu\mathbb E
 \left(\left|\ell S+2\sum_{j\in J}z_j\psi_j(S)\right|
              -\frac{(1-\alpha)D_1}{\nu\tau^3}\right)_+.
                                                               \tag{13}
\end{split}
$$
The inequality holds for every fixed choice of parameters. Consequently
different fixed choices may be used on different $x$-intervals,
or their upper bounds may be minimized pointwise. The proof does not
assume that a numerical minimization finds any optimum. Both correlations
are scalar parameters of valid inequalities on the same actual functions;
there is no requirement to couple two Gaussian pairs with those different
correlations simultaneously.

## 5. Finite verification of a scalar support

For any chosen rational $\rho,\tau,\alpha,\ell,(z_j)$, a finite
scalar certificate for (13) must specify a closed rational $x$-interval
and prove all of the following:

1. $0<\rho,\tau<1$, $0\le\alpha\le1$, the allowed retained
   degrees, and positive denominators $c_j$.
2. A global affine majorant for $F_\tau$, for example using a
   rational tangent location $r$, the exact threshold (8),
   a nonnegative tangent intercept, and rigorous enclosures of (10),
   its slope, and the rounded affine coefficients.
3. A complete rational partition of a finite Gaussian interval for the
   positive-part integral, keeping an upper contribution on every cell
   whose activity or sign is unresolved. For degree $d$, the exact
   Hermite polynomial formula is
   $H_j(s)=j!\sum_{k=0}^{\lfloor j/2\rfloor}
     (-1)^ks^{j-2k}/(2^kk!(j-2k)!)$.
   Sign-certified cells are integrated by
   $\int\psi_j\phi=-\phi H_{j-1}/\sqrt{j!}$, and the constant
   threshold is integrated by Gaussian mass.
4. The full tail bound $\sum_{j=0}^d|b_j|J_j$ for the polynomial
   $\sum b_js^j$, with the recurrence for $J_j$ in Section 4
   of the main note. This also bounds the positive-part tail, since
   $(|P|-B)_+\le|P|$ for $B\ge0$.
5. A rigorous rational upper bound on the total constant in (13).
   The resulting quadratic in $x$ must be bounded over the entire
   specified interval. If its quadratic coefficient is nonnegative,
   endpoint checks suffice. If it is negative, check its vertex as well
   whenever $-\ell/[2(a-\alpha/\rho^2)]$ lies in the interval.
6. Exact coverage by the union of all intervals used, including all
   junctions with the other universal estimates.

No bounds uniform over a search domain of $\rho,\tau,\ell,z$
are required when the certificate specifies finitely many fixed choices.
The universal quantifiers are over the functions, dimension, and assigned
$x$-interval, not over a numerical optimizer's parameter domain.
There is also no required limit as either correlation tends to an
endpoint; every selected correlation lies strictly between zero and one.

## 6. Three explicit affine majorants

The following choices satisfy (9) with rational upper coefficients:
$$
\begin{array}{c|cc|cc}
 &\tau&r&D_0&D_1\\\hline
 1&97/100&1/10&2928311/10^7&114483/200000\\
 2&46029461/46913053&8052339/78719627&20425809/10^8&69791401/10^8\\
 3&22328505/22624636&12172859/98597455&425633/6250000&43188827/50000000
\end{array}
$$
Their exact tangent intercepts lie respectively in the rational intervals
$$
 [0.29283095,0.29283096],\quad
 [0.20425806,0.20425807],\quad
 [0.06810125,0.06810126],
$$
and their exact slopes lie respectively in
$$
 [0.57241486,0.57241487],\quad
 [0.69791399,0.69791400],\quad
 [0.86377652,0.86377653].
$$
All endpoints in these displays are terminating decimal rationals.
For each choice, $r\le1/8$ and (8) gives $r_*^2>1/16$,
so its tangent is on the concave branch. The displayed intervals prove
nonnegative intercepts and that $D_0,D_1$ dominate the exact tangent
coefficients.

Here is a finite integration rule and error proof for these enclosures.
Put $\alpha=\arcsin\tau$. Each choice obeys
$\alpha<143/100$. On every complex disk of radius $1/50$
centered on $[-\alpha,\alpha]$, the real part is between
$-1.45$ and $1.45$. For $z=u+iv$, the exact identity
$$
 |1+\sin z|=\cosh v+\sin u
$$
and $\pi>3.14$ imply
$$
 |1+\sin z|\ge1-\sin1.45>1-\cos0.12>0.12^2/3=3/625.
$$
The last comparison follows from the alternating cosine series. The
integrand $f(z)=\exp(-r^2/(1+\sin z))$ is analytic on these
disks and has modulus less than $e^4<55$. Cauchy's estimate gives
$|f^{(8)}|\le8!\,55\,50^8$ on the real interval.

Partition $[-\alpha,\alpha]$ into $N=1024$ equal intervals.
On the reference interval $[-1,1]$, use the four nodes
$$
 \pm\sqrt{\frac{3\pm(2/5)\sqrt{30}}7},
 \qquad w_\pm=\frac{18\mp\sqrt{30}}{36}.
$$
The node with the corresponding squared value receives weight
$w_\pm$, as does its negative. Their odd moments vanish, and their
moments of degrees $0,2,4,6$ equal $2,2/3,2/5,2/7$,
respectively. Thus this quadrature integrates every polynomial of degree
at most seven exactly. Its monic node polynomial is
$t^4-(6/7)t^2+3/35$, whose squared integral is $128/11025$.

Hermite interpolation of $f$, matching both its values and first
derivatives at the four nodes, bounds the quadrature error on an interval
of length $h$ by
$$
 \frac{\sup|f^{(8)}|}{8!}\frac{h^9}{44100}.
$$
Indeed the interpolation remainder is the eighth derivative divided by
$8!$, times the squared monic node polynomial; the scaled squared
integral is $h^9/44100$. Integrating the degree-seven interpolant
equals its quadrature, which equals the quadrature of $f$.
For (10), summing all panels and dividing by $\pi$ gives the
following rational total error bound:
$$
 E=\frac{55\,50^8}{44100}\left(\frac3{1024}\right)^8.
$$
Here $2\alpha<3$ and $\pi>3$. Every node, weight, panel,
and final sum is enclosed outward. The finite panel bounds and the
intercept and slope enclosures are specified with the other rational
data for the three cases in
[scalar_1.json](certificates_02/scalar_1.json),
[scalar_2.json](certificates_02/scalar_2.json), and
[scalar_3.json](certificates_02/scalar_3.json). Each file contains all
1,024 panel bounds. Their exact lower and upper sums, minus and plus
$E$, respectively, enclose (10). Subtracting the enclosed mass
times the enclosed slope proves the displayed tangent bounds.

## 7. Three complete positive-part certificates

Use the three $(\tau,r,D_0,D_1)$ choices above. The other
parameters in (13) are
$$
\begin{array}{c|ccc}
 &\rho&\alpha&\ell\\\hline
 1&97/100&1/2&-2674301/10^7\\
 2&64951013/67020094&21568362/56492747&-309787492/780329499\\
 3&5050709/5213139&23085713/60976975&-458408705/945444501
\end{array}
$$
For case 1, take $J=\{3,5,7,9\}$ and
$$
 (z_3,z_5,z_7,z_9)
 =(416616,-171151,87194,-111814)/10^7.
$$
For cases 2 and 3, take $J=\{3,5,7,9,11,13,15,17\}$, with
$$
\begin{array}{c|rr}
 j&z_j\text{ in case 2}&z_j\text{ in case 3}\\\hline
 3&35354125/576596594&71528951/959225474\\
 5&-17166393/598918432&-27808442/799018179\\
 7&2472213/134935441&20576857/928627584\\
 9&-7515781/768625502&-11399183/983152295\\
 11&4438747/863585656&3614385/595241402\\
 13&-2798771/893877344&-1764117/456886864\\
 15&1752971/715598683&2388933/802063733\\
 17&-1112959/622056980&-1641267/928239158
\end{array}
$$
Every displayed parameter is rational, both correlations lie strictly
between zero and one, and $0<\alpha<1$. All coefficients $c_j$
used in (13) are positive.

For each case, let
$$
 P(s)=\ell s+2\sum_{j\in J}z_j\psi_j(s),\qquad
 B=\frac{(1-\alpha)D_1}{\nu\tau^3}>0,\qquad
 I=\mathbb E(|P(S)|-B)_+.
$$
The `cells` field in its data file lists triples $(l,u,U_{l,u})$
of exact rationals. These are full consecutive partitions of $[0,12]$:
case 1 has 2,413 cells, and cases 2 and 3 have 60,001 cells each.
The value $U_{l,u}$ bounds
$\int_l^u(|P(s)|-B)_+\phi(s)\,ds$.

These bounds are evaluated at 224-bit precision by the following finite
procedure, implemented in `scalar_calculator` in
[check_scalars.py](check_scalars.py). Convert the finitely many Hermite
polynomials to monomials by the exact formula in Section 5. On each
cell, enclose $P$ by its complete midpoint Taylor expansion and
the sum of absolute nonconstant terms. If the range is within
$[-B,B]$, the integral is zero. If the whole range is above $B$
or below $-B$, integrate respectively $P-B$ or $-P-B$
by the Gaussian monomial recurrence. For every other cell, retain the
upper bound
$$
 \bigl(\sup_{[l,u]}|P|-B\bigr)_+\{\Phi(u)-\Phi(l)\},
$$
with an upper bound on the supremum and a lower bound on $B$.
All range tests and arithmetic are outward, and each final cell bound
is rounded upward to a multiple of $2^{-110}$.

The polynomial $P$ is odd, so the integrand is even. Thus the
full finite-interval upper bound is exactly twice the sum of all cell
bounds. Add the entire two-sided tail at 12 using
$\sum_j|b_j|J_j$ for the monomial coefficients of $P$,
as in Section 5. This gives
$$
\begin{split}
 I_1&\le\frac{238660642243234137343187628781}{2^{110}},\\
 I_2&\le\frac{62107970551708416360755113583}{2^{110}},\\
 I_3&\le\frac{136387505281180401733677637758}{2^{110}}.
\end{split}
$$
The included tail upper bounds, with the same denominator $2^{110}$,
have numerators 1243179, 5849605579, and 5949439754, respectively.
All three nonzero tail bounds are included in the displayed values.

Assemble the constant term
$$
 K=\frac{\alpha\arcsin\rho}{\rho^3}
   +\frac{(1-\alpha)(D_0+D_1)}{\nu^2\tau^3}
   +\sum_{j\in J}\frac{z_j^2}{c_j}+\frac I\nu.
$$
Outward arithmetic and rational upward rounding give the simpler upper
bounds
$$
 K_1\le\frac{1473083457}{10^9},\qquad
 K_2\le\frac{1488416356}{10^9},\qquad
 K_3\le\frac{1505722167}{10^9}.                         \tag{14}
$$
The exact rational sums and the stronger intermediate constant bounds
are also in the three data files. The quadratic coefficients
$b=a-\alpha/\rho^2$ are
$$
 b_1=\frac{82317}{75272},\quad
 b_2=\frac{2323161283219143570418103}{1906577826531806116769944},\quad
 b_3=\frac{15202356005025750366691}{12444015884877954159800}.
$$
All three are positive. Hence the maximum of
$b_i x^2+\ell_i x+K_i$ on each closed interval occurs at an
endpoint. Exact rational endpoint arithmetic using (14) gives
$$
\begin{array}{c|c|c}
 i&x\text{ interval}&C-\max(b_i x^2+\ell_i x+K_i)\\\hline
 1&[1/10,4/25]&>1/5000\\
 2&[4/25,1/5]&>7/10000\\
 3&[1/5,11/50]&>3/10000
\end{array}                                                     \tag{15}
$$
These comparisons cover every point of these intervals, including their
junctions. They follow from the same parameters that define each
integral and its tangent, with no parameter optimization assumption.

On $[0,1/10]$, use only the joint constraint (6), with
$\rho=97/100$. Taking $\alpha=1$ in (12) and discarding
its explicitly nonpositive remainder yields
$$
 L_a/\nu^2\le\frac{\arcsin\rho}{\rho^3}
                   +(a-\rho^{-2})x^2
 \le\frac{58306173}{40000000}
 =C-\frac{33827}{40000000}<C.                           \tag{16}
$$
Here $a-\rho^{-2}>0$, so the last comparison is an upper
enclosure at $x=1/10$. Equations (15)--(16) prove the desired
bound on the entire interval $[0,11/50]$.

All 122,415 finite cells, the three full Gaussian tails, and all 3,072
quadrature panels enter the finite verification. The quadrature formula
and its analytic remainder are implemented in `noise_certificate` in
[check_scalars.py](check_scalars.py), at 224-bit precision with rational
rounding to multiples of $2^{-110}$. The implementation dependencies
and outward arithmetic requirements are those in Section 7 of the main
proof. [check_certificate_arithmetic_02.py](check_certificate_arithmetic_02.py)
checks exact interval coverage, rational sums, tangent comparisons, and
the displayed endpoint margins from the mathematical data. Combining
this small-first-moment bound with $(T_A)$, $(T_B)$, and
Sections 5--6 of the main proof establishes the same finite matrices
and the same lower bound for $K_G^{\mathbb R}$.
