# An upper bound for the complex Grothendieck constant

For every pair of positive integers $m,n$, every matrix $A=(a_{ij})\in\mathbb C^{m\times n}$, every finite-dimensional complex Hilbert space, and all unit vectors $u_i,v_j$ in that space,

$$
\left|\sum_{i,j}a_{ij}\langle u_i,v_j\rangle\right|
\leq \frac{100000000}{71188883}
\sup_{|\alpha_i|=|\beta_j|=1}
\left|\sum_{i,j}a_{ij}\alpha_i\overline{\beta_j}\right|.
\tag{1}
$$

The inner product is linear in its first argument. Thus the classical unrestricted complex Grothendieck constant satisfies

$$
K_G^{\mathbb C}\leq\frac{100000000}{71188883}
=1.4047137107067686\ldots.
$$

The argument below reduces (1) to explicit exact data and outward interval computations. [candidate.json](data/candidate.json) specifies all rational parameters; [plan.json](data/plan.json) specifies the retained modes and their radii and degrees. The numerical part uses the conventional correctness of Python integer arithmetic and FLINT/Arb outward arithmetic. The distinction between checking stored data and recomputing native enclosures is explained in the [reproduction instructions](README.md).

## 1. Gaussian correlations and the coefficient criterion

A standard proper complex Gaussian variable has density $\pi^{-1}e^{-|z|^2}$ and $\mathbb E|Z|^2=1$. Let $W,X$ be independent such variables. For two functions $f,g\in L^2(\mathbb C^2)$, consider their correlation when the only cross covariances are

$$
\mathbb E W\overline{W'}=a,\qquad
\mathbb E X\overline{X'}=b.
$$

For real $a,b\in(-1,1)$, expansion in an orthonormal Gaussian Hermite basis gives

$$
C_{f,g}(a,b)=\sum_{p,q\geq0}h_{pq}a^pb^q,
\qquad
\sum_{p,q}|h_{pq}|\leq\|f\|_2\|g\|_2.
\tag{2}
$$

To see the inequality, expand $f$ and $g$ in the same basis. The correlation of matching basis elements is $a^pb^q$, where $p,q$ are their total degrees in the two complex coordinates, and unmatched elements have zero correlation. Group the products of basis coefficients by $(p,q)$, apply the triangle inequality, and then Cauchy–Schwarz to the full coefficient sequences. This argument also proves continuity and sesquilinearity of the correlation map into the Banach space with norm $\sum|h_{pq}|$.

Let $P(x),B(x)$ be real odd polynomials whose absolute coefficient sums are at most one. Substitution into (2) is a contraction for the absolute coefficient norm. In particular it defines an absolutely summable real-variable series. If both functions have angular character one and commute with complex conjugation, their composed correlation has real coefficients and is odd:

$$
C_{f,g}(P(x),B(x))=\sum_{j\geq0}c_jx^{2j+1}.
\tag{3}
$$

The *physical* correlation for complex input $z$ is

$$
\sum_{j\geq0}c_jz|z|^{2j}.
\tag{4}
$$

Indeed its values on real inputs are (3), and rotating both cross covariances by $e^{i\phi}$ multiplies the correlation by $e^{i\phi}$. Equation (4) follows by taking $z=e^{i\phi}|z|$. It is distinct from the analytic auxiliary series $\sum c_jz^{2j+1}$ that will be used for Cauchy estimates.

The radial monomial $z|z|^{2j}$ has the tensor realization

$$
u\longmapsto u^{\otimes(j+1)}\otimes\overline u^{\otimes j}.
$$

For an odd polynomial $P(x)=\sum_jp_jx^{2j+1}$, take direct sums of these tensors with factors $\sqrt{|p_j|}$, and put the sign of $p_j$ in the column vectors. Their inner products are $P_{\rm rad}(\langle u,v\rangle)$. Add row-only and column-only orthogonal components to make every norm one. The hypothesis $\sum|p_j|\leq1$ makes this possible. Apply the same construction to $B$ in a separate space; independent Gaussian vectors on the two spaces give independent $W,X$.

For a finite matrix $A$, write $S(A)$ for the scalar supremum in (1), and $V(A)$ for the supremum over unit vectors of the absolute matrix value. A vector maximizer may be chosen in dimension at most $m+n$, by its Gram matrix, and a common rotation makes its matrix value positive real. If $|f|=|g|=1$, the expected rounded scalar value has modulus at most $S(A)$. Each nonlinear term in (4) has absolute matrix value at most $V(A)$ by its tensor realization. Absolute convergence permits termwise summation, including correlations on the unit circle. Therefore

$$
\Gamma V(A)\leq S(A),\qquad
\Gamma=c_0-\sum_{j\geq1}|c_j|.
\tag{5}
$$

A lower bound $\Gamma>71188883/100000000$ proves (1). No restriction on the entries of $A$, matrix dimensions, or complex vector dimension occurs in this implication.

## 2. The exact phase functions

Put $q(w,x)=\operatorname{Im}(x\overline w)$ and

$$
E_d(w,x)=d_w|w|^2+d_x|x|^2,\qquad
\Psi(w,x)=\sum_{l=1}^5\varepsilon_l e^{-\eta_l E_d(w,x)}\sin(2\lambda_lq(w,x)).
$$

Use

$$
f(w,x)=\operatorname{ph}(w+\theta x)e^{-2ikq(w,x)+i\Psi(w,x)},
\qquad g(w,x)=f(w,-x).
\tag{6}
$$

The value of the phase at zero can be assigned arbitrarily; that event has Gaussian measure zero. Both functions have modulus one, angular character one, and commute with conjugation. The rational input gives

$$
\begin{aligned}
\theta&=\frac{15869581262852117}{50000000000000000},&
k&=\frac{465518641186931}{250000000000000000},\\
d_w&=\frac{9209367321696771}{10000000000000000},&
d_x&=\frac{2171701844586119}{2000000000000000}.
\end{aligned}
$$

The five triples $(\varepsilon_l,\lambda_l,\eta_l)$ and every coefficient of $P,B$ are specified exactly in `candidate.json`. In particular

$$
\|P\|_1=\frac{9999999460789019878400803}{10000000000000000000000000}<1,
\qquad \|B\|_1=1,
$$

with

$$
B(x)=\frac{4994394845069840469}{5000000000000000000}x^3
-\frac{5605154930159531}{5000000000000000000}x^5.
$$

The degree of $P$ is 159. These finite rational inputs completely determine the construction.

## 3. Complete phase expansion and omitted modes

Set $\ell=999/1000$. Exact arithmetic verifies $d_wd_x\geq\ell^2$. Since

$$
E_d\geq2\sqrt{d_wd_x}|q|\geq2\ell|q|,
$$

and $\sup_{t\geq0}te^{-t}=1/e<2/5$, one has $\|\Psi\|_\infty\leq M$, where

$$
M=\sum_l|\varepsilon_l|
\begin{cases}
1,&\eta_l=0,\\
\min\{1,2|\lambda_l|/(5\eta_l\ell)\},&\eta_l>0,
\end{cases}
=\frac{19007117794841543023}{399600000000000000000}.
\tag{7}
$$

Replace $\Psi$ by $t\Psi$ in both functions in (6), and view their correlation after preprocessing as an element of the absolute coefficient space. Five differentiations, the product rule, and (2) bound its fifth derivative for real $t$ by

$$
\sum_{j=0}^5\binom5j\|\Psi^j\|_2\|\Psi^{5-j}\|_2\leq(2M)^5.
$$

Integral Taylor remainder at $t=1$ therefore bounds the complete error after **combined** row and column phase order four by

$$
T_{\rm phase}=\frac{(2M)^5}{120}
=6.49266007142443\ldots\,10^{-8}.
\tag{8}
$$

This bounds all larger phase orders and all scalar degrees at once.

For each side and each profile introduce two generators. A generator of sign $\sigma\in\{-1,1\}$ has weight $\sigma\varepsilon_l/2$, frequency increment $\sigma\lambda_l$, and damping increment $\eta_l$. A multiset of at most four generators has weight equal to the product of generator weights divided by every multiplicity factorial. If its accumulated frequency increments are $t_1,t_2$ and damping increments are $u,v$, its primitive key is

$$
(r,u,s,v)=(k-t_1,u,k-t_2,v).
\tag{9}
$$

The frequency signs in (9) follow from $i\sin(2\lambda q)=(e^{2i\lambda q}-e^{-2i\lambda q})/2$. The column reflection changes $q'$ to $-q'$, and conjugation changes its sign again. Thus both primitive quadratic phases are negative, and both phase perturbations in the correlation contribute with the same positive sign.

Merge equal keys and exchange the row and column pairs into lexicographic order. The primitive is invariant under this exchange, as verified by the determinant identities below. The result consists of exactly 2,962 nonzero rationally weighted modes. Two different exact constructions, direct multiset enumeration and repeated polynomial multiplication divided by the next order, are included in the checkers.

The schedule retains 2,204 modes and omits 758. For a mode of damping units $u,v$, the two damped phase functions have $L^2$ norm product

$$
\frac1{\sqrt{(1+2ud_w)(1+2ud_x)(1+2vd_w)(1+2vd_x)}}
\leq\frac1{(1+2u\ell)(1+2v\ell)}.
\tag{10}
$$

The inequality follows by expanding each pair and using $d_w+d_x\geq2\ell$ and $d_wd_x\geq\ell^2$. By (2), the omitted mode's entire coefficient norm is at most its absolute weight times (10). Summing over the 758 exact omitted indices gives

$$
T_{\rm omit}=\sum_{i\text{ omitted}}
\frac{|w_i|}{(1+2u_i\ell)(1+2v_i\ell)}
=1.3111185111285402\ldots\,10^{-7}.
\tag{11}
$$

No analytic continuation is needed for any omitted mode.

## 4. The Gaussian primitive and its normalization

For a retained key $(r,u,s,v)$, let the row coordinate dampings be $a_w=ud_w,a_x=ud_x$ and the column dampings be $b_w=vd_w,b_x=vd_x$. The primitive is the expectation of

$$
\operatorname{ph}(W+\theta X)\overline{\operatorname{ph}(W'-\theta X')}
 e^{-2irq(W,X)-2isq(W',X')}
 e^{-uE_d(W,X)-vE_d(W',X')}.
$$

Initially $a,b$ are real correlations near zero. With coordinate order $W,X,W',X'$, put

$$
S=\begin{pmatrix}1&0&a&0\\0&1&0&b\\a&0&1&0\\0&b&0&1\end{pmatrix},\quad
J=\begin{pmatrix}-a_w&-r&0&0\\r&-a_x&0&0\\0&0&-b_w&-s\\0&0&s&-b_x\end{pmatrix},\quad
U=\begin{pmatrix}1&\theta&0&0\\0&0&1&-\theta\end{pmatrix}.
$$

Define the exact polynomials

$$
D=\det(I-SJ),\qquad
\begin{pmatrix}N_1&P_0\\Q_0&N_2\end{pmatrix}
=U\operatorname{adj}(I-SJ)S U^T.
\tag{12}
$$

The file `primitive_polynomials.json` gives their integer coefficients in nine variables. `check_algebra.py` derives (12) by determinant and matrix multiplication, and checks all five polynomials and row/column exchange. Their term counts are 48, 42, 42, 26, 26.

Let $C=(I-SJ)^{-1}S$ and $UCU^T=\left(\begin{smallmatrix}n_1&p\\q_1&n_2\end{smallmatrix}\right)$. The Laplace identity

$$
\frac z{|z|}=\frac z{\sqrt\pi}\int_0^\infty t^{-1/2}e^{-t|z|^2}\,dt
$$

applied to both phase factors, followed by Gaussian integration with two nonnegative rank-one precision updates, gives

$$
\frac p{\pi D}\int_0^\infty\!\int_0^\infty
\frac{(xy)^{-1/2}\,dx\,dy}
 {[(1+xn_1)(1+yn_2)-xypq_1]^2}.
$$

Fubini is justified by absolute integrability: after taking absolute values and performing the two scalar Laplace integrals, the phase denominators cancel their numerators, leaving an integrable quantity at most one. The original Gaussian precision minus $J$ is accretive, since the damping is nonnegative. The Gaussian determinant formula and precision updates are therefore valid.

Near the origin, expand in $pq_1$ and integrate the two beta integrals. This gives

$$
\frac{\pi p}{4D\sqrt{n_1n_2}}
\,{}_2F_1\!\left(\tfrac12,\tfrac12;2;\frac{pq_1}{n_1n_2}\right)
=
\frac{\pi P_0}{4D\sqrt{N_1}\sqrt{N_2}}
\,{}_2F_1\!\left(\tfrac12,\tfrac12;2;Z\right),
\quad Z=\frac{P_0Q_0}{N_1N_2}.
\tag{13}
$$

All square roots start at their positive values at zero. At zero, the damping hypotheses give $D,N_1,N_2>0$, while $P_0=Q_0=0$. Formula (13) defines the local analytic germ. Its continuation beyond the unit disk is proved separately in the next section. In the unmixed undamped specialization, the scalar coefficients are $\pi/4$ times $1,1/8,3/64,25/1024,\ldots$, providing a direct normalization check.

## 5. Full-circle enclosures and branches

Substitute the analytic polynomials $a=P(z),b=B(z)$ in (12). For each retained mode, its specified radius satisfies $R>1$. The circle is parameterized by $z=R e^{i\pi t}$, $0\leq t\leq2$, with 1,024 initial panels and dyadic subdivisions.

For either input polynomial $p(z)=\sum c_nz^n$, the second angular derivative on that circle has modulus at most

$$
L_2(p,R)=\sum_n n^2|c_n|R^n.
$$

On an arc $[a,b]$, let $z_0=R e^{i\pi(a+b)/2}$ and $h=(22/7)(b-a)/2$. The complete image of the arc is contained in the interval evaluation of

$$
p(z_0)+\left(i\sum_n n c_nz_0^n\right)[-h,h]
+[-h^2L_2/2,h^2L_2/2]+i[-h^2L_2/2,h^2L_2/2].
\tag{14}
$$

Taylor's theorem gives (14), since $\pi<22/7$. Native Arb evaluation encloses its midpoint, derivative, and arithmetic. Substitution into the finite polynomials gives rectangles for $D,N_1,N_2,P_0,Q_0$, and division gives a rectangle for $Z$. All stored endpoints are rounded outward. The endpoints $t=0$ and $t=2$ use exactly the same point $R$.

Every accepted arc rectangle for $D,N_1,N_2$ excludes zero. A convex compact rectangle excluding zero lies in an open half-plane through the origin; the argument of a curve within it varies in an interval of width less than $\pi$. Native endpoint evaluations also enclose the three quotients of the right endpoint value by the left endpoint value. Each quotient has positive real part, and its imaginary-to-real ratio lies in $[-1/2,1/2]$. The argument change on the arc is consequently the arctangent of that ratio, with no extra multiple of $2\pi$.

For $|x|\leq1/2$, the exact checker uses

$$
\arctan x=\sum_{j=0}^{11}\frac{(-1)^jx^{2j+1}}{2j+1}+E,
\qquad |E|\leq\frac{|x|^{25}}{25}.
$$

Summing the resulting outward intervals over the full partition puts each of the three total changes in $(-3,3)\subset(-\pi,\pi)$. Each is an integer multiple of $2\pi$, so it is zero. The argument principle proves that $D,N_1,N_2$ have no zero throughout the closed disk. Analytic logarithms and square roots are then uniquely fixed by their positive origin values.

For any complex $Z$, uniformly for $0\leq t\leq1$,

$$
|1-tZ|\geq
\max\!\left\{\min(1,1-\operatorname{Re}Z),\frac{|\operatorname{Im}Z|}{|Z|}\right\}.
\tag{15}
$$

The second term is assigned zero when $Z=0$. The first follows from the real part, and the second by minimizing the squared distance to the full line through $Z$. The rectangle arithmetic gives a strictly positive lower bound $\delta$ for (15) on every arc.

The polynomial family $N_1N_2-tP_0Q_0$, $0\leq t\leq1$, therefore has no boundary zero. At $t=0$ it has no interior zero; its root count remains zero throughout this homotopy. It follows that $Z$ avoids the slit $[1,\infty)$ throughout the disk. Formula (13) thus extends analytically to a neighborhood of the closed disk, with the roots fixed at zero and the hypergeometric branch fixed on the slit domain. No principal square root of a winding product is assumed on the contour.

Euler's normalized beta integral is

$$
{}_2F_1(\tfrac12,\tfrac12;2;Z)
=\frac2\pi\int_0^1t^{-1/2}(1-t)^{1/2}(1-tZ)^{-1/2}\,dt.
$$

The weight integrates to one, so its modulus on an arc is at most $\delta^{-1/2}$. An entire-arc bound for the primitive is therefore

$$
\frac{\pi|P_0|}{4|D|\sqrt{|N_1||N_2|}\sqrt\delta}.
\tag{16}
$$

The exact checker derives (16) from the rectangle endpoints using $\pi<22/7$ and integer square-root bounds. It reconstructs each dyadic panel from its index and subdivision path, checks ordered coverage of $[0,2]$ without gaps or overlaps, and recomputes every minimum, argument sum and modulus maximum. The complete supplied computation has 2,305,576 accepted arc leaves.

## 6. Scalar coefficients and complete scalar tails

Write $y=z^2$, $P(z)=zA(y)$ and $B(z)=zB_0(y)$. The polynomials $D,N_1,N_2$ are even and $P_0,Q_0$ are odd. After extracting the latter factors of $z$, formula (13) is computed as a real formal series in $y$. Its denominators have nonzero constants, $N_1(0)N_2(0)>0$, and $Z(0)=0$. The hypergeometric coefficients used are exactly

$$
\frac{\binom{2j}{j}^2}{16^j(j+1)}.
$$

For degree parameter $d$, coefficients through $y^d$ require no later hypergeometric terms. Inversion and the positive formal square root agree with the branch established above. The native computation uses 768-bit Arb arithmetic and computes exactly 792,700 coefficients across the 2,204 retained modes.

Let $M_i$ be the maximum of (16) over the entire circle for mode $i$, and $d_i$ its scalar degree parameter. Oddness and Cauchy's estimate give the complete remaining coefficient norm

$$
T_i=\frac{M_iR_i^{-2d_i-3}}{1-R_i^{-2}}.
\tag{17}
$$

The weighted sum of (17) is

$$
T_{\rm scalar}=\sum_{i\text{ retained}}|w_i|T_i
=5.34428752849829\ldots\,10^{-12}.
$$

Different radii and degree parameters are allowed: finite coefficient vectors are padded with zeros, and each omitted coefficient is already charged in its mode's full tail. Negative mode weights reverse interval endpoints and use absolute weights in the tail sum.

## 7. Exact final comparison

The bundled coefficient intervals are outward dyadic intervals of the form $[L,L+W]2^{-96}$, $W\geq0$. Each contains its original native interval. Exact integer arithmetic aggregates the signed weights, rounding each contribution outward at 256 bits. This compression changes the margin by far less than $10^{-18}$; the comparison is checked exactly, without relying on that decimal statement.

If $[a_j,b_j]$ encloses the finite signed coefficient sum, the verified lower bound is

$$
\Gamma\geq
 a_0-\sum_{j\geq1}\max(|a_j|,|b_j|)
 -T_{\rm scalar}-T_{\rm phase}-T_{\rm omit}.
\tag{18}
$$

The functional $c_0-\sum_{j\geq1}|c_j|$ is 1-Lipschitz in the complete absolute coefficient norm. Thus every complete error is subtracted once in (18), accounting for both possible linear and nonlinear contributions.

The values, displayed here only for orientation, are

$$
\begin{aligned}
a_0&\approx0.7118890497808265,\\
\sum_{j\geq1}\max(|a_j|,|b_j|)&\approx1.38835525249104\,10^{-8},\\
T_{\rm scalar}&\approx5.34428752849829\,10^{-12},\\
T_{\rm phase}&\approx6.49266007142443\,10^{-8},\\
T_{\rm omit}&\approx1.3111185111285402\,10^{-7}.
\end{aligned}
$$

The standard-library verifier computes (18) as a rational number and checks

$$
\Gamma>\frac{71188883}{100000000}.
$$

Together with (5), this proves (1). A full reproduction consists of the exact polynomial and input checks, the actual outward native computation of every coefficient and complete arc, the complete rational coverage and branch validation, and this final exact comparison. Hashes or arbitrary claimed rectangles alone do not provide the native enclosure step.
