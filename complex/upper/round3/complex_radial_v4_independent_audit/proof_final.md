# A computer-assisted upper bound for the complex Grothendieck constant

For the standard convention with complex scalar phases and complex unit
vectors, we prove
\[
 K_G^{\mathbb C}\leq\frac{100000000}{71185999}
 <\frac{5000000}{3559299}
 =\frac1{0.7118598}.
\]
The proof uses an explicit Gaussian phase rounding, signed tensor
preprocessing, and finite interval certificates. It controls all scalar
degrees and all phase orders. The finite numerical assertions and their
verification rules are specified below; the mathematical certificate is
provided in [certificate_final.json](certificate_final.json), with its
complete covering data referenced there.

**Convention and statement.** The complex inner product is linear in its first
argument. For a finite complex matrix \(M=(m_{ij})\), let

\[
 S(M)=\max_{|\alpha_i|=|\beta_j|=1}
 \left|\sum_{ij}m_{ij}\alpha_i\overline{\beta_j}\right|,
 \qquad
 V(M)=\max_{\|u_i\|=\|v_j\|=1}
 \left|\sum_{ij}m_{ij}\langle u_i,v_j\rangle\right|.
\]

The vector maximum is over complex Hilbert spaces; dimension at most the
number of row and column vectors suffices by Gram matrix factorization.
Define \(K_G^{\mathbb C}=\sup_{M\ne0}V(M)/S(M)\). Using an unconjugated
column scalar in the definition gives the same maximum by replacing each
column scalar with its conjugate.

**Theorem.** For every finite complex matrix \(M\),
\[
 \frac{71185999}{100000000}\,V(M)\leq S(M).
\]
Consequently \(K_G^{\mathbb C}<1/0.7118598\).

We use the exact parameters and polynomials in the appendix, phase order
\(N=4\), Cauchy radius \(R=103/100\), and scalar truncation parameter
\(d=480\). The retained sum has coefficient enclosures through degree
\(2d+1=961\). Its certified linear lower bound \(L\), finite nonlinear
norm upper bound \(Q\), and complete scalar tail upper bound \(T\) satisfy
\[
 L>0.711860000025,\qquad Q<0.000000000103,\qquad
 T<0.000000000013.                                    \tag{1}
\]
The finite certificate establishing these assertions is described after
the analytic reduction. Every terminating decimal denotes an exact
rational number.

**The coefficient algebra.** Let \(\mathcal W\) be the algebra of absolutely
summable power series on the closed unit disk, with norm equal to the sum of
the absolute values of the coefficients. Its two-variable version has norm
\(\sum_{m,n}|c_{mn}|\). Multiplication is submultiplicative. Substitution of
two series \(P,B\), each with zero constant term and norm at most one, is a
contraction from the two-variable algebra into the one-variable algebra:

\[
 \left\|\sum_{m,n}c_{mn}P^mB^n\right\|_{\mathcal W}
 \leq\sum_{m,n}|c_{mn}|\|P\|_{\mathcal W}^m\|B\|_{\mathcal W}^n
 \leq\sum_{m,n}|c_{mn}|.                               \tag{2}
\]

**Gaussian contraction.** A standard proper complex Gaussian has density
\(\pi^{-1}e^{-|w|^2}\), so \(\mathbb E|W|^2=1\). Let \(W,X\) be independent
such variables and let \(W',X'\) have the same marginal distribution, with
real cross correlations \(a,b\), respectively, and all other cross
correlations zero. For any two \(L^2\) functions \(f,g\), their correlation
has an absolutely summable expansion

\[
 C_{f,g}(a,b)=\mathbb E f(W,X)\overline{g(W',X')},
 \qquad \|C_{f,g}\|_{\mathcal W(\mathbb D^2)}\leq\|f\|_2\|g\|_2. \tag{3}
\]

To prove this, expand \(f,g\) in the orthonormal product Hermite basis of
\((\sqrt2\Re W,\sqrt2\Im W,\sqrt2\Re X,\sqrt2\Im X)\). The generating
identity for real standard Gaussians of correlation \(a\) is

\[
 \mathbb E e^{sY-s^2/2}e^{tY'-t^2/2}=e^{ast}.
\]

It implies that the correlation of equal normalized Hermite degrees \(n\)
is \(a^n\), while unequal degrees are orthogonal. If \(f_\nu,g_\nu\) are
the four-coordinate coefficients, the kernel is therefore

\[
 \sum_{\nu\in\mathbb N^4}f_\nu\overline{g_\nu}
 a^{\nu_1+\nu_2}b^{\nu_3+\nu_4}.
\]

Its absolute coefficient sum is at most
\(\sum_\nu|f_\nu g_\nu|\leq\|f\|_2\|g\|_2\) by Parseval and
Cauchy–Schwarz. Grouping equal bidegrees cannot increase the norm. Hermite
polynomials are complete: a function orthogonal to all of them has zero
Gaussian exponential transform, hence zero Fourier transform after
multiplication by the Gaussian density. The transform is justified by
Cauchy–Schwarz. Finite Hermite approximations prove (3), including degenerate
endpoint correlations, without requiring pointwise regularity of \(f,g\).
Combining (2) and (3) gives

\[
 \|C_{f,g}(P,B)\|_{\mathcal W}\leq\|f\|_2\|g\|_2.   \tag{4}
\]

**Signed complex tensor preprocessing.** Put
\(T_j(u)=u^{\otimes(j+1)}\otimes\overline u^{\otimes j}\), where a basis
fixes the conjugate Hilbert space. For unit vectors,

\[
 \langle T_j(u),T_j(v)\rangle
 =\langle u,v\rangle|\langle u,v\rangle|^{2j},\qquad
 \|T_j(u)\|=1.                                        \tag{5}
\]

If \(P(t)=\sum p_jt^{2j+1}\) has real coefficients and
\(\sum|p_j|\leq1\), take the direct sum of
\(\sqrt{|p_j|}T_j(u)\) on the row side and
\(\operatorname{sgn}(p_j)\sqrt{|p_j|}T_j(v)\) on the column side.
Add a row padding coordinate and a different, orthogonal column padding
coordinate, each with length \(\sqrt{1-\sum|p_j|}\). The resulting vectors
are unit vectors and have cross inner product

\[
 \widetilde P(z)=\sum_jp_j z|z|^{2j}.                  \tag{6}
\]

Apply this construction separately to \(P\) and \(B\). Independent standard
complex Gaussian functionals on the two resulting Hilbert spaces give the
primary coordinates \(W,W'\) and auxiliary coordinates \(X,X'\). Their
cross correlations are exactly \(\widetilde P(z),\widetilde B(z)\), and
all primary–auxiliary covariances vanish. Signs on the column side and the
two different padding spaces are essential parts of this construction.

**Phases and the universal coefficient criterion.** Define

\[
 q(w,x)=\Im(x\overline w),\quad E(w,x)=|w|^2+|x|^2,
\]
\[
 \Psi(w,x)=\varepsilon_0\sin(2\lambda_0q)
 +\varepsilon_1 e^{-\eta E}\sin(2\lambda_1q),
\]
\[
 f_t(w,x)=\operatorname{ph}(w+\theta x)
 e^{-2ikq+it\Psi(w,x)},\qquad g_t(w,x)=f_t(w,-x).      \tag{7}
\]

Here all parameters are real and \(\eta\geq0\). Give \(\operatorname{ph}(0)\)
any unit-modulus value. The linear Gaussian argument \(W+\theta X\) has
variance \(1+\theta^2>0\), so this convention affects no correlation.
These are measurable, modulus-one \(L^2\) functions. They have angular
character one under common multiplication of \(w,x\) by a phase and commute
with complex conjugation, outside that irrelevant null set.

Let \(H_t(s)=C_{f_t,g_t}(P(s),B(s))\). Equation (4) implies absolute
summability. Conjugation symmetry makes its coefficients real, and angular
character at the phase \(-1\) makes it odd. Write

\[
 H_1(s)=\sum_{j\geq0}h_js^{2j+1},\qquad
 \Gamma=h_0-\sum_{j\geq1}|h_j|.                       \tag{8}
\]

For complex \(z=\rho e^{i\phi}\), rotate both row Gaussian coordinates by
\(e^{i\phi}\). Equation (6) shows that the actual rounded correlation is

\[
 e^{i\phi}H_1(\rho)=\sum_{j\geq0}h_jz|z|^{2j}.       \tag{9}
\]

The analytic scalar variable in \(H_1(s)\) and the actual complex
correlation \(z\) thus have different extensions: (9) is a radial
extension, not the substitution \(s=z\) into a holomorphic function.

Choose a maximizing vector configuration for \(M\) and rotate all row
vectors so its objective equals the positive real number \(V(M)\).
The expected scalar objective after (7) is bounded in absolute value by
\(S(M)\). Every term with \(j\geq1\) in (9), by (5), has matrix objective
of absolute value at most \(V(M)\). Absolute summability justifies the
exchange of the finite matrix sum and the coefficient series. Therefore

\[
 S(M)\geq h_0V(M)-\sum_{j\geq1}|h_j|V(M)=\Gamma V(M). \tag{10}
\]

This proves the coefficient criterion uniformly for all finite complex
matrices, with complex phases and complex unit vectors.

**The entire phase remainder.** The map in (4) is bounded and sesquilinear.
The \(L^2\)-valued functions \(t\mapsto f_t,g_t\) have derivatives of every
order because \(\Psi\) is bounded. The product rule gives, for real \(t\),

\[
 \|H_t^{(n)}\|_{\mathcal W}
 \leq\sum_{j=0}^n\binom nj\|\Psi^j\|_2\|\Psi^{n-j}\|_2
 \leq(2\|\Psi\|_\infty)^n.                           \tag{11}
\]

Reflection changes the sign of the phase but not these norms. Since
\(2|q|\leq E\),

\[
 e^{-\eta E}|\sin(2\lambda_1q)|
 \leq\min\{1,|\lambda_1|/(e\eta)\}
 \leq\min\{1,2|\lambda_1|/(5\eta)\}                  \tag{12}
\]

for positive \(\eta\); use one when \(\eta=0\). The last inequality follows
from \(e>5/2\), for example by its first four series terms. In this case
the last bound is \(1/3\), so

\[
 A=|\varepsilon_0|+|\varepsilon_1|/3
 =\frac{86385277374369871}{3000000000000000000}.
\]

Banach-space Taylor's formula with its integral remainder now proves

\[
 \left\|H_1-\sum_{n=0}^4\frac{H_0^{(n)}}{n!}\right\|_{\mathcal W}
 \leq E_{\rm phase}:=\frac{(2A)^5}{120}<0.00000000528. \tag{13}
\]

This is a bound in the complete coefficient norm. It includes the linear
coefficient and every nonlinear scalar degree at every omitted phase order.
It is not an estimate restricted to a finite scalar truncation.

**Exact primitive weights and the two omissions.** For one profile, elementary
algebra gives the uniformly absolutely convergent identity

\[
 e^{i\varepsilon e^{-\eta E}\sin(2\lambda q)}
 =\sum_{n\geq0}\sum_{j=0}^n
 \frac{(\varepsilon/2)^n(-1)^j}{j!(n-j)!}
 e^{-n\eta E}e^{2i(n-2j)\lambda q}.                   \tag{14}
\]

The absolute coefficient sum at order \(n\) is
\(|\varepsilon|^n/n!\). Thus all rearrangements are justified. Multiply the
two profiles separately on the row and column sides and retain precisely
the terms whose combined order is at most four. This gives the Taylor
polynomial in (13), not two independent fourth-order truncations.

A row atom with indices \(n_0,j_0,n_1,j_1\) has

\[
 r=k-(n_0-2j_0)\lambda_0-(n_1-2j_1)\lambda_1,\quad
 \alpha=n_1\eta,
\]
\[
 c=\frac{(\varepsilon_0/2)^{n_0}(\varepsilon_1/2)^{n_1}
 (-1)^{j_0+j_1}}
 {j_0!(n_0-j_0)!j_1!(n_1-j_1)!}.                       \tag{15}
\]

Column atoms have the same formula with \(r,\alpha\) replaced by \(s,\beta\).
Both quadratic phase signs in the correlation are negative: reflection
negates \(q\), and conjugation of the column function negates its exponent
again. Both perturbation signs in (14) are consequently positive.
The primitive is

\[
 F_{r,s;\alpha,\beta}(a,b)=\mathbb E\left[
 \operatorname{ph}(W+\theta X)\overline{\operatorname{ph}(W'-\theta X')}
 e^{-2ir q(W,X)-2is q(W',X')-\alpha E(W,X)-\beta E(W',X')}\right]. \tag{16}
\]

For real \(a,b\), the change of variables
\((w,x,w',x')\mapsto(\overline{w'},-\overline{x'},\overline w,-\overline x)\)
preserves the Gaussian law and the phase product and exchanges the two
frequency–damping pairs. Thus \(F_{r,s;\alpha,\beta}=F_{s,r;\beta,\alpha}\).
Combining unordered pairs in (15) with their exact weights gives 210
nonzero modes for the appendix data. Equivalently, these weights can be
defined without a factorial formula by repeated multiplication of the eight
elementary terms, two signs for each profile on each side, dividing the
order-\(n\) product by \(n\) at each step.

Omit the two modes with pair sets

\[
 \{(k-4\lambda_0,0),(k,0)\},\qquad
 \{(k,0),(k+4\lambda_0,0)\}.                           \tag{17}
\]

Each has weight \(\varepsilon_0^4/192\), including the exchange
multiplicity. The two functions in a primitive have \(L^2\) norms
\((1+2\alpha)^{-1}\) and \((1+2\beta)^{-1}\), because
\(\mathbb E e^{-2\alpha E}=(1+2\alpha)^{-2}\). Equation (4) proves

\[
 \|F_{r,s;\alpha,\beta}(P,B)\|_{\mathcal W}
 \leq\frac1{(1+2\alpha)(1+2\beta)}.
\]

The complete cost of (17) is therefore exactly bounded by

\[
 E_{\rm omit}=\varepsilon_0^4/96<0.000000000026.       \tag{18}
\]

This requires no extension of either omitted mode beyond the unit disk.
Let \(G=\sum_{\ell\in I}w_\ell F_\ell(P,B)\), where \(I\) consists of the
remaining 208 modes. Then
\(\|H_1-G\|_{\mathcal W}\leq E_{\rm phase}+E_{\rm omit}\).

**Evaluation of a damped Gaussian primitive.** Initially take real
\(|a|,|b|<1\). In the order \(W,X,W',X'\), put

\[
 S=\begin{pmatrix}1&0&a&0\\0&1&0&b\\a&0&1&0\\0&b&0&1\end{pmatrix},\quad
 J=\begin{pmatrix}-\alpha&-r&0&0\\r&-\alpha&0&0\\0&0&-\beta&-s\\0&0&s&-\beta\end{pmatrix},
 \quad U=\begin{pmatrix}1&\theta&0&0\\0&0&1&-\theta\end{pmatrix}.
\]

The exponent in (16) is precisely \(z^*Jz\), since
\(-2ir\Im(x\overline w)=-r\overline wx+r\overline xw\).
Define

\[
 D=\det(I-SJ),\quad C=(S^{-1}-J)^{-1},\quad
 UCU^T=\begin{pmatrix}n_1&p\\q_1&n_2\end{pmatrix},
\]
\[
 \begin{pmatrix}N_1&P_0\\Q_0&N_2\end{pmatrix}
 =U\operatorname{adj}(I-SJ)SU^T.
\]

Thus \(n_i=N_i/D, p=P_0/D, q_1=Q_0/D\), and the five capital-letter
quantities are exact polynomials in \(a,b,\theta,r,s,\alpha,\beta\).
Their definitions by determinants and minors specify every primitive
without any simplified rational formula.

For a matrix \(A\) with positive Hermitian part, complex Gaussian integration
gives

\[
 \int_{\mathbb C^n}e^{-z^*Az}\,dz=\frac{\pi^n}{\det A},\qquad
 \int_{\mathbb C^n}z_i\overline z_j e^{-z^*Az}\,dz
 =\frac{\pi^n(A^{-1})_{ij}}{\det A}.                   \tag{19}
\]

Here \(dz\) is real Lebesgue measure on \(\mathbb R^{2n}\). The first identity
follows by unitary diagonalization when \(A\) is positive Hermitian. For
general \(A=H+iK\), with \(H>0\) and \(K\) Hermitian, apply the one-variable
identity theorem to \(H+\zeta K\): both sides are analytic on the connected
strip where \(H+(\Re\zeta)K>0\), which contains the imaginary axis and a real
neighborhood of zero. The second identity follows by differentiating in
\(A_{ji}\); a Gaussian integrable majorant justifies differentiation.

For \(y\ne0\),

\[
 \operatorname{ph}(y)=\frac y{\sqrt\pi}
 \int_0^\infty t^{-1/2}e^{-t|y|^2}\,dt.               \tag{20}
\]

Apply (20) to the two phases in (16). Absolute integration is legitimate:
after taking absolute values, the two Laplace integrals cancel the two
phase denominators, and the remaining damping is at most one. Equation
(19) after the two nonnegative rank-one precision updates gives

\[
 F=\frac p{\pi D}\int_0^\infty\int_0^\infty
 \frac{(uv)^{-1/2}\,du\,dv}
 {\big((1+un_1)(1+vn_2)-uvpq_1\big)^2}.              \tag{21}
\]

At \(a=b=0\), \(p=q_1=0\), \(D>0\), and \(n_1,n_2>0\).
In a neighborhood of this point, expand the denominator of (21).
The term of order \(j\) contains

\[
 (j+1)(pq_1)^j
 \prod_{m=1}^2\int_0^\infty
 \frac{u^{j-1/2}\,du}{(1+n_m u)^{j+2}}
 =(j+1)(pq_1)^j(n_1n_2)^{-j-1/2}
 B(j+1/2,3/2)^2.
\]

The geometric expansion has an integrable uniform majorant there. Using
\(B(j+1/2,3/2)=\frac\pi2(1/2)_j/(j+1)!\) proves

\[
 F(z)=\frac{\pi P_0(z)}{4D(z)\sqrt{N_1(z)N_2(z)}}
 \Phi\!\left(\frac{P_0(z)Q_0(z)}{N_1(z)N_2(z)}\right),\quad
 \Phi(w)=\sum_{j\geq0}\frac{(1/2)_j^2}{j!(j+1)!}w^j. \tag{22}
\]

Here \(a=P(z),b=B(z)\), and the square root is fixed by its positive value
at zero. The leading factor is \(\pi/4\), not a real-sign rounding factor.
For example, with \(\theta=r=s=\alpha=\beta=0\), (22) reduces to
\((\pi a/4)\Phi(a^2)\). Its derivative at zero is \(\pi/4\) and its limit
at \(a=1\) is one, since \(\Phi(1)=4/\pi\). These are the standard complex
phase normalizations.

**Branch selection on the open unit disk.** For complex \(z\), the symmetric
matrix \(S(z)\) is a formal analytic matrix, not a probability covariance
matrix. Nonetheless it has positive Hermitian part when \(|z|<1\):

\[
 \|S-I\|=\max\{|P(z)|,|B(z)|\}\leq |z|<1.
\]

For any invertible \(A\),
\(\Re(A^{-1})=A^{-*}(\Re A)A^{-1}\). Thus \(S^{-1}\) is strictly
accretive. The Hermitian part of \(J\) is
\(-\operatorname{diag}(\alpha,\alpha,\beta,\beta)\leq0\), so
\(S^{-1}-J\) and \(C=(S^{-1}-J)^{-1}\) are strictly accretive. In particular
\(D=\det(S)\det(S^{-1}-J)\ne0\).

The real matrix \(U\) has full row rank for every real \(\theta\).
Consequently the compression \(V=UCU^T=UCU^*\) is strictly accretive:
for every nonzero \(v\in\mathbb C^2\),

\[
 \Re(v^*Vv)=\Re((U^*v)^*C(U^*v))>0.
\]

Its two diagonal entries have positive real parts, hence \(N_1,N_2\ne0\).
For every \(0\leq t\leq1\),

\[
 V_t=\begin{pmatrix}n_1&\sqrt t\,p\\\sqrt t\,q_1&n_2\end{pmatrix}
 =\sqrt t\,V+(1-\sqrt t)\operatorname{diag}(n_1,n_2)
\]

is strictly accretive. Therefore
\(n_1n_2-tpq_1\ne0\). If

\[
 Z=\frac{P_0Q_0}{N_1N_2}=\frac{pq_1}{n_1n_2}
\]

were a real number at least one, choosing \(t=1/Z\) would contradict this
nonvanishing. Thus \(Z\notin[1,\infty)\) everywhere on the open unit disk.
This proves slit avoidance there by compression and off-diagonal
interpolation, without assuming that \(V\) is Hermitian or symmetric.

**Extension to the larger Cauchy circle.** The preceding proof uses
\(|z|<1\); it gives no automatic conclusion for \(1\leq|z|\leq R\).
For each retained primitive, the certificate supplies a finite cover of this entire closed
annulus by sets

\[
 \{\rho e^{i\pi t}:\rho\in[r_-,r_+],\ t\in[t_-,t_+]\},
\]

with exact rational endpoints. On each set, interval arithmetic must prove
\(D,N_1,N_2\ne0\) and a strictly positive lower bound for

\[
 \delta(Z)=\max\left\{\min(1,1-\Re Z),\frac{|\Im Z|}{|Z|}\right\}, \tag{23}
\]

where the second entry is zero at \(Z=0\). For any real \(x,y\),

\[
 |1-t(x+iy)|^2=
 \left(|x+iy|t-\frac{x}{|x+iy|}\right)^2+
 \frac{y^2}{|x+iy|^2}
\]

when \(x+iy\ne0\); also \(|1-tZ|\geq\min(1,1-\Re Z)\) for
\(0\leq t\leq1\). Thus (23) is a simultaneous lower bound for
\(|1-tZ|\), and positivity excludes the slit. An interval evaluation may
replace the second entry by zero if division cannot be certified; this only
weakens the test.

Together with the open-unit-disk argument, these annulus conditions prove
nonvanishing and slit avoidance on the whole closed disk of radius \(R\).
Compactness and continuity extend them to a slightly larger disk.
There is an analytic logarithm of \(N_1N_2\) on this simply connected disk,
and hence a unique analytic square root with the prescribed positive value
at zero. Its value need not be the pointwise principal square root at points
on the larger circle. Formula (22) uses its analytic continuation instead.

The Euler integral is

\[
 \Phi(w)=\frac2\pi\int_0^1
 t^{-1/2}(1-t)^{1/2}(1-tw)^{-1/2}\,dt,
 \qquad w\notin[1,\infty).                            \tag{24}
\]

The integrand uses the principal root; for these \(w\), its argument never
lies on the negative real axis or at zero. The positive integration weight
has total mass one. Equation (24) follows by expanding near zero and
integrating beta coefficients, then by analytic continuation; it is also
the specialization of the [Euler hypergeometric integral](https://dlmf.nist.gov/15.6.E1).
If an arc has a certified lower bound \(\delta>0\) for every
\(|1-tZ|\), then \(|\Phi(Z)|\leq\delta^{-1/2}\) on that entire arc.

There is also a complete finite interval cover of the circle \(|z|=R\), with the
same nonvanishing and positive-distance checks. Its local upper bounds

\[
 \frac{\pi |P_0|}{4|D|\sqrt{|N_1||N_2|}\sqrt\delta}   \tag{25}
\]

give an upper bound \(M_\ell\) for primitive \(\ell\) on the whole circle.
The expression is independent of the sign of the continued square root.
The continued expression agrees with the Gaussian primitive near zero by
(21), and with its coefficient series on the unit disk by the identity
theorem. No nonholomorphic Gaussian expectation at a complex mean is used.

**Exact finite coefficients and every scalar tail.** For one primitive write
\(a=zA(y), b=zB_*(y), y=z^2\), where \(A=P(z)/z\) and \(B_*=B(z)/z\).
The determinant and marginal numerators are even polynomials in \(z\);
\(P_0=zL(y), Q_0=zQ_*(y)\). Put

\[
 T(y)=N_1(y)N_2(y),\qquad T_0=T(0)>0,\qquad
 Z(y)=yL(y)Q_*(y)/T(y).
\]

Equation (22) becomes

\[
 \frac{F(z)}z=\frac\pi{4\sqrt{T_0}}
 \frac{L(y)}{D(y)}\left(\frac{T(y)}{T_0}\right)^{-1/2}
 \sum_{j\geq0}\frac{\binom{2j}{j}^2}{16^j(j+1)}Z(y)^j. \tag{26}
\]

Everything following the factor \(\pi/(4\sqrt{T_0})\) is a formal series
over the rationals. The reciprocal coefficients are determined by product
equal to one, and the square-root coefficients by its square equal to
\(T/T_0\), with constant coefficient one. In particular these operations
choose exactly the branch used in (22). The hypergeometric coefficients
satisfy \(c_0=1\) and

\[
 \frac{c_j}{c_{j-1}}=\frac{(2j-1)^2}{4j(j+1)}.
\]

Since \(Z(0)=0\), powers with \(j>d\) contribute no coefficient through
\(y^d\). Consequently (26) through degree \(d\), computed with exact
rationals and rigorous enclosures for the constant factor, or with outward
formal interval arithmetic, encloses the exact coefficients through
\(z^{2d+1}\). This is an algebraic truncation, not an approximation to the
hypergeometric function at a nonzero argument.

The primitive is odd. Cauchy's coefficient bound on the verified circle
therefore gives its entire remaining coefficient norm:

\[
 \sum_{j>d}|[z^{2j+1}]F_\ell(z)|
 \leq M_\ell\sum_{j>d}R^{-2j-1}
 =M_\ell\frac{R^{-2d-3}}{1-R^{-2}}.                   \tag{27}
\]

For retained coefficient enclosures \(I_{\ell j}\), form interval sums
\(I_j=\sum_{\ell\in I}w_\ell I_{\ell j}\). Define

\[
 L=\inf I_0,\quad Q=\sum_{j=1}^d\sup_{x\in I_j}|x|,
 \quad T=\sum_{\ell\in I}|w_\ell|M_\ell
 \frac{R^{-2d-3}}{1-R^{-2}}.                           \tag{28}
\]

For any real coefficient error \(e\),

\[
 (h_0+e_0)-\sum_{j\geq1}|h_j+e_j|
 \geq h_0-\sum_{j\geq1}|h_j|-\sum_{j\geq0}|e_j|.
\]

Thus the complete error norm is subtracted once, including its linear
part. Equations (13), (18), and (27) imply

\[
 \Gamma\geq L-Q-T-E_{\rm phase}-E_{\rm omit}.
\]

By (1), this is strictly greater than

\[
 0.711860000025-0.000000000103-0.000000000013
 -0.00000000528-0.000000000026
 =0.711859994603>0.71185999.                            \tag{29}
\]

Equation (10) proves the theorem, including its strict comparison with
\(1/0.7118598\). All infinite phase orders, both omitted modes, and all
infinite scalar tails have been accounted for in (29).

**Finite certificate and its verification.** The certificate contains 208
retained modes, each with 481 exact dyadic coefficient balls. A ball with
midpoint pair \((m,e)\) and radius pair \((r,f)\) means
\([m2^e-r2^f,m2^e+r2^f]\). The parameters, mode weights, and polynomial
definitions are exact rationals. Each mode also has a finite family of
rational polar rectangles and strictly positive dyadic lower bounds for
\(|D|,|N_1|,|N_2|\), and the Euler distance on every rectangle.

The initial angular intervals
\([2j/8192,2(j+1)/8192]\), \(0\leq j<8192\), cover \([0,2]\),
including the common endpoint representing angles zero and \(2\pi\).
Their radial roots are \([1,R]\) for the annulus and \([R,R]\)
for the circle. The complete certificate consists of 1,703,936 boundary
leaves and 1,832,072 annulus leaves. These counts describe the finite
certificate; coverage follows from the full split trees, not from counting
alone.

For each root, the verifier reconstructs every midpoint split with rational
arithmetic. Each internal vertex must have both children and one consistent
split coordinate. A leaf may have no descendant, and duplicate leaves are
forbidden. All stored endpoints must equal the reconstructed endpoints.
Every root, in both the circle and annulus families, must be present.
Induction down each finite binary tree proves complete coverage, including
every endpoint.

A polar rectangle is enclosed by
\[
 z\in [r_-,r_+]\exp\!\big(i\pi[t_-,t_+]\big).
\]
The midpoint and half-width of each interval are exact rationals before
outward interval conversion. Polynomial evaluation encloses its complete
image. For the interval value \(\mathbf Z\), the arithmetic uses the lower
bounds
\[
 \min\{1,\operatorname{lower}(1-\Re\mathbf Z)\},
 \qquad
 \frac{\operatorname{lower}|\Im\mathbf Z|}
      {\operatorname{upper}|\mathbf Z|}.
\]
The second is replaced by zero if its denominator is zero. If the interval
contains \(Z=0\), its imaginary-modulus lower bound cannot be positive,
so this rule still gives a valid simultaneous bound. The maximum of these
two lower bounds is rounded downward. Strict positivity is checked only
after that rounding. Nonfinite evaluations or unresolved leaves cannot be
accepted.

All accepted local lower endpoints and upper endpoints are further rounded
in the outward direction to dyadics. The signed-integer floor operation is
used for lower endpoints and ceiling for upper endpoints, so negative
values are handled correctly. Coefficient balls instead retain their full
binary midpoints and radii. Thus serialization does not require a decimal
round trip. JSON integer tokens are parsed as arbitrary-precision integers,
without conversion to binary floating-point numbers.

For these certificates, uniform conservative consequences are
\[
 |D|>\frac3{100},\qquad |N_1|>\frac1{10},\qquad
 |N_2|>\frac7{50},\qquad
 |1-tZ|>\frac1{200000}\quad(0\leq t\leq1)              \tag{30}
\]
on every retained annulus and circle cover. Per-mode lower bounds,
which are stronger than these common bounds, are retained for the
scalar-tail calculation.

The circle upper bounds in (25) can be checked without any numerical
square root. If \(d,n_1,n_2,\Delta>0\) are the local lower bounds,
\(p\geq0\) is the numerator upper bound, and \(\pi_+\geq\pi\), a
nonnegative rational \(m\) is certified as a local modulus upper bound if
the sufficient inequality
\[
 16m^2d^2n_1n_2\Delta\ \geq\ \pi_+^2p^2               \tag{31}
\]
holds. This is a rational comparison. Taking a maximum over all circle
leaves gives the mode's global boundary bound.

There is a second, entirely rational calculation of the boundary constants
which uses no per-arc numerator or modulus field. Define
\[
 A_R=\sum_j|p_j|R^{2j+1},\qquad
 B_R=\sum_j|b_j|R^{2j+1}.
\]
For the exact polynomial
\(P_{0,\ell}(a,b)=\sum_{u,v}c_{\ell uv}a^ub^v\), let
\[
 U_\ell=\sum_{u,v}|c_{\ell uv}|A_R^uB_R^v.
\]
This bounds its modulus on the whole circle. Let
\(d_\ell,n_{\ell1},n_{\ell2},\Delta_\ell\) be the certified positive
per-mode minima, and put
\[
 s_\ell=
 \frac{\left\lfloor10^{35}\sqrt{n_{\ell1}n_{\ell2}\Delta_\ell}\right\rfloor}
      {10^{35}}>0,\qquad
 M_\ell=\frac{(22/7)U_\ell}{4d_\ell s_\ell}.             \tag{32}
\]
The integer in \(s_\ell\) is obtained by integer square-root isolation;
\(s_\ell^2\leq n_{\ell1}n_{\ell2}\Delta_\ell\) is also checked exactly.
The elementary inequality \(\pi<22/7\) follows, for example, by expanding
the positive integral
\[
 0<\int_0^1\frac{x^4(1-x)^4}{1+x^2}\,dx=\frac{22}{7}-\pi.
\]
Thus (32) supplies rational boundary bounds by (25).
Their complete weighted tail is
\[
 \frac{R^{-963}}{1-R^{-2}}\sum_{\ell\in I}|w_\ell|M_\ell
 <1.264554344531686\times10^{-11}<13\times10^{-12}.      \tag{33}
\]
The first omitted degree is 963; every higher odd degree is included.

Finally, exact rational addition of the 481 weighted coefficient intervals
per mode gives
\[
 L>0.711860000025,\qquad
 Q<1.021854640\times10^{-10}<1.03\times10^{-10}.
\]
Together with (33), these are (1). Formula (26), exact determinant/minor
polynomials, and outward formal arithmetic establish the individual
coefficient enclosures; geometric coverage and rational aggregation
establish their complete use. The arithmetic is specified in
[replay_independent.py](replay_independent.py), the exact split-tree and
aggregation rules in [validate_replay.py](validate_replay.py), and the
alternative rational bound (32) in
[check_rational_tails.py](check_rational_tails.py).
The local transcendental and formal-series enclosures use Arb's outward
interval operations. The exact tree and aggregation checks verify the
subsequent rational implications; they do not independently implement
those Arb operations.

**Appendix: exact input data.** The parameters are

\[
\begin{aligned}
\theta&=0.32180438477341805\\
k&=0.004478959445491443\\
\varepsilon_0&=0.007029952343109007\\
\lambda_0&=1.0583917494573267\\
\varepsilon_1&=0.06529542034504285\\
\lambda_1&=1/4,\qquad\eta=3/10.
\end{aligned}
\]

Write \(P(t)=\sum_{j=0}^{49}p_jt^{2j+1}\). Its exact coefficients are

\[
\begin{aligned}
p_{0}&=0.99938755437877901944\\
p_{1}&=0.00000856028176808403\\
p_{2}&=0.00012224125349985935\\
p_{3}&=-0.00043143687566471554\\
p_{4}&=-0.00002366553779359520\\
p_{5}&=-0.00000029862574899892\\
p_{6}&=0.00000442933793458406\\
p_{7}&=0.00000714867787922735\\
p_{8}&=0.00000609738892318563\\
p_{9}&=0.00000254980511744320\\
p_{10}&=0.00000172442641254250\\
p_{11}&=0.00000188680068544758\\
p_{12}&=0.00000053919240458309\\
p_{13}&=0.00000037696428991676\\
p_{14}&=0.00000058910851684379\\
p_{15}&=0.00000010812366533416\\
p_{16}&=0.00000009428392941230\\
p_{17}&=0.00000017618442612030\\
p_{18}&=0.00000003948027132926\\
p_{19}&=0.00000002787879670676\\
p_{20}&=0.00000004807324311830\\
p_{21}&=0.00000002832824800300\\
p_{22}&=0.00000000598217111181\\
p_{23}&=0.00000001089263648086\\
p_{24}&=0.00000001873434382704\\
p_{25}&=-0.00000000071835196755\\
p_{26}&=0.00000000227258101415\\
p_{27}&=0.00000001015470841133\\
p_{28}&=-0.00000000121170758900\\
p_{29}&=0.00000000061811875885\\
p_{30}&=0.00000000443956314480\\
p_{31}&=-0.00000000037297142386\\
p_{32}&=0.00000000016055293961\\
p_{33}&=0.00000000158871801423\\
p_{34}&=0.00000000020745930099\\
p_{35}&=-0.00000000004401233185\\
p_{36}&=0.00000000047486524671\\
p_{37}&=0.00000000033497359846\\
p_{38}&=-0.00000000010830721766\\
p_{39}&=0.00000000012157351800\\
p_{40}&=0.00000000024624263248\\
p_{41}&=-0.00000000008623167083\\
p_{42}&=0.00000000002762569817\\
p_{43}&=0.00000000013127267652\\
p_{44}&=-0.00000000004272848749\\
p_{45}&=0.00000000000330337134\\
p_{46}&=0.00000000005595196583\\
p_{47}&=-0.00000000001204584236\\
p_{48}&=-0.00000000000380444644\\
p_{49}&=0.00000000001983139326
\end{aligned}
\]

The auxiliary polynomial is

\[
B(t)=-\frac{3349062818954309}{100000000000000000}t
+\frac{96650937181045691}{100000000000000000}t^3.
\]

Their exact absolute coefficient sums are

\[
\|P\|_{\mathcal W}=\frac{124999960009331519}{125000000000000000}
=0.999999680074652152<1,\qquad\|B\|_{\mathcal W}=1.
\]
