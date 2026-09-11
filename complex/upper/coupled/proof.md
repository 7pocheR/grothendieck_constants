# A coupled Gaussian upper bound for the complex Grothendieck constant

For every finite complex matrix, the construction and the specified complete
outward-arithmetic calculation give

\[
K_G^{\mathbb C}\le \frac{200000000000}{142379302169}
<1.404698554869.
\]

For \(A=(a_{ij})\in\mathbb C^{m\times n}\), use the normalization
\[
S_{\mathbb C}(A)=\sup_{|s_i|=|t_j|=1}
 \left|\sum_{i,j}a_{ij}s_i\overline{t_j}\right|,
\qquad
V_{\mathbb C}(A)=\sup_{\|u_i\|=\|v_j\|=1}
 \left|\sum_{i,j}a_{ij}\langle u_i,v_j\rangle\right|.
\]
The vector supremum allows arbitrary finite complex Hilbert-space dimension.
The least constant in \(V_{\mathbb C}(A)\le K_G^{\mathbb C}S_{\mathbb C}(A)\)
over all such matrices is the classical complex Grothendieck constant. The
unit-phase scalar normalization is equivalent to allowing scalar modulus at
most one: each point of the disk is an average of unit phases, and the
bilinear form is affine in each scalar separately. A common scalar phase
also makes the modulus and real-part normalizations equivalent.

The coefficient estimate below uses every term in the finite phase expansion,
including an entire coefficient-norm cost for every omitted mode. It asserts
an upper bound, not an exact value of the constant. The rational endpoint is
obtained by integer comparison with the complete rational lower bound for the
coefficient gap in Section 7.

All Gaussian variables below are standard proper complex Gaussians: one coordinate has density \(\pi^{-1}e^{-|z|^2}\), and \(\mathbb E|Z|^2=1\). Inner products are linear in their first argument. Matrices in the Grothendieck inequality have arbitrary complex entries. Vector dimension is unrestricted.

**1. Exact parameters and admissible correlations.**

The rational parameters and the three polynomials are specified completely in the appendix. In particular,
\[
 h=\frac{15869581262852117}{50000000000000000},\quad
 k=\frac{465518641186931}{250000000000000000},\quad
 m=-\frac1{50},\quad \epsilon=\frac1{100},\quad\lambda=\frac14.
\]
The damping coefficients are
\[
 d_w=\frac{9209367321696771}{10^{16}},\qquad
 d_x=\frac{2171701844586119}{2000000000000000}.
\]
They are positive and \(d_wd_x\ge g^2\) for \(g=999/1000\). There are five real profiles \((e_j,\ell_j,u_j)\), all with \(u_j\ge0\). Set
\[
 q_X=\operatorname{Im}(X\overline W),\quad
 q_Y=\operatorname{Im}(Y\overline W),\quad
 R=d_w|W|^2+d_x|X|^2,
 \qquad \Psi_0(W,X)=\sum_{j=1}^5 e_j e^{-u_jR}\sin(2\ell_j q_X).
\]
Define measurable functions, of modulus one almost everywhere, by
\[
 f(W,X,Y)=\operatorname{ph}(W+hX+mY)
     e^{-2ikq_X}e^{i\Psi_0(W,X)+i\epsilon\sin(2\lambda q_Y)},
 \qquad g_0(W,X,Y)=f(W,-X,-Y).
\]
The value assigned to \(\operatorname{ph}(0)\) is immaterial: every individual linear argument has variance \(1+h^2+m^2>0\). Both functions have common rotational charge one, and both commute with complex conjugation.

For a real odd polynomial \(Q(x)=\sum_{j\ge0}q_{2j+1}x^{2j+1}\), write
\[
 \widetilde Q(t)=\sum_{j\ge0}q_{2j+1}t|t|^{2j},\qquad
 \|Q\|_1=\sum_j|q_{2j+1}|.
\]
Exact rational addition gives
\[
 \|P\|_1=\frac{1999999015980406524216029}{2000000000000000000000000}<1,
 \qquad \|B\|_1=\|C\|_1=1,\qquad C(x)=x^7.
\]

**Lemma 1 (dimension-independent realization).** If \(\|Q\|_1\le1\), any finite collection of unit vectors \(u_i,v_j\) admits unit vector transforms with cross inner products \(\widetilde Q(\langle u_i,v_j\rangle)\).

*Proof.* The tensor vector
\(T_j(u)=u^{\otimes(j+1)}\otimes\overline u^{\otimes j}\) has norm one and cross inner product \(t^{j+1}\overline t^j\). In mutually orthogonal summands put \(\sqrt{|q_{2j+1}|}T_j(u)\) on the row side and \(\operatorname{sgn}(q_{2j+1})\sqrt{|q_{2j+1}|}T_j(v)\) on the column side. Their squared norms equal \(\|Q\|_1\). Complete norms using two additional mutually orthogonal directions, one for all rows and one for all columns. This does not change any cross inner product. All the polynomials here are finite. Apply this construction independently for \(P,B,C\), using three independent Gaussian linear functionals. For each row/column pair the coordinate correlations are precisely \(\widetilde P(t),\widetilde B(t),\widetilde C(t)\), and each individual triple consists of independent standard proper Gaussians. Thus this construction works simultaneously for every finite matrix, with no dependence of its estimates on vector dimension. ∎

Let \(K(t)\) denote the covariance \(\mathbb E f(W,X,Y)\overline{g_0(W',X',Y')}\) under these three correlations. The primed triple is independent within itself; cross correlations between different coordinate types vanish. Cross correlations of the same type are those just specified.

**2. The entire coefficient norm and the bounded-phase error.**

**Lemma 2 (Gaussian coefficient contraction).** For \(u,v\in L^2\) of three independent standard proper complex Gaussian variables, their covariance kernel, after substitution of any three admissible correlation maps, has absolute coefficient sum at most \(\|u\|_2\|v\|_2\). If both functions have common rotational charge one, its scalar form is
\[
 \sum_{n\ge0}c_n t|t|^{2n},\qquad \sum_n|c_n|\le\|u\|_2\|v\|_2.
\]

*Proof.* Define complex Hermite polynomials by the generating function
\[
 e^{sz+t\overline z-st}=\sum_{p,q\ge0}\frac{H_{pq}(z,\overline z)}{p!q!}s^pt^q.
\]
The Gaussian exponential identity proves that \(H_{pq}/\sqrt{p!q!}\) are orthonormal. Their span is dense: a function orthogonal to all polynomials has Gaussian-weighted Fourier transform with all derivatives zero at zero; Cauchy–Schwarz with Gaussian exponential moments makes this transform entire, so uniqueness of the Fourier transform gives the zero function. Taking products gives a complete orthonormal basis \(h_{\alpha\beta}\) in three coordinates. The same exponential identity, for coordinate correlations \(\rho_j=\mathbb E Z_j\overline{Z'_j}\), gives
\[
 \mathbb E h_{\alpha\beta}(Z)\overline{h_{\mu\nu}(Z')}
 =\mathbf 1_{(\alpha,\beta)=(\mu,\nu)}\rho^\alpha\overline\rho^\beta.
\]
Consequently the covariance is \(\sum u_{\alpha\beta}\overline{v_{\alpha\beta}}\rho^\alpha\overline\rho^\beta\). The sum of absolute weights is bounded by Cauchy–Schwarz and Parseval. This includes the boundary \(|\rho_j|=1\): finite Hermite expansions converge in the covariance by the same uniform contraction estimate, including for singular joint Gaussian distributions.

The algebra of formal series in \(t,\overline t\) with absolutely summable coefficients has a submultiplicative norm. Each substituted map and its conjugate has norm at most one. Thus substitution and collection of equal monomials cannot increase this bound. Charge one restricts the nonzero basis coefficients to \(|\alpha|-|\beta|=1\), so all surviving monomials have the form \(t^{n+1}\overline t^n\). This also proves uniform convergence on the closed disk. ∎

The functions in this construction commute with conjugation, so their covariance on real correlations is real. The scalar coefficients of \(K\), and of each real Fourier–Laplace primitive used below, are therefore real.

Since \(R\ge2g|q_X|\) and \(|\sin v|\le\min(1,|v|)\),
\[
 |\Psi_0|\le M:=\sum_j|e_j|
 \begin{cases}1,&u_j=0,\\
 \min(1,2|\ell_j|/(5u_jg)),&u_j>0.
 \end{cases}
\]
Indeed, \(\sup_{s\ge0}s e^{-2u_jgs}=1/(2u_jg e)\) and \(1/e<2/5\). The latter follows already from \(e>1+1+1/2+1/6>5/2\). For the appendix data, rational arithmetic gives
\[
 M=\frac{19007117794841543023}{399600000000000000000}.
\]

**Lemma 3 (complete rectangular phase remainder).** Leave \(m\) fixed exactly. Retain total row-plus-column Taylor order at most \(L\) in \(\Psi_0\), and independently total row-plus-column order at most \(N\) in the new bounded phase. If \(K_{L,N}\) denotes this finite expansion after substitution, then
\[
 \|K-K_{L,N}\|_1\le
 E_{L,N}:=\frac{(2M)^{L+1}}{(L+1)!}
 +\left(\sum_{j=0}^{L}\frac{(2M)^j}{j!}\right)
       \frac{(2|\epsilon|)^{N+1}}{(N+1)!}. \tag{1}
\]

*Proof.* Multiply the original and new bounded phases by real parameters \(s,t\). Reflection and conjugation make the two copies' bounded phases enter with the same sign. A mixed derivative of orders \(a,b\) is a sum of products of row and column phase powers. The binomial theorem and Lemma 2 bound its entire coefficient norm by \((2M)^a(2|\epsilon|)^b\), uniformly for real \(s,t\), because the unexpanded phase factors remain unitary. Bounded multiplication is differentiable to all orders in Gaussian \(L^2\), and Lemma 2 is a continuous real-bilinear covariance map into the absolute-coefficient Banach space. The integral Taylor formula is therefore valid in that Banach space. First truncate \(s\) at \(t=1\); then truncate \(t\) in each retained \(s\)-derivative at zero. Their respective errors give (1). No differentiability of \(\operatorname{ph}(W+hX+mY)\) in \(m\) is used. ∎

Here \(L=N=4\). The exact error is
\[
 E=\frac{2481861045053572732427524988679905533322860263730665113527571095424564637381739509124313249751343}
 {38208383616191961600000000000000000000000000000000000000000000000000000000000000000000000000000000000000}.
\]
This is an error in the full infinite coefficient norm, not a pointwise error subsequently converted using sampled coefficients. In particular, no fourth-order parity cancellation in \(\epsilon\) is assumed at fixed nonzero \(m\).

**3. Exact mode accounting and the six-coordinate primitive.**

For each original profile, the two generators of \(i\Psi_0\) have weight \(\sigma e_j/2\), frequency shift \(\sigma\ell_j\), and damping \(u_j\), for \(\sigma\in\{-1,1\}\). For the new phase the generators have weight \(\sigma\epsilon/2\), frequency shift \(\sigma\lambda\), and zero damping. Power \(n\) includes its factor \(1/n!\). Thus an original side frequency is \(r=k-\sum\sigma\ell_j\), while a new side frequency is \(\tau=-\sum\sigma\lambda\). Both phases are expanded in total row-plus-column order, as required by (1).

Specify a mode by the ordered side triples \((r,u,\tau),(s,v,\sigma)\). Its Gaussian expectation has integrand
\[
 \operatorname{ph}(W+hX+mY)\overline{\operatorname{ph}(W'-hX'-mY')}
 e^{-uR-vR'}e^{-2i(rq_X+sq'_X+\tau q_Y+\sigma q'_Y)}. \tag{2}
\]
In the coordinate order \(W,X,W',X',Y,Y'\), put
\[
 S=\begin{pmatrix}1&0&a&0&0&0\\0&1&0&b&0&0\\a&0&1&0&0&0\\0&b&0&1&0&0\\0&0&0&0&1&c\\0&0&0&0&c&1\end{pmatrix},
\quad
 J=\begin{pmatrix}
 -ud_w&-r&0&0&-\tau&0\\r&-ud_x&0&0&0&0\\
 0&0&-vd_w&-s&0&-\sigma\\0&0&s&-vd_x&0&0\\
 \tau&0&0&0&0&0\\0&0&\sigma&0&0&0
 \end{pmatrix},
\quad
 U=\begin{pmatrix}1&h&0&0&m&0\\0&0&1&-h&0&-m\end{pmatrix}.
\]
For real \(a,b,c\in(-1,1)\), \(S\) is the actual covariance. Away from the real axis these displayed matrices describe analytic continuation; they are not asserted to be probability covariances.

Define five real rational-coefficient polynomials in \(a,b,c\) by
\[
 D=\det(I-SJ),\qquad
 \begin{pmatrix}N_1&P_0\\Q_0&N_2\end{pmatrix}
 =U\operatorname{adj}(I-SJ)SU^T. \tag{3}
\]
These definitions require neither division by \(\det S\) nor symbolic cancellation of removable singularities. In particular, they remain valid when \(c=\pm1\). Each numerator is also the negative determinant of the bordered matrix with blocks \(I-SJ, SU_j^T; U_i,0\). The supplied source uses these bordered determinants, with exact rational polynomial arithmetic.

**Lemma 4 (analytic Gaussian primitive and normalization).** The value of (2) is the restriction to real correlations of the analytic function
\[
 F(a,b,c)=\frac{\pi P_0}{4D\sqrt{N_1}\sqrt{N_2}}
 {}_2F_1\!\left(\frac12,\frac12;2;\frac{P_0Q_0}{N_1N_2}\right), \tag{4}
\]
with the square roots continued from their positive values at the origin. It is enough to continue where \(D,N_1,N_2\) do not vanish and \(N_1N_2-tP_0Q_0\ne0\) for \(0\le t\le1\).

*Proof.* For nonsingular accretive \(S\), the tilted Gaussian precision is \(A=S^{-1}-J\), also strictly accretive, because the Hermitian part of \(J\) is nonpositive. Gaussian integration gives mass \(D^{-1}\) and compressed covariance
\[
 N=U(I-SJ)^{-1}SU^T=\begin{pmatrix}n_1&p\\q&n_2\end{pmatrix}.
\]
These are identities for integrals against an integrable complex Gaussian density; they do not require \(J\) to be Hermitian. The identity \(\int e^{-z^*Az}\,dz=\pi^d/\det A\), and its differentiated moment identities, follow first for positive Hermitian \(A\) and then by holomorphic continuation on the convex domain of strictly accretive matrices. The same argument proves the compression formula for any fixed bounded function of the two linear arguments.

For a positive Hermitian two-coordinate covariance, the Laplace identity for \(1/(|\xi||\eta|)\) gives the phase covariance
\[
 \frac{p}{\pi}\int_0^\infty\!\int_0^\infty
 \frac{ds\,dt}{\sqrt{st}\big[(1+s n_1)(1+t n_2)-stpq\big]^2}.
\]
Absolute integrability follows by the same Laplace identity applied to the absolute value of the numerator. Substitution \(u=n_1s/(1+n_1s)\), \(v=n_2t/(1+n_2t)\), followed by the expansion of \((1-uvpq/(n_1n_2))^{-2}\), gives
\[
 \frac{\pi p}{4\sqrt{n_1}\sqrt{n_2}}
 \sum_{j\ge0}\frac{(1/2)_j^2}{j!(j+1)!}
               \left(\frac{pq}{n_1n_2}\right)^j.
\]
Here one uses \(B(j+1/2,3/2)=(\pi/2)(1/2)_j/(j+1)!\). This establishes the factor \(\pi/4\).

For an arbitrary strictly accretive \(N\), \(\operatorname{Re}n_i>0\), and the matrix with diagonals \(n_i\) and off-diagonals \(\sqrt t\,p,\sqrt t\,q\) is a convex combination of \(N\) and its diagonal part. It is strictly accretive. Thus \(n_1n_2-tpq\ne0\), and \(pq/(n_1n_2)\) avoids \([1,\infty)\). Both the integral and the formula are holomorphic on this connected domain. Equality on positive Hermitian matrices implies equality there: the Hermitian matrices form a maximal totally real subspace, and the elementary identity theorem applied successively to its four real coordinates supplies the local extension, then connectedness supplies the rest.

Finally \(n_i=N_i/D\), \(p=P_0/D\), \(q=Q_0/D\), and the mass is \(D^{-1}\). At the origin all matrices are real. The homotopy \(I-tJ\) is invertible for \(0\le t\le1\), so its real determinant starts positive and stays positive. The compressed marginal variances are real and strictly positive by accretivity. Hence \(D,N_1,N_2\) are positive there. Cancellation of their analytic square roots gives exactly (4), including its remaining factor \(D^{-1}\). The resulting identity continues uniquely from the origin. ∎

Simultaneously negating \(a,b,c\) negates \(P_0,Q_0\) and preserves \(D,N_1,N_2\), by conjugating the matrices with the diagonal matrix that changes signs on the three column coordinates. Therefore \(F(P(z),B(z),C(z))\) is an odd scalar holomorphic function.

Full exchange of the two side triples preserves \(D,P_0,Q_0\) and exchanges \(N_1,N_2\). To verify this, let \(T\) swap the two coordinate triples and let \(H\) change signs on both \(X\)-coordinates and both \(Y\)-coordinates. Then \(HSH=S\), \(HJH=J^T\), and \(UTH\) is \(U\) with its two rows exchanged. It follows that the compressed covariance under side exchange is the transpose covariance with both rows and columns exchanged, which has entries \(n_2,p;q,n_1\). The asserted polynomial identities follow on an open set and hence identically.

The original unordered expansion has 2962 terms. The new ordered expansion has 41 terms, corresponding to the integer frequency pairs \((i,j)\) with \(|i|+|j|\le4\); their exact weights include the contributions of all compatible orders. For original unequal side data, all 41 ordered pairs remain necessary. For each of the 54 original equal-side terms, full side exchange reduces these to 23 pairs. Thus the complete canonical expansion contains
\[
 (2962-54)\,41+54\,23=120470
\]
nonzero rationally weighted modes. The counts for the original profiles follow by the finite recurrence: start with weight one at shift/damping \((0,0)\), and multiply by each generator with division by the new order at each step. Then multiply the two sides subject to the total-order condition and collect equal keys. This recurrence and a separately written original-mode expansion are supplied as exact executable checks.

Combining the original sides before the new expansion is legitimate only because the *complete ordered new expansion* is invariant under exchanging its two frequencies. Equivalently, canonicalize the entire side triples. Sorting the old pair and the new pair independently term by term is not legitimate. An exact example takes \(h=m=0,d_w=d_x=1\), side data \((1/10,0,1/4),(1/5,1,-1/2)\), and scalar correlations \((a,b,c)=(x,0,0)\). The square of the linear coefficient normalized by \(\pi/4\) is
\[
 \frac{2}{(97383/20000)^3}.
\]
Exchanging only the two new frequencies changes it to
\[
 \frac{2}{(52479/10000)^3}.
\]
Both linear coefficients are positive, so these different rational squares prove that the two primitives differ. This is a restriction on mode accounting, not an obstruction to the specified candidate.

**4. The endpoint issue and existence of a radius greater than one.**

**Lemma 5 (strict bound for the auxiliary polynomial).** For the specified \(B\), let \(b_n\) be its coefficients and \(q=|b_1b_{11}|>0\). Then
\[
 |B(z)|^2\le1-4q,\qquad |B(z)|\le1-2q<1\quad (|z|\le1). \tag{5}
\]

*Proof.* Exact rational inequalities give \(b_1<0,b_3>0,b_{11}<0\) and \(b_3\ge25|b_{11}|\). Since \(\sum|b_n|=1\), expansion of the squared modulus on \(z=e^{i\theta}\) expresses \(1-|B(z)|^2\) as a sum of nonnegative pair deficits. Retain only pairs \((1,3)\) and \((1,11)\); their contribution is
\[
 2|b_1|\big[b_3(1+\cos2\theta)+|b_{11}|(1-\cos10\theta)\big].
\]
The elementary inequality \(|\sin5t|\le5|\sin t|\), obtained from the sum of five unit complex numbers, implies \(1+\cos5x\le25(1+\cos x)\). The bracket is consequently at least \(2|b_{11}|\). This gives the first inequality. Since \(q\le1/4\), \(\sqrt{1-4q}\le1-2q\). The maximum modulus principle extends the bound to the disk. ∎

**Lemma 6 (no forced singularity from \(C(z)=z^7\)).** For every retained mode, the five polynomials in (3), after substitution \((a,b,c)=(P(z),B(z),C(z))\), satisfy on \(|z|\le1\)
\[
 D N_1 N_2\ne0,\qquad N_1N_2-tP_0Q_0\ne0\quad(0\le t\le1). \tag{6}
\]
There is a rational \(R>1\), possibly depending on the mode, for which (6) holds on \(|z|\le R\).

*Proof.* By Lemma 5 and \(\|P\|_1<1\), the \(W\)- and \(X\)-blocks of the Hermitian part of \(S\) are strictly positive. The \(Y\)-block is positive semidefinite, and can fail to be strictly positive only when \(c=1\) or \(c=-1\). If it is strictly positive, the accretivity argument in Lemma 4 applies.

At either exceptional value let \(E=(\ker S)^\perp\). The kernel is a one-dimensional subspace of the two \(Y\)-coordinates. It is also the kernel of the Hermitian part of \(S\). With respect to \(E\oplus\ker S\), the matrix \(S\) is \(S_E\oplus0\), with \(S_E\) strictly accretive. The matrix \(I-SJ\) is block upper triangular, with diagonal blocks \(I-S_EJ_{EE}\) and the identity. Its first block is invertible because \(S_E^{-1}-J_{EE}\) is strictly accretive. Moreover,
\[
 (I-SJ)^{-1}S=(S_E^{-1}-J_{EE})^{-1}\oplus0.
\]
The compression by \(U\) is strictly accretive: for a nonzero two-vector \(v\), \(U^Tv\) has \(W,W'\) coordinates \(v_1,v_2\), so its projection onto \(E\) cannot vanish. The diagonal and off-diagonal homotopy argument of Lemma 4 again applies to the compressed covariance. Multiplication by nonzero \(D\) proves (6).

All the expressions in (6) are polynomials in \(z\) and continuous in the additional real parameter \(t\in[0,1]\). They are nonzero on the compact disk times this compact interval. Uniform continuity therefore extends nonvanishing to some larger closed disk. A rational radius can be chosen inside that extension. There are only finitely many retained modes, so even a common such radius exists. This existence proof does not specify a computationally useful radius or a coefficient count. ∎

**5. A finite outward-arithmetic test and a proved scalar tail.**

Fix one mode and a rational radius \(R>1\). Substitute the exact polynomials \(P,B,C\) in (3), retaining the same names for the resulting scalar polynomials. The required certificate consists of complete, contiguous closed arcs covering \(z=Re^{i\theta}\), outward rectangles for each polynomial on each arc, and outward rectangles for the endpoint ratios of \(D,N_1,N_2\). Every polynomial rectangle for these last three must exclude zero. The principal endpoint argument increments must be justified and their summed enclosures must lie in \((-3,3)\).

Here is a complete sufficient test for the additional homotopy. Enclose
\(Z=P_0Q_0/(N_1N_2)\) on an arc. If the rectangle gives \(\operatorname{Re}Z\le x\), \(|\operatorname{Im}Z|\ge y\ge0\), and \(|Z|\le M_Z\), then
\[
 |1-tZ|\ge\delta:=\max\{\min(1,1-x),y/M_Z\}>0
 \quad(0\le t\le1), \tag{7}
\]
where the second term is zero when \(M_Z=0\). The first term follows by taking real parts. The second follows from the distance from the origin to the entire line \(1-tZ\), equal to \(|\operatorname{Im}Z|/|Z|\). Thus either term can justify a positive lower bound; no sampled maximum is used.

**Lemma 7 (boundary test implies the analytic branch and the entire tail).** Suppose the complete arc tests above succeed. If \(d_0,d_1,d_2>0\) bound the moduli of \(D,N_1,N_2\) below on an arc, \(p_*\) bounds \(|P_0|\) above there, and (7) holds, then
\[
 M_*:=\max_{\text{arcs}}\frac{(22/7)p_*}{4d_0\sqrt{d_1d_2\delta}}
\]
is a bound for the primitive on the circle. For its odd Taylor expansion \(F(z)=\sum_{j\ge0}f_jz^{2j+1}\),
\[
 \sum_{j>d}|f_j|\le
 T_d:=\frac{M_*R^{-2d-3}}{1-R^{-2}}. \tag{8}
\]

*Proof.* A zero-excluding rectangle is convex and lies in a half-plane avoiding zero. Consequently the continuous argument variation of a polynomial within each arc equals the certified principal endpoint increment; it cannot conceal an additional full turn. The implementation requires each endpoint ratio to have positive real part and imaginary-to-real ratio of absolute value at most \(1/2\). Its argument is enclosed by twelve arctangent Taylor terms with error \(|x|^{25}/25\). Summing increments gives the total argument variation, an integer multiple of \(2\pi\). Containment in \((-3,3)\) forces that integer to be zero. The argument principle therefore proves that \(D,N_1,N_2\) have no zeros in the disk.

On the boundary, (7) makes \(H_t=N_1N_2-tP_0Q_0\) nonzero for every \(t\in[0,1]\). Its winding is constant in \(t\); at zero it is the sum of the zero windings of \(N_1,N_2\). The argument principle now proves that \(H_t\) has no zeros in the disk for every \(t\). The disk is simply connected, so the square roots are fixed by the positive origin values. The ratio \(Z\) avoids \([1,\infty)\) throughout the disk, supplying the branch of the hypergeometric function in (4).

Euler's integral reads
\[
 {}_2F_1(1/2,1/2;2;Z)=\frac2\pi\int_0^1
        \sqrt{\frac{1-t}{t}}(1-tZ)^{-1/2}\,dt.
\]
The weight integrates to one. On each boundary arc its modulus is at most \(\delta^{-1/2}\); substituting this in (4) and using \(\pi<22/7\) gives \(M_*\). All inequalities can be rounded outward using rational square-root bounds. Cauchy's coefficient estimate gives \(|f_j|\le M_*R^{-(2j+1)}\). Summing the geometric series for \(j>d\) gives (8). ∎

The arc rectangles themselves are obtained as follows. For any input map \(Q(z)=\sum q_nz^n\), set \(L_2=\sum n^2|q_n|R^n\). At an arc midpoint \(\theta_0\) and half-width \(\eta\),
\[
 Q(Re^{i(\theta_0+s)})
 =Q(z_0)+is\sum_n n q_nz_0^n+E_s,
 \qquad |E_s|\le\eta^2L_2/2\quad(|s|\le\eta).
\]
The remainder may be enclosed in the complex rectangle with each component in \([-\eta^2L_2/2,\eta^2L_2/2]\). Substitution into the exact multivariate polynomials gives outward enclosures for every point of the arc. Failure of an enclosure test at a depth, time, or leaf limit means that the sufficient test is unresolved. It does not prove a singularity.

Formal coefficients are computed with \(y=z^2\), writing \(P(z)=zA(y), B(z)=zB_1(y), C(z)=zC_1(y)\). Even numerators acquire powers \(y^{(i+j+k)/2}\); \(P_0,Q_0\) divided by \(z\) acquire powers \(y^{(i+j+k-1)/2}\). The hypergeometric argument has zero constant term. Its coefficients are exact rationals
\[
 h_j=\frac{\binom{2j}{j}^2}{16^j(j+1)},\quad h_0=1,\quad
 h_j=h_{j-1}\frac{(2j-1)^2}{4j(j+1)}.
\]
Truncation of its composition at power \(j=d\) therefore computes the first \(d+1\) odd scalar coefficients exactly as formal expressions. Arithmetic intervals enclose those expressions; (8), rather than formal truncation, controls the remaining infinitely many coefficients.

**Lemma 8 (complete omission cost).** A weighted mode with coefficient \(w\) may be omitted in its entirety at cost
\[
 O=\frac{|w|}{(1+2ug)(1+2vg)} \tag{9}
\]
in the full substituted coefficient norm.

*Proof.* Its row function has modulus \(e^{-uR}\); all skew phases and the mixed linear phase have modulus one almost everywhere. Its Gaussian squared norm is \([(1+2ud_w)(1+2ud_x)]^{-1}\). This is at most \((1+2ug)^{-2}\) by \(d_w+d_x\ge2g\), \(d_wd_x\ge g^2\). The column bound is identical with \(v\). Lemma 2 gives (9), uniformly in all frequencies and in the mixing parameter. ∎

**6. The scalar coefficient implication for the classical constant.**

For any retained subset of the 120470 modes, sum the computed coefficient intervals *with their signed rational weights before taking absolute values*. At a common index \(j\), let the resulting enclosure be \([l_j,u_j]\), with zero placeholders beyond a mode's retained degree. Set
\[
 \gamma_{\rm cert}=l_0-
 \sum_{j\ge1}\max(|l_j|,|u_j|)
 -\sum_{\text{retained modes}}|w|T_d
 -\sum_{\text{omitted modes}} O-E. \tag{10}
\]
The finite sum in (10) ends at the largest retained degree; each shorter mode's entire missing part is already covered by its own tail. Every finite mode must occur once, either retained or omitted. Equations (1), (8), and (9) imply for the complete physical kernel
\[
 K(t)=\sum_{j\ge0}c_jt|t|^{2j},\qquad
 c_0-\sum_{j\ge1}|c_j|\ge\gamma_{\rm cert}.
\]
Only one copy of each error is subtracted: the functional \(c\mapsto c_0-\sum_{j\ge1}|c_j|\) is one-Lipschitz in the absolute coefficient norm.

**Proposition (sufficient Grothendieck criterion).** If the complete outward computation and accounting establish \(\gamma_{\rm cert}>0\), then
\[
 V_{\mathbb C}(A)\le\gamma_{\rm cert}^{-1}S_{\mathbb C}(A)
\]
for every finite complex matrix \(A\).

*Proof.* For any unit vectors, the simultaneous Gaussian construction in Lemma 1 and the scalar definition imply
\(\left|\sum a_{ij}K(\langle u_i,v_j\rangle)\right|\le S_{\mathbb C}(A)\).
The tensor vectors \(T_j(u_i),T_j(v_j)\) show that the modulus of \(\sum a_{ij}t_{ij}|t_{ij}|^{2j}\) is at most \(V_{\mathbb C}(A)\) for every \(j\). Choose vectors attaining the vector optimum, which is possible in dimension at most the number of rows plus columns. Triangle inequalities and absolute convergence then give
\[
 c_0 V_{\mathbb C}(A)
 \le S_{\mathbb C}(A)+\Big(\sum_{j\ge1}|c_j|\Big)V_{\mathbb C}(A).
\]
The positive gap proves the assertion. There is no fixed-dimension hypothesis, restriction on signs or phases of the entries, exchange of uncontrolled limits, or inverse-series assertion. ∎

**7. Exact finite calculation and reproducibility.**

Write \(O(r,u,\tau,s,v,\sigma;w)\) for (9). Retain precisely the modes
with \(O>10^{-10}\). All other nonzero canonical modes are omitted and charged
by (9). The two finite rational expansions supplied in `check_inventory.py`
produce the same 120,470 canonical modes, with 9,948 retained and 110,522
omitted. Their complete sorted key-and-weight digest, retained records and
exact omission sum are in `data/inventory.json.gz`. The unretained records
are reconstructed by either expansion; no list of selected numerical
successes determines the omitted set.

The allowed exact radii are
\(21/20,51/50,101/100,201/200,501/500\). The boundary cover starts with
64 equal angular panels and uses at most 18 subdivision levels, 16,384
accepted arcs and 32,768 enclosure evaluations for a radius. The finite
coefficient index \(d\) is at least 127 and satisfies (8), subject
to \(d\le1536\). The allowed production precisions are 768, 1,024,
1,536, 2,048 and 3,072 bits. The reproduction uses the successful reference
schedule in `data/reference_schedule.json.gz`: subtract eight from each
recorded replay degree and 128 from its precision to obtain the production
parameters. It uses the recorded successful radius directly. Each regenerated
contour must satisfy the unchanged complete-arc tests, and its prescribed
degree must still meet the tail criterion (8). For each mode the weighted complete
tail and total finite coefficient interval width are each at most
\(10^{-14}\). A retained mode that fails any bound remains required.

The fresh outward replay reconstructs all five determinant polynomials from
the exact mode key. It reevaluates all coefficients with 128 additional bits
and eight additional coefficient indices. It reevaluates every complete arc
on the production partition, including endpoint ratios. The production and
replay intervals must overlap at every common coefficient index; the final
calculation uses the replay intervals themselves. This replay repeats outward
arithmetic using the same mathematical kernels. The integer aggregation is
a separate implementation with its own coefficient decoder.

The coefficient enclosure format is the exact binary ball
\(m2^e\pm r2^f\), with integer mantissas and exponents and \(r\ge0\).
Each signed weighted interval is rounded outward once to the grid
\(2^{-320}\mathbb Z\). The resulting integer endpoints are added before
taking any absolute values. Each weighted tail is rounded up on that same
grid; the complete phase error and exact total omission cost are each
rounded up once. Such rounding can only decrease (10).

The reference calculation has 2,415,418 replay coefficient contributions.
The largest retained coefficient index is 1,010 (ordinary scalar degree
2,021). Its six exact rational terms are supplied in
`data/reference_validation.json` and independently in
`data/reference_aggregation.json`. Exact subtraction gives
\[
\gamma_*=
\frac{190075214760368137642162763178829164462637054815894556859938167637150500996946626676320683057935}
{266998379490113760299377713271194014325338065294581596243380200977777465722580068752870260867072}
>\frac{142379302169}{200000000000}>0.
\]
The proposition therefore gives the displayed upper bound. This comparison
uses the exact complete rational result; the decimal is only an outward
display of its reciprocal endpoint.

There are two execution dependencies. First, the recorded coefficient balls
and complete-arc rectangles must result from correct Python and FLINT/Arb
outward operations, including transcendental and formal-series operations.
Second, the exact contour validator must successfully apply the argument
and tail implications above to every retained mode. A hash identifies the
record on which an execution depends; it does not establish an enclosure.
Neither a report, a checkpoint flag, a sample, nor a progress count supplies
the missing analytic-value enclosure assertion.

The complete generation and fresh replay, complete contour validation, and
complete independent aggregation underlying the reference outputs have
executed. A complete 9,948-mode rerun with the portable entrypoints supplied
here has not executed. `check_exact.py` verifies reference metadata and exact
arithmetic without physical evaluations. The README gives the commands that
regenerate all raw evidence, validate every contour, independently aggregate
every contribution, and check the final rational comparison. Different
machines may require repeated invocations under the fixed work limits. The
prescribed parameters are unchanged, and every rerun must pass all fixed
inequalities and the displayed rational target. Equality of decimal
approximations is not an acceptance condition.

Resuming contour validation rechecks source, input and raw-file dependencies
and reconstructs coefficient widths, tails and integer contributions. It
reuses a previous full contour check only as an execution record of that
specified validator on those same raw bytes. A new empty validation output
directory forces every contour check again. The independent aggregation
similarly binds any reused exact contribution to its decoder, raw records
and complete validation descriptor. An empty output directory forces a
fresh independent decoding of every coefficient. These dependencies are
ordinary executable arithmetic guarantees, not a proof-assistant
formalization.

**Appendix: exact data.**

`sources/candidate_definition.json` specifies every coefficient of \(P,B,C\)
and every phase and damping parameter as an exact rational number. A
polynomial object maps an odd exponent to its coefficient. A profile is the
ordered triple \((e_j,\ell_j,u_j)\). `sources/full_task.json` specifies the
complete finite expansion, fixed omission threshold and every numerical
inequality and work limit used by the portable computation. Neither file
uses floating-point mathematical parameters.
