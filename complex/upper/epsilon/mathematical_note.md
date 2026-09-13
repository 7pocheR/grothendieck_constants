# An amplitude refinement of the coupled complex upper bound

**Unpublished computer-assisted research.** The exact coefficient estimate
below gives
\[
 K_G^{\mathbb C}\le 1/\Gamma_*<1.404698554831.
 \tag{1}
\]
Only one bounded phase amplitude is varied. All scalar degrees, the
complete phase remainder, and every omitted Fourier mode remain accounted
for. The numerical primitives and their contour bounds are unchanged.

## Physical function, normalization and admissibility

Use standard proper complex Gaussians, with density
\(\pi^{-1}e^{-|z|^2}\) and \(\mathbb E|Z|^2=1\). Put
\[
 L=W+hX+mY,\quad q_X=\operatorname{Im}(X\overline W),\quad
 q_Y=\operatorname{Im}(Y\overline W),\quad
 R=d_w|W|^2+d_x|X|^2,
\]
\[
 \Psi_0(W,X)=\sum_{j=1}^{5}e_j e^{-u_jR}\sin(2\ell_jq_X),\qquad
 S(W,Y)=\sin(2\lambda q_Y),
\]
\[
 f_e(W,X,Y)=\operatorname{ph}(L)
           \exp(-2ikq_X+i\Psi_0+i e S),\qquad
 g_e(W,X,Y)=f_e(W,-X,-Y).
 \tag{2}
\]
All rational parameters are in `certificate/baseline_candidate.json`.
That file defines the original amplitude \(e_0=1/100\); in (2) we use
\[
 e=e_*:=e_0+2^{-28}=\frac{67108889}{6710886400},\qquad
 m=-\frac1{50},\quad \lambda=\frac14.
 \tag{3}
\]
In particular, changing \(e\) does not impose the relation \(m=-2e\).
The value of \(\operatorname{ph}(0)\) is irrelevant because
\(\mathbb E|L|^2=1+h^2+m^2>0\). The functions have modulus one
almost everywhere, have common rotational charge one, and commute with
complex conjugation.

The exact real odd polynomials \(P,B,C\) in the candidate remain fixed.
They satisfy \(\|P\|_1<1\), \(\|B\|_1=\|C\|_1=1\), and
\(C(z)=z^7\). An odd polynomial \(Q(z)=\sum q_{2n+1}z^{2n+1}\)
acts on a complex correlation by
\(\widetilde Q(t)=\sum q_{2n+1}t|t|^{2n}\).
The tensor map \(x\mapsto x^{\otimes(n+1)}\otimes\overline x^{\otimes n}\)
has cross inner product \(t|t|^{2n}\). Orthogonal direct sums weighted
by \(\sqrt{|q_{2n+1}|}\), with the sign on the column side, realize
\(\widetilde Q\). Two orthogonal padding directions complete norms
when \(\|Q\|_1<1\). Three independent proper Gaussian fields therefore
realize the correlations \(\widetilde P(t),\widetilde B(t),\widetilde C(t)\)
jointly for every finite vector family.

Let \(K_e(t)\) be the covariance of (2) under these correlations.
Complex Hermite orthogonality, charge one, and conjugation symmetry give
\[
 K_e(t)=\sum_{n\ge0}c_n(e)t|t|^{2n},\qquad
 c_n(e)\in\mathbb R,\quad \sum_n|c_n(e)|\le1.
 \tag{4}
\]
More generally the coefficient absolute norm of the substituted covariance
of two \(L^2\) functions is at most the product of their \(L^2\) norms.
Indeed, the complete complex Hermite basis diagonalizes the three Gaussian
correlations; Cauchy--Schwarz bounds the sum of absolute products of
coefficients. Substitution by each polynomial and its conjugate contracts
that norm. This includes singular endpoint correlations by uniform
approximation with finite Hermite expansions. There is no factor of
\(\pi\), two, or the number of coordinates in this contraction.

For a finite complex matrix, let \(S\) be its optimum over scalar unit
phases and \(V\) its optimum over unit vectors in arbitrary finite
complex dimension, with inner products linear in the first argument.
Choose the vector family with positive real optimum \(V\). Gaussian
rounding with (2) and the preceding joint tensor construction gives
\[
 S\ge \operatorname{Re}\mathbb E\sum a_{ij}\sigma_i\overline{\tau_j}
 \ge \left(c_0(e)-\sum_{n\ge1}|c_n(e)|\right)V.
 \tag{5}
\]
For the last inequality each tensor expression \(t|t|^{2n}\) is again
a unit-vector correlation, so its bilinear form has modulus at most
\(V\). Absolute convergence in (4) justifies summation. Thus any positive
lower bound for the parenthesis in (5) gives \(K_G^{\mathbb C}\le1/\Gamma\).

## Complete finite expansion and six weight classes

Expand the original bounded phase \(\Psi_0\) through total
row-plus-column order four, and independently expand the new phase
\(eS\) through total row-plus-column order four. A primitive is indexed
by the complete side triples
\(\kappa=(r,u,\tau;s,v,\sigma)\), and has integrand
\[
 \operatorname{ph}(W+hX+mY)
 \overline{\operatorname{ph}(W'-hX'-mY')}
 e^{-uR-vR'}e^{-2i(rq_X+sq'_X+\tau q_Y+\sigma q'_Y)}.
 \tag{6}
\]
Every primitive in (6) is independent of \(e\). Only its rational
weight changes. The determinant and hypergeometric functions used to
evaluate it, its odd scalar series, and its analytic contour therefore
remain the same when (3) is applied.

For completeness, the primitive is determined by the proper complex
Gaussian covariance \(\Sigma\) of the six variables and the quadratic
exponent matrix \(J\) in (6). Its diagonal damping entries are
\((-ud_w,-ud_x,-vd_w,-vd_x,0,0)\) in the order
\((W,X,W',X',Y,Y')\). Its nonzero skew entries are
\(J_{12}=-r,J_{21}=r,J_{15}=-\tau,J_{51}=\tau\), and the analogous
column-side entries \(J_{34}=-s,J_{43}=s,J_{36}=-\sigma,J_{63}=\sigma\).
The only cross-side covariance entries of \(\Sigma\) are \(a,b,c\)
for the matching coordinates. With
\[
 U=\begin{pmatrix}1&h&0&0&m&0\\0&0&1&-h&0&-m\end{pmatrix},
\quad D=\det(I-\Sigma J),\quad
 \begin{pmatrix}N_1&P_0\\Q_0&N_2\end{pmatrix}
 =U\operatorname{adj}(I-\Sigma J)\Sigma U^T,
\]
the analytic continuation from the real scalar germ is
\[
 F_\kappa(a,b,c)=\frac{\pi P_0}{4D\sqrt{N_1}\sqrt{N_2}}
 {}_2F_1\left(\frac12,\frac12;2;\frac{P_0Q_0}{N_1N_2}\right).
 \tag{7}
\]
The square roots continue from positive origin values. Gaussian tilting
has determinant mass \(D^{-1}\), and the two-coordinate phase integral
has factor \(\pi/4\), giving (7). Its use requires nonvanishing of
\(D,N_1,N_2\) and of \(N_1N_2-tP_0Q_0\) throughout the relevant
continuation for \(0\le t\le1\). The original full contour certificates
and their analytic branch argument remain numerical and mathematical
dependencies; changing the weight does not establish new contours.

Put \(i=-\tau/\lambda\), \(j=-\sigma/\lambda\), which are integers
with \(|i|+|j|\le4\). The exact ordered new-phase weight is
\[
 N_{ij}(e)=s_i s_j
 \sum_{\substack{p,q\ge0\\|i|+|j|+2p+2q\le4}}
 \frac{(-1)^{p+q}(e/2)^{|i|+|j|+2p+2q}}
      {p!(p+|i|)!q!(q+|j|)!},
 \tag{8}
\]
where \(s_i=1\) for \(i\ge0\) and \(s_i=(-1)^{|i|}\) for \(i<0\).
There are exactly 41 ordered pairs. Direct finite expansion of (8)
shows that the weight ratio is one of six functions
\(r_c(e)=p_c(e)/p_c(e_0)\), with
\[
 \begin{aligned}
 p_{00}(e)&=1-e^2/2+3e^4/32,& p_{10}(e)&=e-3e^3/8,\\
 p_{20}(e)&=e^2-e^4/3,&p_{11}(e)&=e^2-e^4/4,\\
 p_3(e)&=e^3,&p_4(e)&=e^4.
 \end{aligned}
 \tag{9}
\]
The first four classes correspond to the unordered pairs of magnitudes
\(\{0,0\},\{1,0\},\{2,0\},\{1,1\}\); classes 3 and 4 correspond
to the indicated total magnitude. Sign and constant factors cancel in
the ratio. Every polynomial in (9) is positive throughout
\([9/1000,11/1000]\), so no nonzero weight disappears or previously
cancelled term appears on that interval.

Exchange of the complete side triples preserves the primitive and the
ratio. Canonicalization must exchange all three coordinates of each side
together. Sorting the old and new frequency pairs independently is not
valid. The original expansion contains 2,962 old modes, including 54
equal-side cases. The full number after the new expansion is
\((2962-54)41+54\cdot23=120470\). The retained subset has 9,948 modes;
every other 110,522 mode is charged below.

For each class, let
\[
 A_c(z)=\sum_{\kappa\text{ retained in class }c}
               w_\kappa(e_0)F_{\kappa,\mathrm{finite}}(z).
 \tag{10}
\]
Signed weights are summed before taking any absolute values. All of
each mode's validated scalar coefficients are used. The file
`certificate/six_groups.json` supplies interval endpoints for (10), on
the grid \(2^{-320}\), and a complete weighted scalar-tail upper bound
for each class. There are 3,381 stored class intervals in total. The
largest class has radial indices 0 through 1,010, corresponding to
ordinary odd powers through 2,021. A shorter class has zero finite
coefficients beyond its saved length, while its entire omitted scalar
tail remains charged. There is no common truncation that discards a
longer mode's saved coefficients.

The complete finite approximation at \(e\) is exactly
\(\sum_c r_c(e)A_c\). For positive ratios, its endpoint intervals are
obtained by the same positive weighted sums of class endpoints. This is
a finite-amplitude evaluation, not a derivative approximation.

## All infinite errors

Let \(g=999/1000\). The exact parameters satisfy
\(d_wd_x\ge g^2\), hence \(R\ge2g|q_X|\). The bound
\[
 M=\sum_j|e_j|
 \begin{cases}
 1,&u_j=0,\\
 \min(1,2|\ell_j|/(5u_jg)),&u_j>0
 \end{cases}
 =\frac{19007117794841543023}{399600000000000000000}
 \tag{11}
\]
majorizes \(|\Psi_0|\). For a damped term use
\(|\sin v|\le\min(1,|v|)\), maximize \(s e^{-2u_jgs}\), and use
\(1/\mathrm e<2/5\). The new phase satisfies \(|S|\le1\).

The coefficient contraction underlying (4), followed by the integral
Taylor remainder in the coefficient absolute norm, gives the complete
rectangular phase error
\[
 E_{\rm phase}(e)=\frac{(2M)^5}{120}
       +\left(\sum_{j=0}^{4}\frac{(2M)^j}{j!}\right)
                         \frac{(2|e|)^5}{120}.
 \tag{12}
\]
The factors two count the row and column phases. Derivatives with
respect to the two bounded phase amplitudes are bounded in \(L^2\);
the fixed singular phase is not differentiated. First truncate the
old phase at full new amplitude, then truncate the new phase in each
retained old derivative. This proves (12) for all scalar degrees. No
extra quadratic finite-step remainder is required in addition to (12).

For a mode whose ordinary odd series is retained through radial index
\(d\), a complete analytic circle bound \(|F_\kappa(z)|\le M_\kappa\)
at \(|z|=\rho_\kappa>1\) gives
\[
 T_\kappa=\frac{M_\kappa\rho_\kappa^{-2d-3}}
                      {1-\rho_\kappa^{-2}}.
 \tag{13}
\]
The class tail \(T_c\) is the sum of all accepted
\(|w_\kappa(e_0)|T_\kappa\) in that class, rounded upward when
stored on the dyadic grid. Under (9) it becomes \(r_c(e)T_c\).
Each original contour must cover the entire circle and establish the
nonvanishing/branch conditions in (7). For example, Euler's integral
for the hypergeometric function gives a full circle bound
\((22/7)p_*/[4d_0\sqrt{d_1d_2\delta}]\) from lower modulus bounds
\(d_0,d_1,d_2\), an upper bound \(p_*\) for \(|P_0|\), and
\(|1-tP_0Q_0/(N_1N_2)|\ge\delta>0\) for every \(t\in[0,1]\).
The original winding and homotopy checks are retained as premises.

For omitted modes, the coefficient contraction and Gaussian damping give
\[
 \|w_\kappa F_\kappa\|_1\le
 \frac{|w_\kappa|}{(1+2ug)(1+2vg)}.
 \tag{14}
\]
Indeed, the squared \(L^2\) norm of the row damping is
\([(1+2ud_w)(1+2ud_x)]^{-1}\le(1+2ug)^{-2}\), and similarly on
the column side. Let \(O_0\) be the complete sum of (14) for the
110,522 omitted modes at \(e_0\), given as an exact rational in
`certificate/endpoint.json`. Their new total is bounded by
\[
 O(e)=\max_c r_c(e)\,O_0.
 \tag{15}
\]
It would not be valid to infer separate class omission sums from this
single total. Positivity of all ratios is what permits (15).

Let \([L_n(e),U_n(e)]\) be the signed class combination at radial
index \(n\). With \(\operatorname{up}_{320}\) denoting upward
rounding to an integer multiple of \(2^{-320}\), the final sufficient
objective is
\[
 \begin{split}
 \Gamma(e)={}&L_0(e)-\sum_{n\ge1}\max(|L_n(e)|,|U_n(e)|)\\
 &-\sum_c r_c(e)T_c
   -\operatorname{up}_{320}(E_{\rm phase}(e))
   -\operatorname{up}_{320}(O(e)).
 \end{split}
 \tag{16}
\]
All rounding widths remain in the intervals. Since the objective is
1-Lipschitz in coefficient absolute norm, each complete error in (16)
is subtracted once. Every retained mode contributes all its available
coefficients and scalar tail; every omitted mode contributes its full
norm bound. Unvisited modes are not assigned zero without an omission
charge.

At (3), exact rational evaluation of (16) gives the fraction
\(\Gamma_*\) in `certificate/endpoint.json`, approximately
\(0.7118965108641602\). Its reciprocal is approximately
\(1.4046985548308357\), and exact cross multiplication proves the
strict decimal endpoint in (1). The coefficient lower bound at the
original amplitude is also recomputed from the same six arrays. A
rational trial amplitude suffices for this implication; no global
optimization or stationarity claim is needed.

## Computational scope

The original 9,948-mode determinant, coefficient, and complete contour
calculations remain part of the result. A complete grouping execution
used all 2,415,418 validated coefficient contributions in 78 blocks.
Their exact block sums and the selected finite objective were checked
separately. This grouping did not constitute a fresh native decoding of
all original coefficients, and no new physical integral was needed for
the amplitude change.

The small standard-library verifier here reproduces (8)--(9), the
signed combination of all six supplied arrays, (11)--(12), (15)--(16),
and the exact reciprocal comparison. It does not reconstruct the group
arrays from the absent 78 blocks, reprove the primitive enclosures or
reevaluate their contours. The original omission sum is also retained
as a complete numerical/algebraic execution premise, rather than
silently replaced by a partial support calculation.

The omitted data and their identities are listed in
`certificate/external_dependencies.json`. The included exact and Arb
source modules, together with the preceding coupled package's proof and
portable generation/validation interface, specify how to extend the
check to the original numerical components. `reproduction.md` distinguishes
those operations from the endpoint check available in this small bundle.
The refinement is a small improvement of a sufficient upper bound; it
does not determine either Grothendieck constant or establish a substantial
reduction of the remaining gap by itself.
