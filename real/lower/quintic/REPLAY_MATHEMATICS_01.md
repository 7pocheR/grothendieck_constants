# Verification of the six-dimensional inequalities

Let \(\nu=\sqrt{2/\pi}\), \(a=33/20\), and \(\lambda=3/20\).
The exact matrices, positive parameters \(t\), bounds \(D=t+\kappa\), and
rectangles used below are specified in `targets_01.json`.

For each parity matrix \(M\), verify every principal minor of
\(M\), \(M-\operatorname{diag}(1,\lambda,0)\), and
\(M-\operatorname{diag}(-a,1,\lambda)\). The first matrix is positive
definite; the other two are positive semidefinite. Exact rational elimination
gives \(M=N\Delta N^T\), where \(N\) is unit lower triangular and
\(\Delta\) has strictly positive rational entries. Thus
\(N\sqrt\Delta\) is the unique lower triangular Cholesky factor with
positive diagonal. Interleave the even and odd factors to obtain a six by
six matrix \(L\), with \(Q=LL^T\).

Write \(\psi_j=H_j/\sqrt{j!}\), where
\(H_0=1,H_1=s,H_{j+1}=sH_j-jH_{j-1}\), through degree five. Put
\(c_j(s)=\sum_i L_{ij}\psi_i(s)\). The quantities to bound are

\[
 b_j=\nu^{-1}\mathbb E|c_j(S)|,
 \qquad
 F(v)=\frac{\mathbb E|2\sum_jv_jc_j(S)+tS|
                   +\mathbb E|2\sum_jv_jc_j(S)-tS|}{2\nu}-\|v\|^2.
\]

The identity \((|p+q|+|p-q|)/2=\max(|p|,|q|)\) explains the
two integrals. Each polynomial has degree at most five. Neither a numerical
root calculation nor an unbounded numerical quadrature is needed.

On a rational interval \([l,u]\), ordinary outward interval Horner
evaluation contains every value of the polynomial, including uncertainty
in its algebraic coefficients. If the range has a fixed sign, integrate
with that sign using

\[
 J_0(s)=\Phi(s),\quad J_1(s)=-\phi(s),\quad
 J_j(s)=-s^{j-1}\phi(s)+(j-1)J_{j-2}(s).
\]

Direct differentiation gives \(J_j'=s^j\phi\), for each \(0\le j\le5\).
If the sign is undetermined, its integral is between zero and
\(\sup_{[l,u]}|p|[\Phi(u)-\Phi(l)]\). Bisect such intervals; sum every
remaining uncertainty without assuming that the polynomial has simple roots.
The initial 24 intervals partition \([-12,12]\). Signed intervals and
unresolved intervals continue to partition this entire closed interval.

Both infinite tails are bounded by \(\sum_{j=0}^5|p_j|T_j(12)\), where

\[
 T_0(R)=\operatorname{erfc}(R/\sqrt2),\quad
 T_1(R)=\nu e^{-R^2/2},\quad
 T_j(R)=\nu R^{j-1}e^{-R^2/2}+(j-1)T_{j-2}(R).
\]

These are the exact two-sided absolute Gaussian moments beyond \(R\),
as follows by integration by parts. The entire tail, including every
possible change of sign beyond \(R\), is included. A sign test is used
only when the interval comparison proves it. A floating heap priority
affects only which interval is bisected next. Failure to reduce the
uncertainty within 20,000 bisections is an error, not a certificate.

Every numerical operation on algebraic coefficients and elementary functions
uses an outward Arb interval at 256 bits. Endpoints are converted exactly
to rational numbers and rounded outward to a dyadic grid with denominator
\(2^{52}\); no floating conversion enters a reported inequality.
The sum of unresolved integral upper bounds is at most \(2^{-48}\).
The antiderivative intervals and complete tail bounds are additionally
included in the reported upper bound. Thus the tolerance is a subdivision
stopping condition, not an omitted error.

It is enough to prove \(b_j\le B_j\) for all six coordinates and
\(F\le D\) on
\([0,B_0]\times[0,B_1]\times\prod_{j=2}^5[-B_j,B_j]\).
The validity of this compact reduction is established in the full proof.
For a rectangle with side lengths \(w_j\), convexity of
\(F(v)+\|v\|^2\), followed by multilinear interpolation, gives

\[
 F(v)\le\max_{z\text{ a vertex}}F(z)+\frac14\sum_{j=0}^5 w_j^2.
\]

The interpolation error in coordinate \(j\) is exactly
\((v_j-l_j)(u_j-v_j)\le w_j^2/4\). Thus verifying the 64 vertices
of every terminal rectangle proves its full closed interior and boundary.
The rectangles are reconstructed by repeatedly bisecting the longest
side, taking the smallest coordinate index in a tie.

For each of the 16 initial binary paths of length four, the terminal paths
must be distinct, prefix free, and satisfy
\(\sum 2^{-(|p|-4)}=1\). Since the family is finite, expanding all paths
to their common maximum depth proves that every descendant at that depth
is covered. The closed rectangles therefore cover the initial rectangle,
including every shared face. All 16 initial rectangles cover the root
rectangle. The verification uses each distinct required vertex once and
reconstructs every terminal correction as an exact rational number.

The original proposed integral values are not premises. Fresh interval
values at every required vertex are substituted in the 64-vertex maximum.
The final inequality is checked again for every terminal rectangle. A
comparison with proposed vertex values is reported separately and is not
required when the newly verified rectangle bound itself proves \(F\le D\).
