# A finite scalar coefficient calculation with global error control

All formulas below use standard real Gaussian measure and the
orthonormal probabilists' Hermite polynomials
`psi_k=He_k/sqrt(k!)`.  Let

    B(x)=ell*x+sum_l c_l sin(omega_l*x),
    S(x)=sum_{i=1}^d B(x_i),
    A(s)=sum_{l=1}^J a_l sin(2*l*s),  U(s)=s+A(s),
    f(w,x)=sign(w+U(S(x))).

The exact input has d=4 and J=4. The primary preprocessing P
is the exact polynomial in `certificate/physical.json`;
Q(z)=z.  Its absolute coefficient sum is 1-10^-12.  The target scalar
function, with the physical auxiliary reflection, is

    H(z)=(pi/2) sum_{j,beta} A_{j,beta}^2 P(z)^j (-z)^|beta|.  (1)

Here A_{j,beta} denotes the Hermite coefficient of f, rather than the
outer profile A(s).  No unnormalized coefficient table may be inserted
into (1).  Only j+|beta| odd can occur.  Since P(0)=0, coefficients of
H through scalar degree N depend only on j+|beta|<=N.

## A finite approximation of the primary conditional coefficients

Set c=sqrt(2/pi) and phi(u)=exp(-u^2/2)/sqrt(2*pi).  For j>=1,

    q_j(s)=E_W sign(W+s) psi_j(W)
          =c*i^(j-1) integral V_j(u) exp(i*u*s) du,
    V_j(u)=phi(u)*u^(j-1)/sqrt(j!).                         (2)

For j=0 the same formula holds with V_0(u)=phi(u)/u, as a symmetric
principal value, and q_0(s)=erf(s/sqrt(2)).  Differentiating (2)
recovers c*exp(-s^2/2), so its principal-value constant is fixed by
q_0(0)=0.

Choose an integer M>=0 and replace exp(i*u*A(s)) by its Taylor
polynomial of degree M.  Define exact real Laurent coefficients by

    R(z)=sum_l (a_l/2)(z^l-z^(-l)),
    R(z)^m/m! = sum_k r_{m,k} z^k.

The resulting weight and function are

    V_j^[M](u)=sum_{m=0}^M sum_k
                    r_{m,k}(u-2k)^m V_j(u-2k),
    q_j^[M](s)=sum_{m=0}^M A(s)^m q_j^(m)(s)/m!
             =c*i^(j-1) integral V_j^[M](u) exp(i*u*s) du.  (3)

Only the m=0,j=0 term in (3) has a principal value.  Every other term
is a Gaussian times a polynomial.  This finite expansion retains all
signed coefficients r_{m,k}, including cancellations between shifts.
Parity gives V_j^[M](-u)=(-1)^(j-1)V_j^[M](u).

Writing A_*=sum_l |a_l|, the real phase-Taylor remainder in (2) gives

    sup_s |q_j(U(s))-q_j^[M](s)| <= e_j,
    e_j=c*A_*^(M+1) E|Z|^(j+M)/((M+1)! sqrt(j!)).            (4)

The formula is also valid for j=0: the extra power from the remainder
removes the principal-value singularity.  Gaussian Bessel's inequality
therefore bounds the Euclidean error of the complete coefficient array
with j+|beta|<=N by

    delta_phase=(sum_{j=0}^N e_j^2)^(1/2).                  (5)

No factor counting auxiliary multiindices belongs in (5).

## The exact finite polynomial, without enumerating auxiliary indices

Let L>0, Delta=2*pi/L, and u_n=n*Delta for n=1,...,F.  Put

    h_k(u)=E psi_k(X) exp(i*u*B(X))=i^k p_k(u).

Since B is odd, p_k is real.  Define the following polynomials,
truncated after degree N:

    K_{uv}(z)=sum_{k=0}^N p_k(u)p_k(v) z^k,
    W_{uv}(z)=sum_{j=0}^N V_j^[M](u)V_j^[M](v) z^j.

The contribution of all nonzero Fourier nodes to the finite scalar
polynomial is exactly

    4*Delta^2 * Odd_{<=N} sum_{u,v>0}
                       W_{uv}(P(z)) K_{uv}(-z)^d.          (6)

The sum in (6) runs over the F chosen positive nodes.  A triangular
implementation uses weight 1 for u=v and 2 for u<v.  Formula (6)
follows by squaring each finite Fourier sum before summing the
auxiliary multiindices.  The coefficient of z^b in K(z)^d equals
the sum of the corresponding products over all |beta|=b.  Projection
onto odd powers removes exactly the forbidden j+b even terms.  The
normalization c^2*pi/2=1 explains the absence of an additional pi/2
in (6); it does not change the normalization of (1).

There are two zero-node corrections.  For odd j>=1, define

    a_j=2*Delta*sum_{u>0} V_j^[M](u) p_0(u)^d,
    z_j=Delta*V_j^[M](0).

Add sum_{odd j<=N}(2*a_j*z_j+z_j^2)P(z)^j to (6).
For odd k>=1, define

    b_k=2*Delta*sum_{u>0} V_0^[M](u) p_k(u)p_0(u)^(d-1),
    t_k=Delta*phi(0)*p_k'(0),
    p_k'(0)=ell*1_{k=1}
               +sum_l c_l omega_l^k exp(-omega_l^2/2)/sqrt(k!).

Add -d*sum_{odd k<=N}(2*b_k*t_k+t_k^2)z^k.  The minus sign is the
auxiliary reflection in (1).  To see the second correction, the
principal-value quotient at zero has a nonzero limit only when beta
has exactly one positive entry, and that entry is odd.  Its limit is
phi(0)*p_k'(0).  Equivalently, the j=0 Fourier approximation includes
the linear term 2*S/L.  If two entries of beta are positive the limit
is zero; beta=0 is excluded by oddness.  These arguments also cover
vanishing coefficients and k for which p_k'(0)=0.

Composition and products in (6) are finite operations in real ball
power series.  Computing W(P) by power-series composition and K(-z)^d
by multiplication uses only N+1 coefficients per series.  The number
of node pairs, rather than the number of degree-N multiindices in d
coordinates, controls the outer computation.

## Uniform control of the Fourier lattice

The following estimates bound errors of functions of the actual S.
They do not assume S is Gaussian.  Put

    sigma^2=d*ell^2, D=d*sum_l |c_l|.

Then S=G+E, where G is a centered Gaussian of variance sigma^2 and
|E|<=D; E may depend on G and on all auxiliary coordinates.

Fix 0<t<1, set a=1+t, and write C_t=(1-t^2)^(-1/4).  Mehler's
identity gives, for every n>=0 and real x,y,

    |He_n(x+i*y)| exp(-(x^2-y^2)/2)/sqrt(n!)
        <= C_t*t^(-n/2)
                exp(-x^2/(2*a)+y^2/(2*(1-t))).              (7)

Define A_y=sum_l |a_l| cosh(2*l*y) and
A'_y=sum_l 2*l*|a_l| cosh(2*l*y).  For j>=1 put

    K_j(t,A_y)=c*C_t sum_{m=0}^M A_y^m/m!
                         sqrt((j+m-1)!/j!)*t^(-(j+m-1)/2).

For the derivative of q_0^[M], put

    K'_0(t,A_y,A'_y)=c*C_t [
        sum_{m=0}^M A_y^m/m! sqrt(m!)*t^(-m/2)
      + A'_y sum_{m=1}^M A_y^(m-1)/(m-1)!
                                   sqrt((m-1)!)*t^(-(m-1)/2)].

Equation (7) yields Gaussian envelopes for q_j^[M] and (q_0^[M])'
with these constants, multiplied by exp(y^2/(2*(1-t))).

For a real-space alias bound choose 0<beta<1 and define

    C=beta*D^2/(2*a*(1-beta)),
    kappa=beta/(a+2*beta*sigma^2),
    cL=kappa*L^2/2,
    mu=sqrt(beta*sigma^2/(2*a)).

Assume 2*exp(-3*cL)<1.  The L2 errors from replacing the Fourier
integrals by their complete infinite lattices are at most

    alias_j=2*K_j(t,A_*)*exp(C-cL)/(1-exp(-3*cL)), j>=1,
    alias_0=2*K'_0(t,A_*,A'_0)*(D+sigma+mu*L)
                           *exp(C-cL)/(1-2*exp(-3*cL)).      (8)

Here and below a zero variance is interpreted directly, with sigma=mu=0;
the specified input has sigma>0.

For completeness, |S+mL|>=(|G+mL|-D)_+ and
(|x|-D)_+^2>=beta*x^2-beta*D^2/(1-beta).  Gaussian square completion
then gives ||q_j^[M](S+mL)||_2<=K_j exp(C-cL*m^2); dropping its
additional factor at most one gives (8) for j>=1.  For j=0, Poisson
summation applies to (q_0^[M])', whose integral is 2.  The error of its
integrated, odd Fourier series is

    sum_{m!=0} integral_0^1 S*(q_0^[M])'(v*S+mL) dv.

Square completion with v*G+mL gives variance at most sigma^2, mean
of absolute value at most mu*|m|*L, and exponential bound
exp(C-cL*m^2), uniformly for 0<=v<=1.  Minkowski's inequality bounds
the remaining factor |S| by D+sigma+mu*|m|*L.  The ratio of successive
positive summands is at most 2*exp(-3*cL).  This proves (8), including
the linear term 2*S/L.  This bound avoids a split at |S|=L/2.

For a frequency cutoff choose y>0 and any t in (0,1), independently
of the choice in (8).  Set

    C_j(y)=K_j(t,A_y)*sqrt(1+t)/2 * exp(y^2/(2*(1-t))),
    C'_0(y)=K'_0(t,A_y,A'_y)*sqrt(1+t)/2
                                      *exp(y^2/(2*(1-t))).

Shifting the Fourier contour in (3), or in its derivative for j=0,
gives |V_j^[M](u)|<=C_j(y)exp(-y|u|), and
|V_0^[M](u)|<=C'_0(y)exp(-y|u|)/|u| for u!=0.  Thus the omitted
nonzero lattice nodes have uniform function errors

    freq_j <= 2*c*Delta*C_j(y)*exp(-y*(F+1)*Delta)
                                      /(1-exp(-y*Delta)), j>=1,
    freq_0 <= 2*c*Delta*C'_0(y)*exp(-y*(F+1)*Delta)
                          /((F+1)*Delta*(1-exp(-y*Delta))).  (9)

Every shifted Gaussian in (3) is included in (7)-(9); no finite shift
window or fitted frequency decay is being assumed.  Taking minima of
different rigorously evaluated choices of t and y is permitted.

## Computing the Hermite characteristic functions

For a common x lattice h_x*Z, truncate at |x|<=A_x=K_x*h_x.  A bin
implementation may use any integer period M_b and bin spacing L/M_b.
Assign each B(k*h_x) to an integer bin j_k and rigorously check
|B(k*h_x)-j_k*L/M_b|<=delta_b.  Keep the signed Hermite weights in
each bin.  Taylor expansion in the real residual through degree J_b
followed by the positive-frequency DFT gives the finite characteristic
values.  The DFT index is n modulo M_b, while the Taylor factor uses
the full frequency u_n; n may exceed M_b.  There is no 1/M_b factor
in this DFT.  All branch choices for bin integers are checked against
whole Arb intervals before they are used.

The following bounds control every degree k<=N at every |u|<=T=F*Delta.
For a strip height y_x>0, let

    D_k(y_x)=sum_{r=0}^k binom(k,r)y_x^(2r)/r!,
    B_im(y_x)=|ell|y_x+sum_l |c_l|sinh(|omega_l|y_x).

The exact Hermite translation identity gives
E|psi_k(X+i*y_x)|^2=D_k(y_x).  The infinite x-lattice error is at most

    ex_k=2*sqrt(D_k(y_x))*exp(y_x^2/2+T*B_im(y_x))
                                      /(exp(2*pi*y_x/h_x)-1). (10)

Using (7) at y=0, the discarded x-lattice tail is at most

    tx_k=C_t*t^(-k/2)*sqrt(1+t)
                          *erfc(A_x/sqrt(2*(1+t))).          (11)

The decreasing Gaussian envelope proves the lattice tail bound in
(11), without requiring A_x^2>=k or expanding Hermite coefficients
in absolute value.  To bound the sum of absolute bin weights, use
Cauchy--Schwarz on the complete x lattice and the same strip bound
for psi_k^2.  It is at most

    mass_k <= sqrt((1+g_x)*(1+g_x*D_k(y_x))),
    g_x=2*exp(y_x^2/2)/(exp(2*pi*y_x/h_x)-1).

The phase-bin remainder is therefore

    bx_k=mass_k*(T*delta_b)^(J_b+1)/(J_b+1)!.                (12)

Equations (10)-(12), together with the explicit finite Arb radius,
enclose every h_k.  The positive formulas for D_k and t^(-k/2) are
nondecreasing in k, so the k=N error also bounds all lower degrees.
Every transform here is a transform of psi_k, with exactly one
factorial normalization. Transforms of differentiated boundary weights
cannot be substituted for these transforms.

## Propagating all errors to the finite scalar coefficients

Let eta(u) bound the Euclidean error of the vector
(h_0(u),...,h_N(u)); for example use sqrt(N+1) times the largest
of (10)+(11)+(12).  Bessel's inequality gives ||h(u)||_2<=1.
The error of its d-fold tensor power is consequently at most

    (1+eta(u))^d-1.                                        (13)

For each j let b_j be the bound

    b_j=2*c*Delta*sum_{n=1}^F |V_j^[M](u_n)|
                                  *((1+eta(u_n))^d-1).

The zero-node corrections are evaluated separately by finite Arb
arithmetic. Alternatively, for a uniform positive bound one
may replace each |V_j^[M](u_n)| by C_j(0), or C'_0(0)/u_n for j=0.
This gives explicit b_j before computing any Hermite transforms.

The total Euclidean error between the actual Hermite coefficient
array of total degree at most N and the finite array defining (6)
with its two corrections is bounded by

    delta <= delta_phase
               + (sum_{j=0}^N (alias_j+freq_j+b_j)^2)^(1/2). (14)

For the true coefficient array its Euclidean norm is at most 1.
Hence sum|A_{j,beta}^2-Atilda_{j,beta}^2|<=delta*(2+delta).
The absolute coefficient norms of P and Q are at most 1.  Composition,
total-degree restriction, and restriction to scalar degree N are
contractions for the relevant absolute coefficient sums.  Therefore
the total error in the scalar coefficients through N is at most

    E_head=(pi/2)*delta*(2+delta).                           (15)

Finite ball arithmetic in (6) and the corrections supplies its own
coefficient intervals and is not omitted from the final bound.
In (14)-(15) delta controls only the additional analytic remainders.

The specified calculation uses N=403, M=31, L=40, F=768,
h_x=1/256, K_x=5120, y_x=1/8, t=99/100 for the real-space
envelopes, and beta=1/3 for the aliases. The auxiliary dimension is four.
The additional conservative bin parameters are M_b=2048, delta_b=1/100,
and J_b=32; their error is retained for the direct spatial calculation.
Every frequency-contour choice is provided as an exact rational t and y
in certificate/head_error_input.json. There are 404 such choices.

The complete formulas (4), (8), (9), and (10)--(15), evaluated with
outward arithmetic, give the exact E_head in certificate/scalar_bounds.json.
The final coefficient objective subtracts this allowance once. It then
subtracts the separately computed localized tail from mathematical_note.md;
no preliminary fourth-order coefficient tail is added to that tail.

For direct spatial summation, replace the finite bin implementation by

    p_k(u) = h_x sum_(v=-K_x)^K_x phi(v*h_x) psi_k(v*h_x)
                 * i^(-k) exp(i*u*B(v*h_x)).

Oddness of B and Hermite parity make this expression real: use
(-1)^(k/2)cos for even k and (-1)^((k-1)/2)sin for odd k.
The positive half lattice has multiplicity two except at zero. This is
exactly the finite expression in source/head/primitives.py. Its integral
error is bounded by (10) and (11), and retaining the nonnegative (12)
remains valid. Interval roundoff in these finite sums is carried by the
coefficient intervals and is separate from E_head.
