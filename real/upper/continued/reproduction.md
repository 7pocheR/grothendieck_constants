# Reproduction scope

This package reports an unpublished computer-assisted bound. Its executable
standard-library check verifies the endpoint consequences of supplied
numerical enclosures. It is not a complete portable numerical reproduction
of those enclosures.

## The check available from this package

Run `python3 -B verify.py`. The same command is supplied in `commands.txt`.
It checks the immutable file inventory and exact physical candidate,
all 404 coefficient intervals, all 404 frequency-choice indices, the
complete scalar-error subtraction, the localized two-square-root bound,
the strict objective margin, and the rational endpoint comparison using
Machin's identity. It also checks the completeness of the names in the
external data inventory. It does not read or validate the contents of
those absent external files.

Every coefficient interval is a pair of exact rational endpoints. The
coefficient file contains the finite Fourier approximation before the
separate analytic allowance. The scalar file contains exact rational
upper bounds for that allowance and both localized squared norms. The
checker recomputes the absolute-value objective from all coefficients;
it does not trust a stored decimal or a success flag. Hash checks establish
file identity, not the truth of the numerical enclosures.

## Numerical software and source modules

The numerical calculation used Python 3.11.14, python-flint 0.8.0 and
FLINT 3.3.1 on Linux. Exact integers and rational values, Arb real balls
and Acb complex balls are used. The native libraries and extension hashes
are listed in `certificate/software.json`. Python 3.11 or later suffices
for the standard-library checker; that checker does not need python-flint.
For numerical reproduction, keep assertions enabled and set the power-series
capacity explicitly to 404. A constructor's requested precision alone
does not override FLINT's global series cap.

The following source modules are included without numerical changes:

- `source/head/primitives.py` evaluates complete direct spatial Hermite
  vectors, shifted primary profiles, degree-403 pair products and both
  zero corrections. Its `exact.py` dependency is alongside it.
- `source/head/error_replay.py` evaluates every analytic head error at a
  requested precision, using the complete supplied frequency choices.
- `source/tail/directional_density.py` evaluates the full Gaussian-density
  derivative norm through a product Laguerre basis.
- `source/tail/cutoff.py` computes the exact rational cutoff derivatives
  and their Bernstein bounds.
- The remaining tail modules give the conditional Fourier terms, the
  mixed-energy calculation, its error estimates, and the whole-circle
  derivative calculation. They retain their original input-file interface.

These mathematical modules have a portable arithmetic interface. Some
original command-line entrypoints require manifest and data files absent
from this small package. Do not interpret their presence as a tested
one-command full reproduction system. No installation, allocation, process
supervisor or remote access program is supplied here.

## Complete finite-head reproduction

The following is a finite reproduction specification using the included
source functions. It requires a driver with explicit output and resource
management. A new complete execution of this specification has not been
performed as part of this portable package.

1. Load `certificate/physical.json` and the parameters in
   `certificate/head_error_input.json`. Set `ctx.prec=512`, `ctx.cap=404`
   and `ctx.threads=1`; `primitives.configure()` checks this setting. Form
   the exact Laurent rows with `exact.laurent(data,31,fmpz_poly)`.
   A call with a FLINT integer polynomial backend changes only the exact
   polynomial arithmetic. Each row is defined independently by the
   corresponding power, not by a supplied numerical profile table.
2. For every frequency index 1 through 768, calculate all 404 values from
   `primitives.direct_features(data,parameters,nodes)`. Its batches contain
   at most eight frequency indices. The 10,241 signed spatial nodes are
   included by parity, with separate treatment of zero. Save every vector
   with outward endpoints. Calculate all 404 primary values with
   `primitives.profile(n*Delta,laurent_rows)` for every frequency. Save the
   zero vectors from `primitives.zero_vectors` separately. All values
   contain exactly one factorial normalization.
3. Form the exact rational polynomial `P` as a checked ball series. For
   every pair `1 <= u <= v <= 768`, call `primitives.pair` with the full
   two characteristic and two primary vectors. Retain all 404 coordinates,
   weight diagonals by one and other pairs by two, and sum signed terms.
   There are exactly 295,296 pairs. Multiplication, composition, addition
   and scaling must all preserve capacity 404. Serialize outward; a
   continuation must retain all old interval enlargements and include
   each committed pair exactly once.
4. Multiply the complete sum by `4*Delta**2`, add
   `primitives.corrections(p,v,z,t,Delta,P)`, and project onto the odd
   powers. The correction is not multiplied by `4*Delta**2` again. A
   separate implementation can sum row endpoints exactly and recompute
   the zero corrections at 640 bits. Agreement or overlap between two
   valid interval calculations is a consistency test, not permission to
   replace them by an unjustified intersection.
5. Call `error_replay.replay(data,plan,768,bits=640)` to recompute all
   analytic bounds. Its input `plan` is exactly the object in
   `certificate/head_error_input.json`. This checks all 404 contour
   choices and the declared complete allowance. The planning target
   stored in the parameters is not a theorem hypothesis and need not
   be attained.
6. Use every resulting odd coefficient through 403 in equation (13) of
   the note, subtract the complete analytic error once, and combine with
   a separately valid infinite tail. Fresh interval endpoints may differ
   across numerical runtimes. Prove the final inequalities for the fresh
   result rather than importing this package's last digits.

The final saved exact aggregation was independently recomputed, but all
295,296 transient pair products were not independently rerun. Full numerical
reproduction must repeat them from newly computed or otherwise proved
primitive enclosures. A sum of saved rows alone does not perform that step.

## Complete tail reproduction

The original mixed-energy route has these finite components: 5,625 signed
conditional terms, 981 signed scalar files, 126 energy enclosures, and
2,048 complete first-quadrant panels with 14 mixed derivative enclosures
per panel. The full physical model is the same rational input as the head.
The conditional Fourier coefficients are assembled before taking the
relevant squared Gaussian norms. Every phase, lattice and omitted-frequency
allowance from the energy error budget is retained. Nonnegative radial
series give the interpolation between the energy-grid radii. The
whole-circle calculation combines these data with the complete chain rule.

The included original tail entrypoints require the external
`fourier_circle_manifest.json`, `fourier_circle_error_budget_final.json`,
the mixed-energy grid and its scalar records. The dominant-energy
entrypoint also uses `fourier_pilot_manifest.json` and
`fourier_error_budget.json`. These execution input files are not included
in this small bundle. The final circle input is a 6,104,017-byte file
identified by SHA-256 in `certificate/external_dependencies.json`.
The complete mixed-energy proof data must be restored with their original
bindings or regenerated; an identity hash does not supply the data.

For the localized tail, recompute each of the 2,048 intervals
`[i*pi/4096,(i+1)*pi/4096]`. Use the mixed estimates near real auxiliary
correlations, and use the complete six-order density function at 512 bits
only on intervals satisfying its strict denominator conditions. The
derivative with respect to the auxiliary correlation is evaluated with
`b=-z`; retain all four auxiliary pairs. Take the minimum only of complete
valid bounds. All fifth and sixth derivatives needed by a nonzero cutoff
term must be supplied. Set `M0=pi/2` from the Hermite norm bound.

For cutoff indices `(592,1248)`, use the exact Bernstein bounds from
`cutoff.py`, multiply their derivative order `k` by
`(4096/((1248-592)*pi))**k`, and retain every term in equation (18).
Sum the squared complete-panel upper endpoints exactly and divide by
2,048. This gives fresh `Q4,Q6` for equation (20). Their final square-root
comparison can be checked using integer/rational arithmetic, but producing
the bounds requires the physical mixed-energy and density calculations.

## Bulk data and remaining portability work

`certificate/external_dependencies.json` lists the hashes of all 1,536
primitive vector files, 768 signed row files and 2,048 localized panel
files. These records are absent from this bundle. The larger mixed-energy
grid, circle file, source manifest and error budget are also identified.
The transient pair results are not separately archived as 295,296 raw
files; their completed signed row sums are the retained consequences.

This inventory makes the omitted dependencies explicit. It is not a
download mechanism or a claim that an external reader can execute all
original programs from this directory. A complete portable numerical
release still needs either these data with a tested verifier, or a tested
driver that regenerates every required component on allocated hardware.
The included source modules and exact physical input provide a concrete
starting point for that work. No runtime, memory or successful-completion
guarantee for a full new calculation is made by the endpoint checker.
