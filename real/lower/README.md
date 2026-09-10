# The cubic lower-bound construction for the real Grothendieck constant

The stronger [quintic construction](quintic/) proves
$K_G^{\mathbb R}\ge1375\pi/2454=1.760264832390369\ldots$.
This directory preserves the earlier cubic proof of

$$
K_G^{\mathbb R}\ge\frac{1625\pi}{2917}>\frac74.
$$

[Mathematical proof](proof.md) gives the analytic argument, the finite
certificates, and an explicit finite matrix whose ratio is greater than
$7/4$. The matrix is specified by a formula; its entries are not enumerated.
The proof uses an analytic circle inequality and two inequalities for ternary
functions. [Additional noise argument](noise-proof.md) gives an additional
Gaussian noise argument for one part of the same bound.

Original Markdown sources: [main proof](MATHEMATICS.md) and
[noise argument](NOISE_MATHEMATICS.md).

The mathematical certificate data are in `certificates_02/`. The original
subdivision inputs for the verification programs are in `inputs/fibers/` and
`inputs/scalar_campaign/`. The latter preserves the relative directory layout
required by `check_scalars.py`; it has no external directory dependency.

## Dependencies

The exact rational checks require CPython and its standard library. Fresh
integral verification also requires `python-flint`, whose FLINT/Arb arithmetic
provides rigorous real-ball enclosures. CPython 3.11 and `python-flint` 0.8.0
are the runtime versions used for the supplied certificate values.
`requirements.txt` specifies that package version.

Assertions are part of the verification. Do not use Python's `-O` or `-OO`
options or set `PYTHONOPTIMIZE`. Run the commands below from this directory.
The programs use explicit paths for inputs and outputs and require no network
access or cluster scheduler.

## Exact arithmetic checks

    python3 -B check_final_exact_02.py
    python3 -B check_certificate_arithmetic_02.py
    python3 -B check_fibers.py inventory --evidence inputs/fibers --output verification/fibers

The first command checks the polynomial identities, matrix determinants,
rational bounds for pi, and finite-matrix margins. The second checks all
extracted subdivision and scalar partitions, rational sums, and comparison
margins. The third checks the original ternary subdivision inputs, including
complete coverage and their vertex sets. These commands do not recompute
any vertex or scalar-cell integral.

The first two commands write `final_exact_02.json` and
`certificate_arithmetic_02.json` beside their scripts. The inventory command
writes to its explicit output directory. `final_arithmetic_02/scalar.json`
is a supplied input to the rational assembly check; it is not a substitute
for fresh integral verification.

## Fresh verification of all finite integrals

    python3 -B run_fiber_replay.py --evidence inputs/fibers --output verification/fibers_fresh --workers 2 --seconds-per-part 1200
    python3 -B run_scalar_replay.py --campaign inputs/scalar_campaign --output verification/scalars_fresh --workers 2 --seconds-per-part 1200

Use new output directories for each run. The first command recomputes all
5,638 vertex bounds for the two ternary certificates and checks all 1,966
leaves. It is the finite integral computation required by the main proof.
The second recomputes all 122,415 scalar-cell bounds, three Gaussian tail
bounds, and 3,072 quadrature panels used by the additional noise argument.
Both commands perform exact coverage and final comparison checks after the
individual computations. A successful run exits with status zero and writes
`verified.json` or `scalar_verified.json`, respectively, in its output
directory. A directory containing only some parts is incomplete.

The per-part time limit bounds each worker's computation. If any part is
incomplete, the final checks fail; no incomplete collection establishes the
result. The `--workers` option permits between one and sixteen worker
processes. Runtime depends on the machine and arithmetic libraries.

## Scope of the verification

The proof quantifies over all measurable Gaussian sign pairs in every finite
dimension. The analytic arguments reduce these quantifiers to explicit
finite certificates. The programs check complete finite domains, retain
unresolved integration remainders, and include the full Gaussian tails.
The result relies on CPython's exact integer arithmetic and FLINT/Arb's
enclosure guarantees; it is not a proof-assistant formalization.
