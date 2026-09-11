# A certified upper bound for the complex Grothendieck constant

The [stronger construction](winding/) proves $K_G^{\mathbb C}\le 100000000/71188883<1.404713710707$. The earlier construction below and its complete verification package remain available.

For complex scalar phases and complex unit vectors, the construction proves

$$
K_G^{\mathbb C}\leq\frac{100000000}{71185999}
<1.404770620695.
$$

The [mathematical proof](round3/complex_radial_v4_independent_audit/proof.md)
derives the Gaussian phase kernel, its realization for arbitrary complex
correlations, and a sufficient coefficient inequality. The calculation uses
two sine phase perturbations, one with radial damping, and signed tensor
preprocessing. All omitted phase orders and scalar degrees are included in
the estimates.

The [original Markdown source](round3/complex_radial_v4_independent_audit/proof_final.md) is also available.

The [finite certificate](round3/complex_radial_v4_independent_audit/certificate_final.json)
contains 208 retained modes and 481 exact dyadic coefficient intervals per
mode. Its relative links identify all 208 compressed files containing the
complete circle and annulus subdivision trees. Those trees contain
1,703,936 boundary leaves and 1,832,072 annulus leaves.

## Check the supplied data and exact final comparison

Run from this directory with Python 3.11 or later:

    python3 -B verify_certificate.py

This uses only the Python standard library. It checks every file listed in
`SHA256SUMS.json`, all 100,048 coefficient records, their exact signed
aggregation, the complete rational scalar tail, the omitted phase orders,
and the two omitted primitive terms. In particular it verifies

$$
\Gamma>0.711859994603>0.71185999.
$$

The supplied leaf files are hashed at this stage. Their subdivision trees
and local transcendental enclosures are checked separately below. A `PASS`
from this first command therefore records integrity and exact aggregation;
it is not a fresh evaluation of every local interval bound.

The frozen candidate's `target_gamma` field is the weaker value
`0.7118598`. The final exact aggregation and the complete rational-tail
checker explicitly establish the theorem's stronger coefficient
`71185999/100000000`.

## Validate every saved subdivision tree

The full validation commands run on Linux. The frozen scripts require the
literal `--cluster-resident-handoff` flag and reject macOS. The flag enables
computation on the machine running the command; the scripts use its CPU and
filesystem without SSH, a scheduler, or a cluster account.

Set the following shell variables from this directory. Use a fresh output
directory for each repetition; the scripts refuse to overwrite results.

    complex_audit=round3/complex_radial_v4_independent_audit
    complex_records=cluster/resident_20260907/results/complex_upper_radial_independent/cluster_attempt/independent_512_480_8192
    complex_verification=verification
    mkdir -p "$complex_verification"

    python3 -B "$complex_audit/check_validator.py"
    python3 -B "$complex_audit/validate_replay.py" \
      --cluster-resident-handoff \
      --records "$complex_records" \
      --candidate "$complex_audit/source_snapshot/radial_low_candidate.json" \
      --formulas "$complex_audit/primitive_polynomials.json" \
      --output "$complex_verification/saved_validation.json"
    python3 -B "$complex_audit/check_rational_tails_v2.py" \
      --candidate "$complex_audit/source_snapshot/radial_low_candidate.json" \
      --formulas "$complex_audit/primitive_polynomials.json" \
      --records "$complex_records" \
      --validation "$complex_verification/saved_validation.json" \
      --output "$complex_verification/saved_rational_tails.json"

These commands use standard-library rational arithmetic. The validator
checks every root interval and subdivision, complete coverage, leaf hashes,
the positive dyadic witnesses, squared boundary inequalities, and the full
Cauchy tails. It does not reimplement Arb's transcendental interval
operations. The final tail checker independently constructs a rational
boundary bound from the determinant polynomials and the positive per-mode
minima, and checks the stronger coefficient required by the theorem.

## Recompute all local interval bounds and coefficients

The complete interval computation was run with Python 3.11.14 and
`python-flint==0.8.0`. The exact symbolic check uses `sympy==1.14.0`.
These two Python packages are listed in `requirements.txt`. Install them
into a Python environment before running this section.

First reconstruct the determinant polynomials and phase weights:

    python3 -B "$complex_audit/check_identities.py" \
      --source "$complex_audit/source_snapshot" \
      --output "$complex_verification/identities"
    cmp "$complex_verification/identities/primitive_polynomials.json" \
        "$complex_audit/primitive_polynomials.json"

Then recompute all retained coefficients and all covering rectangles using
512-bit Arb arithmetic, degree parameter 480, and 8,192 initial angular
panels:

    python3 -B "$complex_audit/replay_independent.py" \
      --cluster-resident-handoff \
      --candidate "$complex_audit/source_snapshot/radial_low_candidate.json" \
      --formulas "$complex_audit/primitive_polynomials.json" \
      --output "$complex_verification/recomputed" \
      --workers 4 --bits 512 --degree 480 --panels 8192 \
      --radius 103/100 --max-leaves 200000
    python3 -B "$complex_audit/validate_replay.py" \
      --cluster-resident-handoff \
      --records "$complex_verification/recomputed" \
      --candidate "$complex_audit/source_snapshot/radial_low_candidate.json" \
      --formulas "$complex_audit/primitive_polynomials.json" \
      --output "$complex_verification/recomputed_validation.json"
    python3 -B "$complex_audit/check_rational_tails_v2.py" \
      --candidate "$complex_audit/source_snapshot/radial_low_candidate.json" \
      --formulas "$complex_audit/primitive_polynomials.json" \
      --records "$complex_verification/recomputed" \
      --validation "$complex_verification/recomputed_validation.json" \
      --output "$complex_verification/recomputed_rational_tails.json"

The last command must finish successfully and report a coefficient lower
bound strictly exceeding `71185999/100000000`. All preceding commands
must also finish successfully. A completion file from the interval producer
alone is insufficient.

The recorded calculation took approximately 229 seconds for the interval
producer with four workers and 215 seconds for the subsequent exact tree
validation on Linux x86-64. Hardware and filesystem differences affect
these times. Peak memory was not recorded. Each validation streams the
compressed leaves; their combined uncompressed size is approximately
707 MB. A full recomputation writes approximately 186 MB of new data.

The computation depends on Python's exact integers and rational arithmetic
and on Arb's outward interval arithmetic guarantees. It is not a formal
proof-assistant development. The complete mathematical argument and all
data used by these commands are included.

`source_snapshot/certify_radial.py` is used by the symbolic identity check
to compare the compact formulas with independently expanded determinants.
The full interval producer is `replay_independent.py`. The mathematical
proof for this package is `proof_final.md`.
