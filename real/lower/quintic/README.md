# A quintic Gaussian lower bound for the real Grothendieck constant

This package proves

$$
K_G^{\mathbb R}\ge\frac{1375\pi}{2454}
=1.76026483239036907208\ldots.
$$

The definition permits every nonzero finite real matrix and every finite
real vector dimension. The scalar row signs and column signs are chosen
independently. The exact constant is not determined.

[Complete proof](proof.md) derives a universal scalar-norm bound for
$(33/20)P_1-P_3-(3/20)P_5$, where $P_j$ is Gaussian Hermite projection
of total degree $j$. The proof reduces all measurable sign pairs in every
finite dimension to three explicit six-dimensional inequalities, together
with exact agreement and scalar-envelope inequalities. Section 9 specifies
one finite matrix with ratio strictly greater than $44/25$. The limiting
bound above is not asserted to be attained by that single finite matrix.
The matrix is given by a formula; its entries are not enumerated.

The original [proof source](PROOF_01.md) is retained unchanged. Its explicit
finite hypothesis is established by the complete certificates distributed
here. [Integration mathematics](integration-mathematics.md) describes the
outward Gaussian calculation, all unresolved remainders, and the complete
infinite tails; its [original source](REPLAY_MATHEMATICS_01.md) is retained
unchanged. These are unpublished computer-assisted research results
provided for independent examination; no proof-assistant formalization is
claimed.

## Supplied finite certificates

The three certificates have the following complete counts and positive
minimum leaf margins:

- A: 54,542 vertices, 11,093 leaves, margin
  $6160026411051/112589990684262400$.
- B: 90,401 vertices, 17,460 leaves, margin
  $670337411673/112589990684262400$.
- C: 107,078 vertices, 18,846 leaves, margin
  $48409282727/281474976710656000$.

In total, the supplied outward execution covers all 252,021 distinct
vertices and all 47,399 terminal rectangles. Each rectangle has 64 corners.
All eighteen inequalities establishing the compact parameter domains are
included. No supplied fresh vertex bound exceeds its original proposed
bound. The scalar package checks 1,007 agreement supports and 200 endpoints
of the combined closed-interval subdivision.

`certificate_manifest.json` records source hashes and both compressed and
uncompressed data hashes. The six archives in `certificates/` preserve the
original mathematical inventory and outward execution records exactly.
They occupy about 47.3 MB compressed and 363.9 MB uncompressed. The complete
package includes additional scalar data and proof sources. Reproduction
also needs working space and memory for the expanded inventory.

`bindings_01.json` retains the historical hashes required by the unchanged
integration kernel. Historical producer files are not needed: the public
verifier reconstructs every mathematical rectangle and corner directly
from the inventory. No private filesystem path or cluster scheduler is
required.

## Dependencies and commands

Use CPython 3.11 or later. Exact scalar, geometry, and stored-result checks
use only its standard library. Fresh Gaussian integration additionally
requires `python-flint==0.8.0`. The supplied execution used CPython 3.11.14,
python-flint 0.8.0, and FLINT 3.3.1.

Run commands from this directory. [COMMANDS.txt](COMMANDS.txt) contains the
complete commands without Markdown indentation. Do not use `-O`, `-OO`, or
`PYTHONOPTIMIZE`. The public entry points reject optimization explicitly;
assertions are part of the unchanged integration and scalar programs.

1. `python3 -B verify_certificate.py --mode integrity --output verification/integrity.json` checks distributed
   source and compressed-file hashes. It makes no mathematical claim.
2. `python3 -B run_replay.py --stage scalar --work verification/scalar`
   independently checks rational exponential bounds, agreement supports,
   exact envelopes, the circle estimate, and finite-matrix arithmetic.
3. `python3 -B verify_certificate.py --mode stored` reconstructs all closed
   rectangles and corners using exact arithmetic, proves complete coverage,
   checks every supplied vertex record and its source binding, and assembles
   all domain and leaf margins. It computes no new Gaussian integrals.
4. `python3 -B run_replay.py --stage all --work verification/fresh --workers 4`
   performs a complete fresh outward Gaussian replay and its final exact
   validation. Use a new, empty working directory for an independent run.

A small smoke check is available with `--stage smoke`; it recomputes six
vertices and eighteen compact-domain column integrals. It is not a full
replay. A bounded replay invocation with `--rounds 1` retains checkpoints
and can be repeated with the same working directory. The partition count
is always four; `--workers 1` through `--workers 4` control concurrency.
Each worker invocation is bounded by time and by the number of newly
computed vertices. A full run may take substantial time, depending on the
machine and arithmetic libraries. Sending SIGINT or SIGTERM to the driver
terminates and reaps its active worker processes. Run only one driver per
working directory. If interruption leaves an incomplete final JSON record,
resumption fails; start in a new empty working directory.

Exit status 2 means that a bounded run is incomplete. Repeat the command
against the same working directory. An incomplete directory does not
establish the theorem. A successful complete replay returns zero and writes
`replay/verified.json`, `final_validation.json`, and `run_status.json` under
the chosen work directory; the status must state
`complete_universal_replay: true`. The driver rechecks actual record counts,
geometry, and inequalities before producing this status. Existing
checkpoints are reused rather than recomputed.

## What the checks establish

File integrity, exact assembly of supplied numerical outputs, and fresh
Gaussian integration are different operations. In particular, checking
arbitrary JSON endpoint values would not prove that those endpoints contain
the Gaussian integrals.

The distributed endpoint values were obtained by a complete execution of
the reviewed outward integration program whose exact source hash is bound
to every record. Its precision is 256 bits; reported endpoints are rounded
outward to a dyadic grid of denominator $2^{52}$. The unresolved-cell
tolerance is $2^{-48}$, and the reported upper bounds additionally include
all remaining unresolved cells and both infinite tails. A fresh replay
executes this same source at every required vertex.

The stored-result verifier checks complete geometric coverage, exact
source/task bindings, rational normalization and interpolation, all target
comparisons, and the retained domain-panel partitions. Its output explicitly
records zero new Gaussian integrals. Gaussian containment depends on correct
execution of the reviewed integration control flow and FLINT/Arb's outward
enclosure guarantees. The dimension-independent reductions and the transfer
to finite matrices remain mathematical parts of the accompanying proof.

## Comparison

This endpoint improves the earlier repository lower bound
$1625\pi/2917$ by approximately $0.010148938635342$.
It exceeds the higher reported endpoint $51\pi/92$ in
[Li et al., Table 2](https://arxiv.org/html/2608.11195v3) by approximately
$0.018729774422114$. Those authors label their Table 2 claims as system-tested
results not yet independently verified, rather than theorems. The formal
lower theorem $6\pi/11$ in
[Saha et al.](https://arxiv.org/html/2608.11158v2) is a separate comparator.
