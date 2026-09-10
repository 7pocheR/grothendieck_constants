Reproducibility of the real Grothendieck upper bound

The numerical certificate supports
  K_G^R <= pi/(2*0.8818276493988) < 1.781296297373 < 1.78130.
The mathematical proof is in ../paper.tex and ../paper.pdf. The certificate
alone does not replace the infinite-series and finite-matrix arguments.

Software

Use Python 3.11 or 3.12 with python-flint 0.8.0 and FLINT 3.3.1. The
certificate was also checked with Python 3.12. Exact rational calculations
use fractions.Fraction and FLINT fmpq. No numerical NumPy/SciPy operation
is needed for certification. SciPy, if present, supplies untrusted initial
root guesses; the alternative uses the included old node midpoints solely
as guesses. Every root and weight is subsequently certified again.
Do not run Python with -O: mathematical assertions must remain enabled.
No script installs packages, contacts a remote service, or submits a job.

Read-only verification of the supplied certificate

From this directory, run:
  python -B -u verify.py --output /tmp/real-upper-verification

Choose an output directory outside this package. New evidence goes there.
All 56 immutable input/program/certificate files are hashed. The verifier
checks the exact index sets, all 32 block hashes and all 128 panels, re-sums
all 2,024,064 partial moment intervals, and reproduces all 63252 merged
moment balls exactly. It independently recomputes the degree-501 uniform
quadrature error from the unchanged source and checks containment in the
initialized saved ball. It proves all 112 node and weight enclosures by
exact rational arithmetic, reruns both complete coefficient assemblies,
and recomputes every smoothing integral and every whole-circle interval.
The integral cache starts empty; supplied integral values are not used as
the input to the new integral calculation. Finally, rational arithmetic
checks the head, norm, tail, positive margin, Machin enclosure of pi, and
strict reciprocal comparison.

This command does not repeat the expensive evaluation of each moment at
each quadrature node. Complete repetition of that calculation is below.
Individual short checks are available with --stage integrity, merge, head,
norm, or final. The default --stage all runs all of them. A failure exits
nonzero and must be investigated; a missing file is not silently ignored.

Complete reproduction from new quadrature calculations

Use reproduce.py with a new work directory and computational resources
already allocated to you. This program preserves the numerical engines
byte for byte. It copies only code, exact input, and untrusted node guesses
to a new calculation; it does not copy moments, coefficient outputs,
integral caches, or circle outputs into that calculation.

  python -B -u reproduce.py full --workdir /scratch/real-upper-reproduction --workers 2

The default moment-engine guard requires an inherited SLURM allocation.
On other assigned hardware, explicitly add --allow-long-local-computation.
This flag passes the existing --probe flag to the unchanged block engine;
it changes the execution guard and metadata, not its numerical algorithm.
Choose the number of workers to fit the resources assigned to the process.
Each worker executes one moment block at a time. The default is one worker.

The phases are:
  prepare: copy and hash the unchanged programs and exact input;
  init: certify 112 new nodes and weights, then recompute the uniform error;
  block: calculate one of 32 four-panel moment blocks;
  merge: validate every block and add the complete uniform error once;
  head: run both complete degree-501 coefficient assemblies;
  norm: start an empty cache and compute all integrals and 512 circle panels;
  final: check exact rational bounds from the newly produced outputs.

The full command runs these phases in order. Each completed block is an
atomic compressed JSON checkpoint with a configuration hash. Rerunning
full reuses completed blocks only after matching their configuration.
To resume a specific block, use stage block with --block INDEX. Every
individual command and exit code is saved with a hashed output log.
Commands for each phase are also supplied in commands.txt.

Certificate layout

head/sources: unchanged exact input and original head/assembly sources.
head/checks: unchanged checkpoint engine, rational assembly, and exact nodes.
head/run: complete original 501 nodes, configuration, blocks, and moments.
head/evidence: complete rational assembly and exact finite-algebra outputs.
head/prior/gauss_nodes_112.json: guesses only, never accepted without proof.
norm/checks/independent.py: unchanged independent integral and circle engine.
norm/evidence: all 266 integrals, 512 circle rows, and exact Euler polynomials.
rational_certificate.json: exact bounds and all 512 rational panel caps.
immutable_manifest.json: SHA-256 inventory of all preceding immutable files.
verify.py, reproduce.py: portable verification and orchestration programs.

The norm engine retains its original N=251 setting and auxiliary tail251
output. Its C4 computation is independent of N: only Hermite polynomials
through degree 9, the degree-81 polynomial, and derivative order 4 enter it.
The portable final calculation always uses 14*501**7 for the tail.
There are 251 positive odd indices through 501, not 252. All are present.

Numerical dependence and interpretation

Every ball (m,r,e) means the exact real interval [(m-r)*10**e,(m+r)*10**e].
Reads, arithmetic, and re-serialization enlarge outward when necessary.
The stored success flags, source comments, and diagnostic candidate values
are not proof premises. The mathematical derivation and executed enclosing
arithmetic establish their validity. Hashes establish identity, not truth.

Adaptive integral enclosures can differ across supported interpreters or
subdivision decisions. The verifier requires valid new bounds proving
C4 < 24744.587741, and separately checks the exact rational sum Q in the
supplied certificate. New rational panel sums need not equal Q bit for bit.
Each complete reproduction writes its own exact rational norm sum and
final directional comparisons to final.json.

The standalone manuscript and the certificate specify all mathematical
dependencies. No historical priority or determination of the exact
Grothendieck constant is asserted.
