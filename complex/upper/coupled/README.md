# A coupled complex Grothendieck upper bound

The [amplitude refinement](../epsilon/) gives the stronger endpoint $K_G^{\mathbb C}<1.404698554831$ using the same primitive functions. This directory retains the preceding coupled construction.

For arbitrary finite complex matrices and arbitrary complex unit vectors,

\[
K_G^{\mathbb C}\le\frac{200000000000}{142379302169}<1.404698554869.
\]

The [proof](proof.md) gives the scalar and vector normalization, complex tensor
construction, Gaussian phase integral, damping, complete phase and omission
errors, and the contour argument controlling every scalar tail. The displayed
rational endpoint follows by exact integer comparison with the complete
coefficient gap, not from rounding a numerical estimate.

The complete reference calculation retained 9,948 of 120,470 nonzero canonical
modes, charged all 110,522 omissions, validated every retained contour, and
independently aggregated all 2,415,418 coefficient contributions. The complete
9,948-mode rerun using these portable entrypoints has **not** executed.

## Quick metadata and arithmetic check

With Python 3.11 or later, run from this directory:

    python3 -B check_exact.py

This standard-library command checks package integrity, the rational candidate,
reference completion metadata, retained identities and costs, the complete
raw-file hash inventory, and the six exact rational terms and final comparison.
It performs no physical evaluations, does not visit raw contours, and does not
reconstruct the full 120,470-mode expansion. The latter has a separate exact
command:

    python3 -B check_inventory.py

The supplied metadata and hashes do not replace correct outward computation or
complete raw validation. The source needed for both is included.

## Generate, replay, validate and aggregate

Use Linux, Python 3.11 or later, and `python-flint==0.8.0`. The exact-only
commands need no third-party library. Install the dependency into an existing
environment, or use the commands in [commands.txt](commands.txt). Choose an
explicit external output directory; the programs reject output inside the
package. In the following commands, set `coupled_output` to that directory.

    python3 -B run_generation.py --output "$coupled_output/raw" --stage conditioning --workers 4
    python3 -B run_generation.py --output "$coupled_output/raw" --stage complete --workers 4 --resume
    python3 -B run_validation.py --evidence "$coupled_output/raw" --output "$coupled_output/validation" --stage conditioning --cpus 4
    python3 -B run_validation.py --evidence "$coupled_output/raw" --output "$coupled_output/validation" --stage complete --cpus 4
    python3 -B independent_aggregate.py --raw "$coupled_output/raw" --validation "$coupled_output/validation" --output "$coupled_output/aggregation" --stage complete
    python3 -B check_exact.py --validation "$coupled_output/validation" --aggregation "$coupled_output/aggregation"

Every command must finish successfully. Generation includes both production
and fresh outward replay; it cannot issue a theorem result. The complete exact
validator writes `validation/complete.json` only on complete success. The
independent aggregation writes `aggregation/complete.json` only after every
required mode and every exact contribution agrees with that validation. The
last command requires the complete executions and a strict rational margin
above `142379302169/200000000000`.

Numerical conditioning uses a deterministic spectrum of at most 26 retained
modes at their full tail-driven degrees. Contour-validation conditioning uses
12 fixed retained indices and requires the complete raw tree to be available.
Neither conditioning result is sufficient for the upper bound. An optional
generation `--indices` selection and aggregation `--stage pilot` are diagnostic
only; the complete checks still require all 9,948 retained modes.

Generation uses the supplied reference schedule of successful radii, degrees
and precisions. It regenerates each contour and coefficient with the unchanged
mathematical kernels and checks the same tail and width budgets. It does not
spend time reproducing unsuccessful radius attempts. Failure at a prescribed
parameter leaves that mode required. The fresh replay still uses eight more
coefficient indices and 128 more precision bits than production.

## Resume or independently repeat

Repeat an interrupted generation command with `--resume`. Repeat an interrupted
validation command on the same validation directory. Repeat an interrupted
aggregation command with `--resume`. Do not edit the package or raw inputs
during execution. Source, rational input, environment and raw dependency
changes invalidate the affected continuation. Validation uses the complete
specified source and verifies the production-to-replay links as well as each
coefficient, contour and mode identity. A completed mode is never replaced by
an omission. Incomplete generation attempts are preserved for inspection.

Validation checkpoints are records of a prior complete execution of the
specified exact contour checker. On reuse their raw hashes and reconstructed
arithmetic descriptors must still agree. To repeat every contour check, choose
a new empty validation directory. To independently decode every coefficient
again, choose a new empty aggregation directory. A partial report, a progress
count, or a saved success flag cannot satisfy complete acceptance.

All work directories are locked against concurrent writers. Numerical workers
and validation subprocesses use Linux time and memory limits and parent-death
protection. Run without `-O` and with `PYTHONOPTIMIZE` unset. These programs use
the local CPU and filesystem; they do not use a scheduler, remote connection,
or account-specific command.

## Validate supplied original raw records

The approximately 5 GB production and replay tree is not stored in Git. There
is no assigned download URL. Fresh generation above needs no external evidence
download. If the original raw files are available, set `coupled_raw` to a
readable directory with `production/mode_000000/...`, `replay/mode_000000/...`,
and the two input metadata files. Actual labels are those in the complete
inventory, not a contiguous range of 9,948 integers.

`data/reference_raw_sha256.json.gz` lists all 39,792 required paths and SHA-256
hashes. `data/reference_execution_inputs.json` is the original execution input
file. `data/inventory.json.gz` preserves the original `complete_inventory.json`
bytes, as bound by `data/reference_binding.json`. If the two input metadata
files are missing, `python3 -B prepare_reference_inputs.py --raw "$coupled_raw"`
materializes them without replacing existing files. File byte sizes are not
asserted by the hash manifest. The references are metadata, not raw enclosures.

    python3 -B run_validation.py --reference --evidence "$coupled_raw" --output "$coupled_output/saved_validation" --stage conditioning --cpus 4
    python3 -B run_validation.py --reference --evidence "$coupled_raw" --output "$coupled_output/saved_validation" --stage complete --cpus 4
    python3 -B independent_aggregate.py --raw "$coupled_raw" --validation "$coupled_output/saved_validation" --output "$coupled_output/saved_aggregation" --stage complete
    python3 -B check_exact.py --validation "$coupled_output/saved_validation" --aggregation "$coupled_output/saved_aggregation"

This route binds all original raw files to the supplied reference hashes and
performs every exact contour implication. It does not rerun their outward
transcendental and formal-series evaluations. The generation route supplies
that separate obligation. The original raw tree is always read-only to these
validation and aggregation programs.

## Resources and files

Use fast local storage and allow at least 16 GiB for a new full work directory;
32 GiB is preferable when retaining repeated interrupted attempts. Numerical
generation allows one to eight workers, each with a 3,000,000,000-byte address
space limit. Four workers therefore need approximately 12 GB plus coordinator
and operating-system memory. Each mode has a 150-second limit; each radius
attempt has 20 seconds. The full generation invocation allows 10,320 seconds.
The 16 GiB generation disk limit stops further scheduling at periodic checks;
active bounded workers may still write, so it is not a strict total disk cap.

Exact validation allows one to four workers with the same per-worker memory
and time limits. Its separate output is capped at 1 GiB, apart from the small
final report copy, and a complete invocation allows less than two hours.
Independent aggregation uses one CPU, a 2 GiB address-space limit and 2,400
seconds per complete invocation. Its resumable integer contributions may use
approximately 0.5–1 GB. More than one invocation may be necessary on slower
machines; incomplete runs retain their required modes.

The complete reference contour validation took approximately 1,775 seconds,
and the independent aggregation took approximately 905 seconds. These timings
do not benchmark the portable entrypoints; hardware, filesystem and resume
verification costs affect them. Generation must first be timed with the
conditioning command. The reference native environment was Python 3.11.14 and
python-flint 0.8.0; further version information is recorded in the supplied
execution metadata.

- [proof.md](proof.md): the mathematical argument and complete error estimates.
- [candidate_definition.json](sources/candidate_definition.json) and
  [full_task.json](sources/full_task.json): exact candidate and fixed settings.
- `data/inventory.json.gz`: exact retained list, complete inventory digest and
  full omission cost; `check_inventory.py` reconstructs all modes twice.
- `data/reference_*.json*` and [endpoint.json](data/endpoint.json): complete
  reference metadata, schedule, raw dependency hashes and rational endpoint.
- `sources/`: determinant, Gaussian formal-series, complete-arc and exact
  coefficient routines; every mathematical dependency is included.
- `validation/`: exact contour checker, checkpoint management and Linux
  subprocess controls, including their executable conditioning checks.
- [independent_aggregate.py](independent_aggregate.py): independent binary-ball
  decoding and exact integer aggregation with resumable contributions.
- [check_exact.py](check_exact.py): strict metadata and final arithmetic checks.
- [SHA256SUMS.json](SHA256SUMS.json): package integrity; nested manifests bind
  the mathematical sources and exact validator separately.

The result depends on correct Python exact arithmetic and FLINT/Arb outward
arithmetic. It is not a proof-assistant formalization.
