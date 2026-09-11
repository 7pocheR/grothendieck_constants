# A certified upper bound for the complex Grothendieck constant

For arbitrary finite complex matrices and arbitrary complex unit vectors,

$$
K_G^{\mathbb C}\leq\frac{100000000}{71188883}
<1.404713710707.
$$

The [proof](proof.md) derives the Gaussian correlation, its complex tensor realization, the complete phase and omitted-mode errors, and the circle winding argument that controls every scalar tail. All parameters are exact rational numbers. No matrix size or vector dimension restriction is imposed.

The small Git bundle contains all mathematical inputs, all 792,700 coefficient intervals, the exact 2,204-mode schedule, and the sources needed to regenerate every interval and circle enclosure. The coefficient records occupy about 7.1 MB compressed. The optional complete native records occupy about 1.67 GB; they are described by an exact [per-file manifest](data/raw_records_manifest.json), and are not included in the small bundle. There is no assigned download URL in this package. A full fresh computation requires no external certificate download.

## Check the compact supplied data

From this directory, with Python 3.11 or later:

    python3 -B verify_certificate.py

This uses only the Python standard library. It checks package hashes, derives the Gaussian determinant polynomials, reconstructs all 2,962 nonzero phase modes in two different ways, checks the complete retained/omitted partition, and aggregates every supplied coefficient interval with exact outward rational rounding. It recomputes the complete phase error, all 758 omitted-mode costs, and the Cauchy tails implied by the supplied circle maxima. It must establish the strict rational comparison

$$
\Gamma>\frac{71188883}{100000000}.
$$

A `PASS` here establishes integrity and exact aggregation of the supplied enclosures and summaries. It does not establish that a supplied rectangle encloses the true function, and does not validate the individual saved arcs. Those obligations are addressed by the following two procedures. Hashes identify data; they are not mathematical enclosure proofs.

The coefficient encoding is `[lower_integer, width_integer]`, with both endpoints multiplied by $2^{-96}$. Its intervals contain the original 768-bit native enclosures. The signed aggregation uses 256-bit outward dyadic arithmetic. Every truncation error is charged independently of this coefficient compression.

## Validate all optional saved native records

If the complete records are available, place them in a directory with this layout:

    raw-records/mode_0014/mode.json
    raw-records/mode_0014/geometry.jsonl.gz
    ...

The manifest lists exactly two files for each retained mode, with their byte counts and SHA-256 hashes. It does not describe an existing release archive. To validate all supplied records:

    python3 -B verify_certificate.py --raw-records raw-records \
      --output verification/saved_validation.json

This command uses only the standard library. It checks every native record hash, verifies that every compact coefficient interval contains its original native interval, reconstructs every dyadic circle subdivision, and recomputes all rational arc inequalities, argument increments, minima, circle maxima and scalar tails. It visits all 2,305,576 supplied leaves and repeats the final aggregation.

This validates the exact implications of the recorded rectangles. It does not rerun Arb's polynomial, exponential, division or formal-series operations. The original complete native computation and complete exact validator were executed; their correctness uses the usual outward-arithmetic guarantee. To reproduce the native enclosure step independently, run the next section.

## Recompute every coefficient and complete circle

The full replay requires Linux, Python 3.11 or later, and `python-flint==0.8.0`. No scheduler, network connection, or account is used. The exact algebra and compact-data checks are portable. A small native smoke test is also available on platforms supported by python-flint.

    python3 -m venv .venv
    .venv/bin/python -m pip install -r requirements.txt
    .venv/bin/python -B smoke_test.py

The smoke test independently derives twelve coefficients by exact rational series and computes three small arc enclosures. It is deliberately partial and proves no complete circle bound.

For the full computation, choose a new output directory:

    .venv/bin/python -B run_replay.py --output verification/recomputed --workers 4
    .venv/bin/python -B verify_certificate.py --replay verification/recomputed \
      --output verification/recomputed_validation.json

Both commands must finish successfully, and the second must report the strict target comparison. The replay computes all 792,700 coefficients at 768-bit precision. Each of the 2,204 retained modes starts with 1,024 angular panels and refines until all required complete-arc and winding inequalities hold. The source and exact inputs determine every radius and scalar degree; no search or external producer data is used.

Up to four worker processes run concurrently, each limited to 2 GiB of address space and 600 seconds per mode. The driver enforces an additional outer process timeout, records failures, and terminates and reaps its children on interruption. Reserve at least 4 GB of disk for the replay and its checks; more may be needed if attempts fail and are retained. The supplied native records occupy about 1.67 GB. The previously completed full exact validation took about 3,861 seconds on Linux x86-64; compute and validation times depend on hardware. A numerical pass from the producer alone is insufficient.

Completed modes are reused only if the source, mathematical inputs, environment and record hashes are unchanged. To resume:

    .venv/bin/python -B run_replay.py --output verification/recomputed --workers 4 --resume

Incomplete attempts remain on disk, and a new attempt receives a new directory. A process interruption cannot promote an incomplete attempt to a complete checkpoint. Two drivers cannot operate on the same output directory concurrently. Optional `--indices` selects a partial computation for diagnosis; the complete verifier always requires every retained mode and cannot turn that selection into a theorem result.

Run these programs without `-O` and with `PYTHONOPTIMIZE` unset. All explicit mathematical guards use exceptions, and the public entrypoints also reject optimized Python execution. The sources do not constitute a proof-assistant formalization.

## Files

- [proof.md](proof.md): mathematical argument and all infinite-tail formulas.
- [candidate.json](data/candidate.json), [plan.json](data/plan.json), and [primitive_polynomials.json](data/primitive_polynomials.json): exact construction and complete schedule.
- `data/coefficients.jsonl.gz`: deterministic compressed coefficient intervals and circle summaries.
- [raw_records_manifest.json](data/raw_records_manifest.json): optional full native coefficient and geometry records, with no assigned download location.
- [verify_certificate.py](verify_certificate.py): compact aggregation and full saved-record/replay validation.
- [check_algebra.py](check_algebra.py): exact determinant and independent phase expansion.
- [run_replay.py](run_replay.py), [native.py](sources/native.py): resumable full Arb computation and mathematical interval routines.
- [exact_core.py](sources/exact_core.py), [geometry_check.py](sources/geometry_check.py), [package_core.py](sources/package_core.py): exact interval implications, full circle coverage and aggregation.
- [smoke_test.py](smoke_test.py): small independent series and native arithmetic checks.
- [commands.txt](commands.txt): the reproduction commands in plain text.

The driver uses relative paths and handles process limits and checkpoint management. The complete saved native records are optional because the small bundle contains everything needed for fresh execution.
