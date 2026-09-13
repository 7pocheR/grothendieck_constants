# An upper bound for the real Grothendieck constant

A stronger upper bound, $K_G^{\mathbb R}<1.779754412112$, is available in the [five-variable construction](continued/). This directory retains the earlier polynomial-threshold result.

The accompanying proof establishes

$$
K_G^{\mathbb R}\leq\frac{\pi}{2(0.8818276493988)}
<1.781296297373.
$$

Read the [manuscript PDF](paper.pdf) or the [LaTeX source](paper.tex) for the mathematical argument. The finite certificate accompanies the proof; its verification does not replace the infinite-series and finite-matrix arguments in the manuscript.

## Verify the supplied certificate

Use Python 3.11 or 3.12 with `python-flint==0.8.0` and FLINT 3.3.1. The Python dependency is listed in [requirements.txt](reproducibility/requirements.txt). Assertions must remain enabled; do not use Python's `-O` option.

From this directory, run:

```sh
python -B -u reproducibility/verify.py --output /tmp/real-upper-verification
```

Choose an output directory outside `reproducibility`. The command reads the supplied certificate and writes new verification evidence only to the output directory. It verifies all 56 immutable files, re-sums every supplied moment block, checks the node and weight enclosures, reruns both coefficient assemblies, recomputes all smoothing integrals and circle enclosures, and checks the final rational inequalities. A successful run ends with `VERIFICATION COMPLETE: all`; a failed check exits nonzero.

This command does not regenerate the moment blocks by evaluating the quadrature nodes. Complete reproduction instructions, software details, individual verification stages, and resource requirements are in [reproducibility/README.txt](reproducibility/README.txt) and [commands.txt](reproducibility/commands.txt). All source programs, exact input, supplied blocks, and numerical certificates required by those instructions are included.

## Build the manuscript

With a standard LaTeX installation, run `pdflatex paper.tex` twice from this directory. Its five local input files—`sobolev.tex`, `conclusion.tex`, `polynomial.tex`, `certificate.tex`, and `references.tex`—are included.

The result is an upper bound. No exact value or historical priority is asserted.
