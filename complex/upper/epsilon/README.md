# An amplitude refinement of the coupled complex upper bound

**Unpublished computer-assisted research.** Our work gives
\[
 K_G^{\mathbb C}\le\frac1{\Gamma_*}<1.404698554831,
\]
where the exact rational \(\Gamma_*\) is in
`certificate/endpoint.json`. The selected amplitude is
\(\epsilon_*=1/100+2^{-28}=67108889/6710886400\).

The [mathematical note](mathematical_note.md) derives the exact reweighting
of six signed coefficient groups and includes all phase, scalar-tail and
omitted-mode errors. The linear mixing remains exactly \(-1/50\), and
the three preprocessing polynomials remain unchanged. This is a small
parameter refinement, not an exact-value or historical-priority claim.

Run `python3 -B verify.py` from this directory with Python 3.11 or later.
On POSIX, use `python3`; substitute the executable name for your Python 3.11-or-later installation if different. Only the standard library is needed. The checker recomputes both the
baseline and selected full signed objectives from 3,381 supplied group
intervals, checks all 41 ordered phase-weight identities, and verifies
the strict reciprocal comparison. Commands are in [commands.txt](commands.txt).

The groups were constructed from all 9,948 retained modes and 2,415,418
coefficient contributions. This small package does not supply or re-decode
all original coefficient records, reproduce their Gaussian primitives,
or reevaluate the complete contour certificates. Their validity and the
construction of the six groups remain explicit numerical execution premises.
[Reproduction scope](reproduction.md) describes the included source modules
and omitted data. The preceding coupled package remains a separate result.
