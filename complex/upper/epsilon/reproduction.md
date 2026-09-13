# Reproduction scope

The command `python3 -B verify.py` is complete for the exact consequences
of the six supplied signed groups. It uses Python 3.11 or later and no
native numerical library. It writes no files by default; `commands.txt`
also gives an option to save a new report outside the package.

The checker verifies all package hashes, the fixed baseline candidate,
admissibility of all three preprocessing polynomials, the 41 ordered
weight-ratio identities, all 3,381 saved group intervals, the complete
finite objective at both epsilon values, every scalar and phase error,
the transported complete omission bound, and the exact strict reciprocal
comparison. The coefficient grid is 2^-320. A stored integer pair `(a,b)`
means the exact interval `[a/2^320,b/2^320]`; a stored tail integer is
an upper bound divided by 2^320. Radial index n means the coefficient of
`t*|t|^(2n)`, not ordinary degree n.

These checks do not prove the supplied group enclosures from physical
integrals. The six groups are consequences of the original complete
9,948-mode calculation. Each group includes every coefficient of each
assigned retained mode, giving 2,415,418 contributions in total. Their
78 original blocks and the primitive coefficient and contour files are
outside this small bundle. The external inventory records all 78 block
hashes and exact half-open index ranges, plus identities of the original
retained inventory and the 39,792 mode/circle dependency records. Checking
that these names partition the modes does not validate absent contents.

## Reconstruction from the original coefficient records

For each of the 9,948 original retained modes, validate its exact physical
key, baseline signed weight, complete coefficient interval list and whole
scalar tail against the original coefficient and contour certificate.
Multiply the coefficient interval by the signed weight, swapping endpoint
roles for negative weights. Round its lower endpoint down and upper
endpoint up to the 2^-320 grid. Round the absolute weight times its entire
scalar tail upward to that grid. Assign the contribution to the class
determined by the two new-phase frequency magnitudes. Sum these integers
within each class before taking coefficient absolute values. The result
must contain all retained modes and agree with the six supplied arrays,
or independently produce valid new enclosing arrays.

An independent full grouping verifier must check all original dependencies
and each completed block's relation to its own mode interval. Reading only
the six final arrays cannot establish that relation. The original accepted
grouping and its separate exact block summation remain execution premises
for this small package. Neither is described as a fresh evaluation of
every original native coefficient.

The full finite expansion has 120,470 canonical modes. To reconstruct the
omission sum, preserve the same 9,948 retained keys and sum the entire norm
bound in equation (14) of the note over all 110,522 other keys, using the
original epsilon 1/100. Reweight this complete original total by the
largest of the six positive ratios. Reapplying a selection threshold at
the new epsilon without checking the changed support would define a
different calculation.

## Numerical sources and the preceding coupled package

The included `source/coupled_exact.py`,
`source/source_snapshot/winding_exact_core.py`, and `source/coupled_arb.py`
are unchanged mathematical implementations of the six-coordinate
determinants, full phase expansion, outward scalar coefficients and
complete contour estimates. The sibling `candidate_definition.json`
is the unchanged baseline candidate because the primitives reused in
this calculation were generated at that candidate. The selected epsilon
is applied by the exact weight ratios in `verify.py`; it must not also
be inserted into the old weights before reweighting them.

Numerical generation used Python 3.11.14, python-flint 0.8.0 and FLINT
3.3.1. Keep assertions enabled and use the original explicit precision,
degree, zero-count and complete-circle conditions. The original
[coupled construction](../coupled/README.md) supplies the longer mathematical
proof and a portable generation and validation interface. That package
remains a separate result and is not overwritten by this amplitude
refinement. Its reference raw data also have the scope stated there.

A full independent numerical reproduction would generate and validate
every retained coefficient and complete contour, reconstruct the full
omission sum, form the six signed groups from those new enclosures, and
then apply this note's exact reweighting. Fresh intervals can differ from
the supplied ones; prove the endpoint using the fresh bounds rather than
requiring the old last digits by assumption. The current small package
does not include a tested adapter joining those full stages, and its
endpoint command performs none of the native calculations.

The amplitude change requires no new primitive integrals mathematically,
because all functions in the finite expansion are epsilon-independent.
This reuse does not eliminate the original numerical proof obligations.
The coefficient and omission enclosures are explicit premises; file hashes
establish identity rather than validity.
