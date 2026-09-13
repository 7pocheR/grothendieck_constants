# A real Grothendieck upper bound from a five-variable sign function

**Unpublished computer-assisted research.** Our work gives

\[
K_G^{\mathbb R}
\le \frac{\pi}{2\Gamma},\qquad
\Gamma>\frac{17651832366373}{20000000000000},\qquad
K_G^{\mathbb R}<1.779754412112.
\]

The [mathematical note](mathematical_note.md) specifies the sign function,
normalization, finite coefficient calculation, complete infinite tail, and
implication for arbitrary finite real matrices. The [head estimates](head_estimates.md)
give the full analytic error formulas. This is an upper bound; neither an
exact determination of the constant nor a claim of historical priority is made.

Run `python3 -B verify.py` from this directory using Python 3.11 or later.
On POSIX, use `python3`; substitute the executable name for your Python 3.11-or-later installation if different. The command needs only the standard library and writes no files. It checks
the package hashes, all 404 supplied coefficient intervals, the single
analytic-error subtraction, the rational comparison for the sum of two
square roots in the tail, and a strict rational enclosure of pi. Commands
are also available in [commands.txt](commands.txt).

This small package does **not** provide a complete independent numerical
proof by running that command. The primitive Gaussian calculations,
295,296 transient pair products, and derivative bounds remain the specified
computer-assisted premises. Their retained file inventories are included;
the bulk primitive, row, mixed-energy, and circle-panel records are not.
The final exact checker does not independently rerun the transient pair
products or prove the primitive integrals.

[Reproduction scope and numerical sources](reproduction.md) describes the
included mathematical kernels, the additional data, and the finite steps
needed for reproduction. It distinguishes a proposed portable calculation
from the standard-library check that is executable from this directory.
No complete portable rerun of this new package is claimed.
