# Grothendieck constants

Proofs, numerical certificates, and verification code for lower and upper bounds on the classical real and complex Grothendieck constants.

These are **unpublished computer-assisted research results**. The accompanying arguments and certificates are provided for independent examination. Neither exact constant is determined.

## Bounds in this snapshot

The definitions allow every finite matrix size and every finite vector dimension. Scalar variables are real signs for the real constant and unit complex phases for the complex constant.

- **Real lower:** $K_G^{\mathbb R}\ge 1375\pi/2454=1.760264832390369\ldots$. [Proof and verification](real/lower/quintic/).
- **Real upper:** $K_G^{\mathbb R}\le \pi/[2(0.8818276493988)]<1.781296297373$. [Paper and verification](real/upper/).
- **Complex lower:** $K_G^{\mathbb C}>1.373$. [Paper and verification](complex/lower/).
- **Complex upper:** $K_G^{\mathbb C}\le 100000000/71185999<1.404770620695$. [Proof and verification](complex/upper/).

The [two-page comparison](docs/grothendieck_bounds_comparison.pdf) gives exact statements, literature comparisons, numerical improvements, and comments on the methods. This snapshot contains the four results listed above.

## Comparison with the literature

For the real constant, [Saha et al., arXiv:2608.11158v2](https://arxiv.org/html/2608.11158v2) prove the lower bound $6\pi/11$. Their certified estimates in Proposition 6.3, together with Lemma 6.2 and Proposition 6.1, imply an upper bound of $1.78184132423347648\ldots$, slightly stronger than their displayed Theorem 6.4. These are the formal real comparison endpoints used here.

For the complex constant, [Guo, Fang and Lu, arXiv:2609.07000v1](https://arxiv.org/html/2609.07000v1) prove $K_G^{\mathbb C}>1.35584631827168$. Their result improves the classical lower bound of approximately $1.33807$ attributed to Davie. [Haagerup's upper bound](https://doi.org/10.1007/BF02790792) is approximately $1.404909132736$.

Guo, Fang and Lu also bound the optimum of their common radial-weight Gram/Schur criterion below $1.35584697425050$. This restricts that criterion, not the Grothendieck constant or the full Gaussian multiplier family. The complex lower proof here uses different weights for real and imaginary components and retains the exact joint midpoint phase constraint. The Gaussian projection framework and the proposed cubic perturbation have prior antecedents, including [Heilman, Section 1.8](https://arxiv.org/html/2603.22616v1).

The companion [Li et al., arXiv:2608.11195v3, Section 3](https://arxiv.org/html/2608.11195v3) reports stronger internally tested real claims, including $51\pi/92$ and $1.7813319810625639$, whose certificates the authors have not independently verified and which they do not state as theorems. Our lower endpoint exceeds that reported value by approximately $0.01872977$, and our upper endpoint is approximately $0.00003568$ smaller. The distinction in verification status is retained in the comparison document.

The quintic construction improves this repository's earlier lower bound $1625\pi/2917$ by approximately $0.01014894$. With the real upper endpoint fixed, it removes **32.55%** of the former remaining real interval. The earlier [cubic construction](real/lower/) remains available.

## Methods

- The real lower proof bounds $(33/20)P_1-P_3-(3/20)P_5$, where $P_j$ is the Hermite projection of degree $j$. Full three-by-three transverse majorants, agreement estimates, a circle inequality, and one-dimensional ternary certificates give a scalar bound in every dimension. Finite Gaussian partitions give explicit matrix witnesses.
- The real upper proof uses polynomial Gaussian threshold rounding and signed tensor preprocessing. Certified finite coefficients and an analytic bound on the entire remaining coefficient tail give a universal inequality for finite matrices.
- The complex lower proof combines a cubic Gaussian correction with anisotropic radial multiplication inequalities and exact midpoint orthogonality. A finite family of auxiliary weights covers the complete scalar parameter range for one fixed operator.
- The complex upper proof uses nonlinear complex phase rounding, including a radially damped perturbation, and signed tensor preprocessing. The certificate includes all omitted scalar degrees and phase orders.

## Reproduction and interpretation

Each result directory contains its mathematical statement and reproduction instructions. The complex lower certificate uses Python's standard library; the other numerical packages use exact rational arithmetic and/or FLINT/Arb through `python-flint`. Large integrations are separate from checking supplied certificate data. Package instructions distinguish these operations and state any platform requirements.

Passing a numerical verifier establishes its stated finite or interval-arithmetic assertions. The dimension-independent reductions, infinite remainder estimates, and transfer to the classical constants are mathematical parts of the accompanying proofs and must also be assessed.

The certificate data and proof sources are retained with their original numerical precision. Existing proof documents may discuss earlier comparison bounds; the literature comparison above and the two-page document describe the baseline for this snapshot. With Python 3.9 or later, run `python3 verify_files.py` to check the files listed in `SHA256SUMS`; this integrity check is separate from mathematical verification. No proof-assistant formalization is claimed.
