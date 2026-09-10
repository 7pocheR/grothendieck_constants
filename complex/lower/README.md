# A lower bound for the complex Grothendieck constant

The accompanying proof establishes

$$
K_G^{\mathbb C}>\frac{1373}{1000}=1.373.
$$

It specifies a finite complex matrix and bounds its scalar and vector values. Read the [manuscript PDF](manuscript.pdf) or the [LaTeX source](manuscript.tex) for the mathematical argument.

## Verify the certificate

The programs require Python 3 and use only standard-library integer and rational arithmetic. They need no numerical optimizer or data generator. Assertions must remain enabled; do not use Python's `-O` option.

From this directory, run:

```sh
cd certificate
python3 -B -E verify_gram.py
python3 -B -E verify_scalar.py
```

The first program verifies all 16 auxiliary weight tuples. The second verifies the full scalar-parameter interval and all radial tails. Both programs check the same exact data in [family.json](certificate/family.json). A successful run ends with a `PASS` statement; a failed assertion exits nonzero.

The programs also write `verify_gram.log` and `verify_scalar.log` beside their sources. To preserve an unchanged checkout, copy the `certificate` directory to a temporary directory and run the same commands there.

## Build the manuscript

With a standard LaTeX installation, run `pdflatex manuscript.tex` twice from this directory. The source is self-contained and uses no external bibliography or figure files.
