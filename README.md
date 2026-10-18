# LaughlinGapCertificate

Verifier and certificate for the thermodynamic limit SDP for the Laughlin gap, accompanying
"A Hierarchy of Exact Lower Bounds on the Laughlin Spectral Gap" by Michael G. Scheer and Eslam Khalaf
([arXiv:2610.18854](https://arxiv.org/abs/2610.18854)).

## Usage

From the root directory of the repository, install all dependencies once with

```
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then, within the `certificate` directory, extract the certificate data with

```
unzip certificates_15_22.zip
```

yielding the file `certificates_15_22.json`. One can verify the `sha256sum` with

```
sha256sum certificates_15_22.json
```

which should print

```
2495c0add7ea0fc40b440e0ec1bcb0d5869e3fab38c0d4d38fcd9a67722a7b99  certificates_15_22.json
```

After that, run

```
julia --project=.. certify.jl
```

to certify the lower bound. The code has been tested on Julia 1.13.1 with dependencies pinned by `Manifest.toml`.
The commands `unzip` and `sha256sum` are standard on Linux. Equivalent commands exist on other operating systems.

## Description

The `certify.jl` program constructs the SDP as a function of `t` using the code in the `src` directory. It sets `t = 0` and
rescales the resulting rational SDP using a saved preconditioner. It then certifies a saved exact rational
interior point by checking exact equality and positivity. Finally, it computes the corresponding rational
objective value, evaluates the tail bound with interval arithmetic, and certifies `ϵ > 1/100`. The algorithms
used to produce the exact certificate are described in the paper.

The saved data is stored in the JSON text file `certificates_15_22.json`. This file additionally contains fields
`rank_stability.indep_row_inds` and `rank_stability.indep_col_inds` which are row and column indices (1-based) for a
maximal rank invertible square submatrix of the scaled `t = 0` equality condition matrix `csdp.R` in `certify.jl`.
The invertibility of this submatrix can be verified by computing the rank modulo the prime in `rank_stability.invertible_prime`,
which divides none of the denominators of `csdp.R`. This submatrix was used in the SDP certification algorithm and also
in the rank stability algorithm as described in the paper, though rank stability certificates are not provided here.

The code's successful output (printed to `stdout`) ends with:
```
t = 0 SDP and tail bound exactly certify ϵ > 1/100: true
```

## Performance

On a laptop with 8 CPUs and 15.3 GB of memory (as reported by `htop`), using 3 worker processes,
the verification completes in around seven minutes (this does not include installation and precompilation time for dependencies).
The number of worker processes is set in `certify.jl` by the variable `num_workers` which defaults to 3.
One can lower this to reduce memory usage at the expense of greater runtime.

## Citation

If you use this code or certificate in a publication, please cite:

```bibtex
@misc{2610.18854,
  author        = {Michael G. Scheer and Eslam Khalaf},
  title         = {A hierarchy of exact lower bounds on the Laughlin spectral gap},
  year          = {2026},
  eprint        = {2610.18854},
  archivePrefix = {arXiv},
  primaryClass  = {cond-mat.str-el},
  url           = {https://arxiv.org/abs/2610.18854},
}
```
