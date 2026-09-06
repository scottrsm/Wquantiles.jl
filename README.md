# Wquantiles.jl
A few versions of weighted quantiles with examples and comparisons 
using data parallelization. See `src/WquantileTest.ipynb` for details.

## Quick Start
```julia
using Wquantiles

x = [1.0, 2.0, 3.0, 4.0]
w = [1.0, 1.0, 1.0, 1.0]           # Non-negative weights (normalized internally).
q = [0.25, 0.5, 0.75]

wquantile(x, w, q)                  # Weighted quantiles of a vector.

X = rand(100, 5)                    # One series per column.
Wquantile(X, w100, q)               # Column-wise weighted quantiles (threaded; start Julia with -t).
wquantile(X, W, q)                  # Per-column weights given as a matrix `W`.
```

Exported: `wquantile`, `Wquantile`. Keyword arguments `chk`, `norm_wgt`, and `sort_q`
control input checking, weight normalization, and quantile sorting.
See the Documenter site (or `docs/src/index.md`) for details.
