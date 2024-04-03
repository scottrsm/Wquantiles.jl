# Wquantiles.jl Documentation

```@meta
CurrentModule = Wquantiles
```
# Overview
This module contains functions to compute weighted quantiles of data series.
One version of this, `Wquantile` takes advantage of parallelization.
To use it make sure that you are running Julia with multiple threads.

There is an associated Jupyter notebook at 
/src/WquantileTest.ipynb.

## Functions

```@docs
wquantile(::AbstractVector{T}, ::AbstractVector{S}, ::AbstractVector{V}; ::Bool = true, ::Bool = true, ::Bool = true) where {T, S <: Real, V <: Real}
```


```@docs
Wquantile(::AbstractMatrix{T}, ::AbstractVector{S}, ::AbstractVector{V}; ::Bool = true, ::Bool = true, ::Bool = true) where {T, S <: Real, V <: Real}
```

```@docs
wquantile(::AbstractMatrix{T}, ::AbstractMatrix{S}, ::AbstractVector{V}; ::Bool = true, ::Bool = true, ::Bool = true) where {T, S <: Real, V <: Real}
```


## Index

```@index
```

