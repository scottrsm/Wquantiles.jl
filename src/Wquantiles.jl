module Wquantiles

export wquantile, Wquantile

import Folds

# Exception for non-sortable types.
struct NotSortable <: Exception
    var::String
end

# How to print NotSortable Exception.
Base.show(io::IO, e::NotSortable) = print(io, "Type, \"$(e.var)\", is NOT sortable.")

# Make a sortable Trait: does the type have an `isless` method?
# (For a non-concrete element type we cannot tell up front; let `sortperm` decide.)
function isSortable(::Type{T}) where {T}
    return !isconcretetype(T) || hasmethod(isless, Tuple{T, T})
end

# The element type used for (normalized) weights: the promotion of the weight and
# quantile types, closed under division (so integer weights become floating point,
# while rational weights stay rational).
function _weight_type(::Type{S}, ::Type{V}) where {S, V}
    P = promote_type(S, V)
    return typeof(one(P) / one(P))
end

#= Definition of the weighted quantile used throughout this module.

   Let `xs` be `x` sorted (ascending) with the (normalized) weights `ws` permuted alongside,
   and let `F(i) = ws[1] + ... + ws[i]` be the cumulative weight up to and including `xs[i]`.
   The `q` quantile is `xs[i*]` with `i* = min{ i : F(i) >= q }`; if rounding keeps every
   `F(i)` below `q` (only possible when `q` is 1), the largest element is returned.
   This is the inverse of the weighted empirical distribution function; with equal weights
   it agrees with the "inverted CDF" (type 1) quantile of `Statistics`/`StatsBase`.
=#


"""
    wquantile(x, w, q[; chk=true, norm_wgt=true, sort_q=true])

Finds the `q` weighted quantile values from the vector `x`.

# Type Constraints
- `S <: Real`
- `V <: Real`

# Arguments
- `x  ::AbstractVector{T}`: Vector(n) of values from which to find quantiles.
- `w  ::AbstractVector{S}`: Vector(n) of weights to use.
- `q  ::AbstractVector{V}`: Vector(l) of quantile values.

# Keyword Arguments
- `chk=true     :: Bool`: If `true`, check the input contract described below.
- `norm_wgt=true:: Bool`: If `true`, normalize the weights.

   **NOTE:** If `norm_wgt` is `false`, it is *ASSUMED* that `w` is already normalized.
- `sort_q=true  ::Bool`: If `true`, sort the quantile vector, `q`.

   **NOTE:** If `sort_q` is `false`, it *ASSUMED* that `q` is already sorted.

# Input Contract
- The type of `x` implements sortable.
-    `|x| == |w|`     -- Length of `x` matches length of weights.
- `∀i,  w[i] >= 0`    -- The weights are non-negative.
-    `Σ w[i] >  0`    -- The sum of the weights is positive.
- `∀i,  q[i] <= 1`    -- The quantile values are in ``[0,1]``.
- `∀i,  q[i] >= 0`

# Return
The vector(l) of weighted quantile values from `x`.
Letting `qs` be the sorted quantiles of `q`.
The entry `i` is the ``i^{\\rm th}`` quantile (in `qs`) of `x`.

# Definition
Let `xs` be `x` sorted (ascending) with the normalized weights `ws` permuted alongside,
and `F(i) = ws[1] + ... + ws[i]` the cumulative weight up to and including `xs[i]`.
The quantile for `q` is `xs[i]` for the smallest `i` with `F(i) >= q` (the inverse of the
weighted empirical distribution function). With equal weights this is the "inverted CDF"
quantile: for `x = [1, 2, 3, 4]`, `q = 0.5` gives `2` and `q = 0.51` gives `3`.
"""
function wquantile(x::AbstractVector{T} , 
                   w::AbstractVector{S} , 
                   q::AbstractVector{V} ;
                   chk::Bool = true     , 
                   norm_wgt::Bool = true, 
                   sort_q::Bool = true   ) where {T, S <: Real, V <: Real}
    # We report back the quantiles of `x` in sorted `q` order, so we need to sort `q`.
    # **NOTE:** If we don't explicitly sort `q`, it means that you are *ASSUMING* `q` is sorted.
    qs = sort_q ? sort(q) : q
    m  = length(qs)
    n  = length(x)

    # Get 0 and 1 for `q` types.
    zeroq = zero(V)
    oneq  = one(V)

    # Check input contract...
    if chk
        !isSortable(T)  && throw(NotSortable(string(T)))
        n != length(w)  && throw(DimensionMismatch("`x` and `w` do not have the same length: $n != $(length(w))."))
    end

    # Get the permutation of indices that sort `x`.
    idx = sortperm(x)

    # Sorted weights, converted to a type closed under division (a fresh array: `w` is never modified).
    W   = _weight_type(S, V)
    wsc = W[w[i] for i in idx]
    zerow = zero(W)

    # Check input contract...
    if chk
        !all(>=(zerow), wsc)         && throw(DomainError(0, "`w`: Some weights are negative."))
        !(sum(wsc) >  zerow)         && throw(DomainError(0, "`w`: The sum of the weights is NOT > 0."))
        !all(v -> zeroq <= v <= oneq, qs) && throw(DomainError(0, "`q`  : Some quantiles are NOT in the interval, [0,1].")) 
    end

    # Normalize sorted weights?
    if norm_wgt
      wsc ./= sum(wsc)
    end

    # Apply permutation to `x`.
    xs = x[idx]

    #= Create an index vector to get the list of quantiles of `x`.
       Default the indices to the largest element of `x`.
       Why is this important: If one chooses 1.0 as a quantile, it could easily be
                              the case, due to numeric inaccuracy, that we do not
                              reach the 1.0 threshold. For this reason we want
                              to pick the default index value to
                              be the largest index in `x`, so if the threshold
                              is not reached we do the right thing.
	=#
    qxsi = fill(n, m)

    # Using the fact that the quantile values are in sorted order,
    # find the index for each associated value in `xs`, placing them in `qxsi`.
    # `s` is the cumulative weight up to and including `xs[i]`.
    j = 1
    s = zerow
    m == 0 && @goto done
    @inbounds for i in 1:n
        s += wsc[i]
        # If we reach the current quantile threshold, `qs[j]`,
        # set the index at `j` of the index vector (and of any further quantiles also reached).
        if s >= qs[j]
            qxsi[j] = i
            while true
                j += 1
                if j == m+1
                    @goto done
                end
                if qs[j] > s
                    break
                end
                qxsi[j] = i
            end
        end
    end
    @label done

    # Return the quantile values (in quantile sorted order).
    return(xs[qxsi])
end


"""
    Wquantile(X, w, q; chk=true, norm_wgt=true, sort_q=true)

Finds the `q` weighted quantile values from the columns of the matrix `X`.

# Type Constraints
- `S <: Real`
- `V <: Real`

# Arguments
- `X  ::AbstractMatrix{T}`: Matrix(n,m) of values from which to find quantiles.
- `w  ::AbstractVector{S}`: Vector(n) of weights to use.
- `q  ::AbstractVector{V}`: Vector(l) of quantile values.

# Keyword Arguments
- `chk=true     ::Bool`: If `true`, check the input contract described below.
- `norm_wgt=true::Bool`: If `true`, normalize the weights.

   **NOTE:** If `norm_wgt` is `false`, it is *ASSUMED* that `w` is already normalized.
- `sort_q=true  ::Bool`: If `true`, sort the quantile vector, `q`.

   **NOTE:** If `sort_q` is `false`, it is *ASSUMED* that `q` is already sorted.

# Input Contract
-  The type of `X` implements sortable.
- `∀i, |X[:, i]|  == |w|` -- Length of each column of `X` matches length of weights.
- `∀i,      w[i]  >= 0`   -- Weights are non-negative.
- `Σ w[i]  >  0`          -- The sum of the weights is positive.
- `∀i,      q[i]  <= 1`   -- The quantile values are in ``[0,1]``.
- `∀i,      q[i]  >= 0`

# Return
The `(l,m)` matrix of weighted quantile values from `X`.
Letting `qs` be the sorted quantiles of `q`.
The entry `(i,j)` is the ``i^{\\rm th}`` quantile (in `qs`) from the ``j^{\\rm th}`` column of `X`.
See `wquantile` for the definition of the weighted quantile.

"""
function Wquantile(X::AbstractMatrix{T} , 
                   w::AbstractVector{S} , 
                   q::AbstractVector{V} ;
                   chk::Bool = true     , 
                   norm_wgt::Bool = true, 
                   sort_q::Bool = true   ) where {T, S <: Real, V <: Real}

    n, m = size(X)

    # Check the parts of the input contract that are shared by all columns, once.
    if chk
        !isSortable(T)  && throw(NotSortable(string(T)))
        n != length(w)  && throw(DimensionMismatch("The columns of `X` and `w` do not have the same length: $n != $(length(w))."))
    end

    # Normalize the weights if needed -- into a new array; `w` is never modified.
    # (Fresh names, `wn`/`qs`, are captured below so the closure is type stable.)
    W  = _weight_type(S, V)
    wn = norm_wgt ? W.(w) ./ sum(w) : W.(w)

    # Sort the quantiles if needed.
    qs = sort_q ? sort(q) : q

    # Nothing to do for a matrix with no columns.
    m == 0 && return(Matrix{T}(undef, length(qs), 0))

    #= Create a closure that will be threaded -- computing the weighted quantiles of the columns of `X`.
       If `chk` is true, only do the (remaining) input check for the first column
       as checking the rest of the columns is redundant.
	=#
    wquant_vec_func = p -> wquantile(p[1], wn, qs, 
                                     chk=(p[2] == 1 ? chk : false), 
                                     norm_wgt=false, sort_q=false)

    #= Computation: (from right to left)
       - Zip up the columns of `X` along with the column number.
       - Use Folds.map to apply multiple threads to compute the weighted quantiles on each column of `X`.
       - Place them back as an array using reduce hcat.
    =#
    return(reduce(hcat, Folds.map(wquant_vec_func, collect(zip(eachcol(X), 1:m)))))

end


"""
    wquantile(X, W, q; chk=true)

Finds the `q` weighted quantile values from the columns of the matrix `X`.

# Type Constraints
- `S <: Real`
- `V <: Real`

# Arguments
- `X  ::AbstractMatrix{T}`: Matrix(n,m) of values from which to find quantiles.
- `W  ::AbstractMatrix{S}`: Matrix(n,m) of weights to use.
- `q  ::AbstractVector{V}`: Vector(l) of quantile values.

# Keyword Arguments
- `chk=true::Bool`     : If `true`, check the input contract described below.

# Input Contract
-  The type of `X` is sortable.
-  size(X) == size(W)
- `∀i∀j, W[i, j] >= 0`   -- Weights are non-negative in each column.
- `∀i,Σ_j w[i,j] >  0`   -- The column sum of the weights is positive.
- `∀i, 0 <= q[i] <= 1`   -- The quantile values are in ``[0,1]``.

# Return
The `(l,m)` matrix of weighted quantile values from `X`.
Letting `qs` be the sorted quantiles of `q`.
The entry `(i,j)` is the ``i^{\\rm th}`` quantile (in `qs`) from the ``j^{\\rm th}`` column of `X`.
See `wquantile` for the definition of the weighted quantile.

"""
function wquantile(X::AbstractMatrix{T}, 
                   W::AbstractMatrix{S}, 
                   q::AbstractVector{V}; 
                   chk::Bool = true     ) where {T, S <: Real, V <: Real}
    #= ----------------------------------------------------------------
       --- Destructure inputs and potentially check input contract. ---
       ----------------------------------------------------------------
	=#
    
    # We report back the quantiles of `X` in sorted `q` order, so we sort `q`.
    qs     = sort(q)
    l      = length(qs)
    n, m   = size(X)
    nw, mw = size(W)

    # Get 0 and 1 for `q` types.
    zeroq = zero(V)
    oneq  = one(V)

    # Check input contract...
    if chk
        !isSortable(T)     && throw(NotSortable(string(T)))
        (n, m) != (nw, mw) && throw(DimensionMismatch("`X` and `W` do not have the same size: $((n, m)) != $((nw, mw))."))
    end

    # Convert weights (to a type closed under division).
    WT = _weight_type(S, V)
    WC = WT.(W)
    zerow = zero(WT)

    # Check input contract...
    if chk
        !all(>=(zerow), WC)                            && throw(DomainError(0, "`W`: Some weights are negative."))
        !all(>(zerow), sum(WC, dims=1))                && throw(DomainError(0, "`W`: At least one weight column does NOT have a positive sum."))
        !all(v -> zeroq <= v <= oneq, qs)              && throw(DomainError(0, "`q`: Some quantiles are NOT in the interval, [0,1].")) 
    end

    #= ----------------------------------------------------------------
       --- Sort X, W based on X and compute quantiles.
       ----------------------------------------------------------------
	=#
    
    # Get the permutation of indices that sort the columns of `X`. This will be used to create
    # a matrix of permutations of `wc` that align with this sorting.
    Idx = sortperm(X; dims=1)

    # `Idx` holds linear indices into `X` (and therefore into `WC`, which has the same shape).
    Wsc = WC[Idx]

    # Normalize sorted weights by column.
    Wsc ./= sum(Wsc, dims=1)

    # Apply permutation to `X` -- sorting each column of `X`.
    Xs = X[Idx]

    #= Create an index matrix to get the list of quantiles of `X`.
       Default the indices to the index of the largest element of X for each column.
       Why is this important: If one chooses 1.0 as a quantile, it could easily be
                              the case, due to numeric inaccuracy, that we do not
                              reach the 1.0 threshold. For this reason we want
                              to pick the default index value for each column to
                              be the largest index in the column, so if the threshold
                              is not reached we do the right thing.
	=#
    Qxsi = Array{Int}(undef, l, m)
    for i in 1:m
      @inbounds Qxsi[:, i] .= i * n
    end

    # Using the fact that the quantile values are in sorted order,
    # find the index for each associated value in `Xs`, placing them in `Qxsi`.

    # For each column...
    l == 0 && return(Xs[Qxsi])
    for k in 1:m
        j = 1
        s = zerow   # Cumulative weight up to and including row `i` of column `k`.

        # For each row...
        for i in 1:n
            @inbounds s += Wsc[i, k]
            # If we reach the current quantile threshold, `qs[j]`,
            # set the index at (`j`,`k`) of the index matrix (and of any further quantiles also reached).
            if s >= qs[j]
                @inbounds Qxsi[j,k] = i + (k-1) * n
                while true
                    j += 1
                    if j == l+1
                        @goto column_done
                    end
                    if qs[j] > s
                        break
                    end
                    @inbounds Qxsi[j,k] = i + (k-1) * n
                end
            end
        end

        # We've finished off a column, onto the next.
        @label column_done
    end

    # Return the quantile values as an (`l`,`m`) matrix in quantile sorted order.
    return(Xs[Qxsi])
end


end # module Wquantiles

