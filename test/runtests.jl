using Wquantiles
import Random
using Test


@testset "Wquantiles (Fidelity)                             " begin
    @test length(detect_ambiguities(Wquantiles)) == 0
end

@testset "Wquantiles (definition: inverse weighted CDF)     " begin
    # Equal weights: the "inverted CDF" quantile.
    @test wquantile([1,2,3,4], [1,1,1,1], [0.0, 0.01, 0.25, 0.26, 0.5, 0.51, 0.75, 0.76, 1.0]) == [1, 1, 1, 2, 2, 3, 3, 4, 4]
    @test wquantile([1,2], [1,1], [0.5]) == [1]
    @test wquantile([1,2], [1,1], [0.500001]) == [2]

    # Unequal weights, hand computed: sorted x = 1,3,5,7,7,9,9,11,11,15 with cumulative
    # weights 1,2,4,5,8,9,14,16,17,18 (out of 18).
    x = [1,3,5,7,9,11,15,7,9,11]
    w = [1,1,2,1,1,2,1,3,5,1]
    @test wquantile(x, w, [1/18, 2/18, 4/18, 5/18, 8/18, 9/18, 14/18, 16/18, 17/18, 1.0]) == [1, 3, 5, 7, 7, 9, 9, 11, 11, 15]
    @test wquantile(x, w, [0.3, 0.45, 0.6, 0.9]) == [7, 9, 9, 11]

    # Unsorted input, and the result is in sorted-q order.
    @test wquantile([5, 1, 3], [1, 1, 1], [0.9, 0.1]) == [1, 5]
end

@testset "Wquantiles (wquantile (vector x, vector w) Types) " begin
    @test wquantile([1,3,5,7,9,11,15,7,9,11], [1,1,2,1,1,2,1,3,5,1], [0.25, 0.5, 0.6, 0.75])             == [7, 9, 9, 9]
    @test wquantile([1,3,5,7,9,11,15,7,9,11], [1.,1.,2.,1.,1.,2.,1.,3.,5.,1.], [1//4, 1//2, 3//5, 3//4]) == [7, 9, 9, 9]
    @test wquantile([1,3,5,7,9,11,15,7,9,11], [1,1,2,1,1,2,1,3,5,1], [1//4, 1//2, 3//5, 3//4])           == [7, 9, 9, 9]
    @test wquantile(UInt8[1,3,5,7,9,11,15,7,9,11], [1,1,2,1,1,2,1,3,5,1], [0.25, 0.5, 0.6, 0.75])        == [0x07, 0x09, 0x09, 0x09]
    @test wquantile([2,4,6,8,10,14,32,4,10,21], [1,1,2,1,1,2,1,3,5,1], [0.25, 0.5, 0.6, 0.75])           == [4, 10, 10, 10]
    @test wquantile(UInt8[1,2,3,4,5,6,7,8,9,10], [1, 1, 2, 1, 1, 2, 1, 3, 5, 1], [0.25, 0.5, 0.6, 0.75]) == [0x04, 0x07, 0x08, 0x09]

    # Rational weights and quantiles stay exact.
    @test wquantile([1,2,3,4], [1//1, 1//1, 1//1, 1//1], [1//2, 3//4]) == [2, 3]

    # Integer quantiles with integer weights.
    @test wquantile([1,2,3], [1,1,1], [0, 1]) == [1, 3]

    # Non-numeric but sortable element types.
    @test wquantile(["b", "a", "c"], [1, 1, 1], [0.5]) == ["b"]
    @test wquantile(Any[3, 1, 2], [1, 1, 1], [0.5]) == [2]
    @test wquantile([[2], [1]], [1, 1], [0.5]) == [[1]]                     # vectors sort lexicographically
    @test_throws Wquantiles.NotSortable wquantile([Set([1]), Set([2])], [1, 1], [0.5])
    @test_throws MethodError wquantile([Set([1]), Set([2])], [1, 1], [0.5]; chk=false)

    # Empty quantile vector.
    @test wquantile([1, 2], [1, 1], Float64[]) == Int[]

    # The return type follows the element type of `x`.
    @test eltype(wquantile([1.5, 2.5], [1, 1], [0.5])) == Float64
end

const TOL = 1.0e-6
Random.seed!(1)
X = rand(10,3)
W = rand(10,3)
qs = [0.1, 0.25, 0.5, 0.75, 0.9]

@testset "Wquantiles (wquantile (matrix X, matrix W))       " begin
    # Consistent with the vector version, column by column.
    R = wquantile(X, W, qs)
    @test size(R) == (5, 3)
    for k in 1:3
        @test R[:, k] == wquantile(X[:, k], W[:, k], qs)
    end

    # Hand computed columns.
    Xm = [1 10; 2 20; 3 30; 4 40]
    Wm = [1 1; 1 1; 1 6; 1 1]
    @test wquantile(Xm, Wm, [0.5, 0.75, 0.9]) == [2 30; 3 30; 4 40]
    @test wquantile(Xm, Wm, [0, 1]) == [1 10; 4 40]
    @test size(wquantile(Xm, Wm, Float64[])) == (0, 2)

    # Error paths.
    @test_throws DimensionMismatch wquantile(Xm, Wm[1:3, :], [0.5])
    @test_throws DomainError wquantile(Xm, [1 1; -1 1; 1 1; 1 1], [0.5])
    @test_throws DomainError wquantile(Xm, [0 1; 0 1; 0 1; 0 1], [0.5])
    @test_throws DomainError wquantile(Xm, Wm, [1.5])
end

@testset "Wquantiles (error paths)                          " begin
    @test_throws DimensionMismatch wquantile([1,2,3], [1,2], [0.5])   # mismatched lengths
    @test_throws DomainError wquantile([1,2,3], [-1,1,1], [0.5])      # negative weight
    @test_throws DomainError wquantile([1,2,3], [0,0,0], [0.5])       # zero-sum weights
    @test_throws DomainError wquantile([1,2,3], [1,1,1], [1.5])       # quantile > 1
    @test_throws DomainError wquantile([1,2,3], [1,1,1], [-0.1])      # quantile < 0
    @test_throws DimensionMismatch Wquantile([1 2; 3 4], [1, 1, 1], [0.5])
    @test_throws DomainError Wquantile([1 2; 3 4], [1, -1], [0.5])
    # With `chk=false` a length mismatch is still a (bounds) error, never a silent read.
    @test_throws BoundsError wquantile([5,4,3,2,1,9,8,7,6,0], [1,1], [0.5]; chk=false)
end

@testset "Wquantiles (edge cases)                           " begin
    @test wquantile([42], [1], [0.0]) == [42]
    @test wquantile([42], [1], [1.0]) == [42]
    @test wquantile([42], [1], [0.5]) == [42]
    @test wquantile([1,2,3], [1,1,1], [0.0]) == [1]
    @test wquantile([1,2,3], [1,1,1], [1.0]) == [3]
    # Zero-weight elements are skipped over.
    @test wquantile([1,2,3], [1,0,1], [0.5, 0.51]) == [1, 3]
end

@testset "Wquantiles (keyword arguments)                    " begin
    x = [1,3,5,7,9,11,15,7,9,11]
    w = [1,1,2,1,1,2,1,3,5,1]
    q = [0.25, 0.5, 0.6, 0.75]
    expected = wquantile(x, w, q)

    # sort_q=false with pre-sorted q
    @test wquantile(x, w, q, sort_q=false) == expected

    # norm_wgt=false with pre-normalized weights
    w_norm = Float64.(w) ./ sum(w)
    @test wquantile(x, w_norm, q, norm_wgt=false) == expected

    # chk=false
    @test wquantile(x, w, q, chk=false) == expected
end

@testset "Wquantiles (Wquantile matrix X, vector w)         " begin
    X_mat = [1 10; 3 30; 5 50; 7 70; 9 90]
    w_vec = [1, 1, 2, 1, 1]
    q_vec = [0.25, 0.5, 0.75]
    result = Wquantile(X_mat, w_vec, q_vec)
    @test size(result) == (3, 2)
    for k in 1:size(X_mat, 2)
        @test result[:, k] == wquantile(X_mat[:, k], w_vec, q_vec)
    end
    @test result == [3 30; 5 50; 7 70]

    # Integer quantiles, and a matrix with no columns.
    @test Wquantile(X_mat, w_vec, [0, 1]) == [1 10; 9 90]
    @test size(Wquantile(zeros(3, 0), [1, 1, 1], [0.5])) == (1, 0)

    # The caller's weight and quantile vectors are not modified.
    w_f = [1.0, 1.0, 2.0, 1.0, 1.0]
    q_u = [0.75, 0.25]
    Wquantile(X_mat, w_f, q_u)
    @test w_f == [1.0, 1.0, 2.0, 1.0, 1.0]
    @test q_u == [0.75, 0.25]
    wquantile(X_mat[:, 1], w_f, q_u)
    @test w_f == [1.0, 1.0, 2.0, 1.0, 1.0] && q_u == [0.75, 0.25]
end

@testset "Wquantiles (wquantile matrix X, matrix W vs vec)  " begin
    Random.seed!(42)
    X2 = rand(20, 4)
    W2 = rand(20, 4)
    q2 = [0.1, 0.25, 0.5, 0.75, 0.9]
    result = wquantile(X2, W2, q2)
    for k in 1:size(X2, 2)
        @test result[:, k] == wquantile(X2[:, k], W2[:, k], q2)
    end
    # Consistent with the threaded version when all weight columns are the same.
    W3 = repeat(W2[:, 1], 1, 4)
    @test wquantile(X2, W3, q2) == Wquantile(X2, W2[:, 1], q2)
end
