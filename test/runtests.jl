using Wquantiles
import Random
using Test


@testset "Wquantiles (Fidelity)                             " begin
    @test length(detect_ambiguities(Wquantiles)) == 0
end

@testset "Wquantiles (wquantile (vector x, vector w) Types) " begin

    @test wquantile([1,3,5,7,9,11,15,7,9,11], [1,1,2,1,1,2,1,3,5,1], [0.25, 0.5, 0.6, 0.75])             == [7, 9, 11, 11]
    @test wquantile([1,3,5,7,9,11,15,7,9,11], [1.,1.,2.,1.,1.,2.,1.,3.,5.,1.], [1//4, 1//2, 3//5, 3//4]) == [7, 9, 11, 11]
    @test wquantile([1,3,5,7,9,11,15,7,9,11], [1,1,2,1,1,2,1,3,5,1], [1//4, 1//2, 3//5, 3//4])           == [7, 9, 11, 11]
    @test wquantile(UInt8[1,3,5,7,9,11,15,7,9,11], [1,1,2,1,1,2,1,3,5,1], [0.25, 0.5, 0.6, 0.75])        == [0x07, 0x09, 0x0b, 0x0b]
    @test wquantile([2,4,6,8,10,14,32,4,10,21], [1,1,2,1,1,2,1,3,5,1], [0.25, 0.5, 0.6, 0.75])           == [6   , 10  , 14  , 14  ]
    @test wquantile(UInt8[1,2,3,4,5,6,7,8,9,10], [1, 1, 2, 1, 1, 2, 1, 3, 5, 1], [0.25, 0.5, 0.6, 0.75]) == [0x05, 0x08, 0x09, 0x0a]
end

const TOL = 1.0e-6
Random.seed!(1)
X = rand(10,3)
W = rand(10,3)
qs = [0.1, 0.25, 0.5, 0.75, 0.9]
Result = [0.0491718  0.138227   0.196832;
          0.119079   0.327726   0.661425;
          0.767518   0.347737   0.778149;
          0.855718   0.592041   0.941894;
          0.89077    0.977264   0.983662]

@testset "Wquantiles (wquantile (matrix X, matrix W))       " begin
    @test all(isapprox.(vec(wquantile(X, W, qs)), vec(Result), atol=TOL))
end

@testset "Wquantiles (error paths)                          " begin
    @test_throws DomainError wquantile([1,2,3], [1,2], [0.5])         # mismatched lengths
    @test_throws DomainError wquantile([1,2,3], [-1,1,1], [0.5])      # negative weight
    @test_throws DomainError wquantile([1,2,3], [0,0,0], [0.5])       # zero-sum weights
    @test_throws DomainError wquantile([1,2,3], [1,1,1], [1.5])       # quantile > 1
    @test_throws DomainError wquantile([1,2,3], [1,1,1], [-0.1])      # quantile < 0
end

@testset "Wquantiles (edge cases)                           " begin
    @test wquantile([42], [1], [0.0]) == [42]
    @test wquantile([42], [1], [1.0]) == [42]
    @test wquantile([42], [1], [0.5]) == [42]
    @test wquantile([1,2,3], [1,1,1], [0.0]) == [1]
    @test wquantile([1,2,3], [1,1,1], [1.0]) == [3]
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
    for k in 1:size(X_mat, 2)
        @test result[:, k] == wquantile(X_mat[:, k], w_vec, q_vec)
    end
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
end

