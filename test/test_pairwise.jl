@testset "pairwise operations" begin
    @testset "pairwise_idxs" begin
        @test pairwise_idxs(3) == [(2, 1), (3, 1), (3, 2)]

        # n=1: empty
        @test pairwise_idxs(1) == []

        # n=2: single pair
        @test pairwise_idxs(2) == [(2, 1)]

        # Verify count = n*(n-1)/2
        @test length(pairwise_idxs(5)) == 10
    end

    @testset "pairwise_idx" begin
        @test pairwise_idx(2, 1, 3) == 1
        @test pairwise_idx(3, 1, 3) == 2
        @test pairwise_idx(3, 2, 3) == 3

        # i==j throws
        @test_throws ArgumentError pairwise_idx(2, 2, 3)

        # Reversed order (j > i) gives same result
        @test pairwise_idx(1, 3, 3) == pairwise_idx(3, 1, 3)
    end

    @testset "map_pairwise" begin
        # Single-vector form
        @test map_pairwise(-, [1, 2, 3]) == [1, 2, 1]

        # Two-vector form: returns matrix with correct dimensions
        as = [1, 2, 3]
        bs = [10, 20]
        result = map_pairwise(-, as, bs)
        @test size(result) == (2, 3)
        @test result[1, 1] == 1 - 10
        @test result[2, 1] == 1 - 20
        @test result[1, 3] == 3 - 10
    end

    @testset "imap_product" begin
        # Collect returns a matrix from Iterators.product
        result = collect(imap_product(+, [1, 2], [10, 20]))
        @test result == [11 21; 12 22]
    end
end
