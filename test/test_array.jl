@testset "array utilities" begin
    @testset "rev_view" begin
        B = [1, 2, 3, 4, 5]
        rv = rev_view(B)
        @test rv == [5, 4, 3, 2, 1]

        # Verify it's a view (mutating original reflects)
        B[1] = 99
        @test rv[end] == 99
    end

    @testset "weighted_mean" begin
        # Scalar result, equal weights = mean
        @test weighted_mean([2.0, 4.0], [1.0, 1.0]) ≈ 3.0

        # Unequal weights
        @test weighted_mean([0.0, 10.0], [1.0, 3.0]) ≈ 7.5

        # Mismatched sizes throws
        @test_throws ArgumentError weighted_mean([1.0, 2.0], [1.0])
    end

    @testset "weighted_mean_dim" begin
        A = ones(3, 3)
        B = [1.0, 2.0, 3.0]

        # Assert the result value
        result = weighted_mean_dim(A, B)
        @test all(result .≈ 1.0)

        # Test along default dim (last dim for 2D = dim 2)
        A2 = [1.0 2.0; 3.0 4.0]
        w = [1.0, 1.0]
        r2 = weighted_mean_dim(A2, w)
        @test size(r2) == (2, 1)
        @test r2[1, 1] ≈ 1.5
        @test r2[2, 1] ≈ 3.5

        # Test along dim=1
        r1 = weighted_mean_dim(A2, w, 1)
        @test size(r1) == (1, 2)
        @test r1[1, 1] ≈ 2.0
        @test r1[1, 2] ≈ 3.0
    end

    @testset "simple_summary_stats" begin
        # Known input
        m, s, sem = simple_summary_stats([2.0, 4.0, 6.0])
        @test m ≈ 4.0
        @test s ≈ 2.0
        @test sem ≈ 2.0 / sqrt(3)

        # Single element
        m1, s1, sem1 = simple_summary_stats([5.0])
        @test m1 ≈ 5.0
    end

    @testset "find_closest" begin
        @test find_closest([1.0, 2.0, 3.0], 2.1) == 2
        @test find_closest([1.0, 2.0, 3.0], 0.5) == 1

        # Eligible-mask form
        eligible = [false, true, true]
        @test find_closest([1.0, 2.0, 3.0], 1.5, eligible) == 2

        # Function form
        @test find_closest(abs, [-3.0, -1.0, 2.0], -1.5) == 2
    end

    @testset "find_subseq" begin
        # Single-element subseq
        @test find_subseq([3], [1, 2, 3, 4, 3]) == [3, 5]

        # Multi-element
        @test find_subseq([1, 2], [1, 2, 3, 1, 2]) == [1, 4]

        # No match
        @test find_subseq([9], [1, 2, 3]) == []

        # Overlapping occurrences
        @test find_subseq([1, 1], [1, 1, 1]) == [1, 2]

        # Empty subseq throws
        @test_throws ArgumentError find_subseq(Int[], [1, 2, 3])
    end

    @testset "subselect" begin
        v = collect(1:10)

        # Basic extraction
        result = subselect(v, [(1, 3), (5, 7)])
        @test length(result) == 2
        @test result[1] == [1, 2, 3]
        @test result[2] == [5, 6, 7]

        # Multiple ranges
        result2 = subselect(v, [(1, 2), (3, 4), (9, 10)])
        @test length(result2) == 3
    end

    @testset "trailing_zeros_idx" begin
        # Trailing zeros
        @test trailing_zeros_idx([1, 2, 3, 0, 0]) == 3

        # No zeros
        @test trailing_zeros_idx([1, 2, 3]) == 3

        # All zeros
        @test trailing_zeros_idx([0, 0, 0]) == 0
    end

    @testset "clipsize!" begin
        v = collect(1:10)
        clipsize!(v, 5)
        @test length(v) == 5
        @test v == [1, 2, 3, 4, 5]
    end

    @testset "invert_perm" begin
        # Identity permutation
        @test invert_perm([1, 2, 3]) == [1, 2, 3]

        # Known permutation roundtrip
        p = [3, 1, 2]
        ip = invert_perm(p)
        @test ip[p] == 1:3
    end

    @testset "copy_length_dest_check" begin
        # Fits
        @test copy_length_dest_check(10, 1, 5)

        # Too small
        @test !copy_length_dest_check(3, 1, 5)

        # Exact boundary
        @test copy_length_dest_check(5, 1, 5)
    end
end
