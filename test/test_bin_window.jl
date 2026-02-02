@testset "bin, window, and slicing" begin
    @testset "bin_bounds" begin
        @test bin_bounds(1, 1024) == (1, 1024)
        @test bin_bounds(2, 1024) == (1025, 2048)
        @test bin_bounds(1, 1024, 1023) == (1, 1023)
        @test bin_bounds(1:2, 1024) == (1:1024:1025, 1024:1024:2048)
    end

    @testset "bin_center" begin
        @test bin_center(1, 1024) == 512.5
        @test bin_center([(1, 1024)]) == [512.5]
        @test bin_center(1:2, 1024) == 512.5:1024.0:1536.5
    end

    @testset "expand_selection" begin
        # Expand within bounds
        @test expand_selection(5, 10, 20, 3) == (2, 13)

        # Clamped at left edge
        @test expand_selection(2, 10, 20, 5) == (1, 15)

        # Clamped at right edge
        @test expand_selection(5, 18, 20, 5) == (1, 20)

        # Clamped both sides
        @test expand_selection(1, 10, 10, 5) == (1, 10)

        # Type promotion
        @test expand_selection(Int32(5), Int32(10), Int64(20), Int64(3)) == (2, 13)
    end

    @testset "view_trailing_slice" begin
        # 2D array
        A = reshape(1:12, 3, 4)
        @test view_trailing_slice(A, 2) == A[:, 2]

        # 3D array
        B = reshape(1:24, 2, 3, 4)
        @test view_trailing_slice(B, 3) == B[:, :, 3]

        # Range index returns sub-array
        @test view_trailing_slice(A, 1:2) == A[:, 1:2]
    end

    @testset "make_slice_idx" begin
        @test make_slice_idx(2, 1, 1) == (1, :)
        @test make_slice_idx(2, 1, 2) == (2, :)
        @test make_slice_idx(3, 1, 1:2) == (1:2, :, :)

        # 1D edge case
        @test make_slice_idx(1, 1, 3) == (3,)
    end

    @testset "make_expand_idx" begin
        @test make_expand_idx(2, 1) == (:, 1)
        @test make_expand_idx(2, 2) == (1, :)
    end

    @testset "copy_length_check" begin
        @test copy_length_check(5, 1)
        @test !copy_length_check(1, 5)
        @test copy_length_check(rand(5), rand(1))
        @test !copy_length_check(rand(1), rand(5))

        # Offset variants
        @test copy_length_check(5, 3, 1, 1, 3)
        @test !copy_length_check(2, 5, 1, 1, 5)

        # Array-with-offset form
        @test copy_length_check(rand(5), 1, rand(3), 1, 3)
    end

    @testset "copy_length_dest_check" begin
        # Fits
        @test copy_length_dest_check(10, 1, 5)

        # Doesn't fit
        @test !copy_length_dest_check(3, 1, 5)

        # Edge-exact fit
        @test copy_length_dest_check(5, 1, 5)
    end
end
