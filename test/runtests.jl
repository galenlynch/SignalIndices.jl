using SignalIndices
using Test

@testset "SignalIndices.jl" begin
    @testset "types" begin
        @test div_type(Int) == Float64
        @test div_type(Int, Int) == Float64
        @test div_type(Float32) == Float32
        @test div_type(Float32, Int) == Float32
        @test div_type(Float32, Float64) == Float64
        @test div_type(2, 5) == Float64
        @test div_type(1.0f0, 1) == Float32
    end

    @testset "indices" begin
        @test clip_ndx(-1, 2) == 1
        @test clip_ndx(3, 2) == 2
        @test clip_ndx(1, 2) == 1
        @test clip_ndx(2, 2) == 2
        @test clip_ndx(Int32(1), Int64(2)) == 1

        @test ndx_to_t(1, 1, 0) == 0
        @test ndx_to_t(30000, 30000, 0) == 29999 / 30000
        @test ndx_to_t(1, 1, 20) == 20
        @test ndx_to_t(1:2, 1, 0) == 0.0:1.0:1.0
        @test ndx_to_t(1:1, 1, 20) == 20.0:1.0:20.0
        @test ndx_to_t([1, 2], 1, 0) == [0, 1]

        @test t_to_ndx(1, 30000, 0) == 30001
        @test t_to_ndx(1, 30000, 1) == 1
        @test t_to_ndx(1:2, 1, 0) == [2, 3]
        @test t_to_ndx([1, 2], 1, 0) == [2, 3]

        @test n_ndx(1, 1) == 1
        @test n_ndx(1, 2) == 2

        @test duration(2, 1) == 1.0
        @test duration(30001, 30000) == 1.0

        @test ndx_offset(1, 1) == 1
        @test ndx_offset(1, 2) == 2
        @test ndx_offset(1, 0) == 0
        @test ndx_offset(3, -3) == 1

        @test bin_bounds(1, 1024) == (1, 1024)
        @test bin_bounds(2, 1024) == (1025, 2048)
        @test bin_bounds(1, 1024, 1023) == (1, 1023)
        @test bin_bounds(1:2, 1024) == (1:1024:1025, 1024:1024:2048)

        @test bin_center(1, 1024) == 512.5
        @test bin_center([(1, 1024)]) == [512.5]
        @test bin_center(1:2, 1024) == 512.5:1024.0:1536.5

        @test make_slice_idx(2, 1, 1) == (1, :)
        @test make_slice_idx(2, 1, 2) == (2, :)
        @test make_slice_idx(3, 1, 1:2) == (1:2, :, :)

        @test make_expand_idx(2, 1) == (:, 1)
        @test make_expand_idx(2, 2) == (1, :)

        @test copy_length_check(5, 1)
        @test !copy_length_check(1, 5)
        @test copy_length_check(rand(5), rand(1))
        @test !copy_length_check(rand(1), rand(5))

        @test ndx_wrap(1, 5) == 1
        @test ndx_wrap(6, 5) == 1
        @test ndx_wrap(5, 5) == 5

        # Internal function test
        SignalIndices.view_trailing_slice_impl(Array{Int,2})
    end

    @testset "array" begin
        A = rand(3, 3)
        B = rand(3)

        @test all(B[end:-1:1] .== rev_view(B))

        weighted_mean_dim(A, B)

        C = [2, 1, 2, 3, 2]
        @test local_extrema(C) == [4]
        @test local_extrema(C, <) == [2]
        @test local_extrema(C[1:4]) == []
    end

    @testset "moving_sum" begin
        @test moving_sum([1, 2, 3, 4, 5], 3) == [6, 9, 12]
        @test moving_sum([1, 2, 3], 1) == [1, 2, 3]
        @test moving_sum([1, 2, 3], 0) == [1, 2, 3]
    end

    @testset "pairwise" begin
        @test pairwise_idxs(3) == [(2, 1), (3, 1), (3, 2)]
        @test pairwise_idx(2, 1, 3) == 1
        @test pairwise_idx(3, 1, 3) == 2
        @test pairwise_idx(3, 2, 3) == 3
        @test map_pairwise(-, [1, 2, 3]) == [1, 2, 1]
    end

    @testset "skipnothing" begin
        @test sum(skipnothing([1, nothing, 2, nothing, 3])) == 6
        @test collect(skipnothing([1, nothing, 2])) == [1, 2]
    end

    @testset "find_closest" begin
        @test find_closest([1.0, 2.0, 3.0], 2.1) == 2
        @test find_closest([1.0, 2.0, 3.0], 0.5) == 1
    end

    @testset "uniformhist" begin
        @test uniformhist([0.5, 1.5, 2.5], 0:1:4) == [1, 1, 1, 0]
        @test uniformhist([0.5, 0.6, 2.5], 0:1:3) == [2, 0, 1]
    end

    @testset "edge_triggers" begin
        sig = [0, 0, 1, 1, 0, 0, 1, 1]
        @test find_all_edge_triggers(sig, 1) == [3, 7]
        @test find_first_edge_trigger(sig, 1) == 3
        @test find_first_edge_trigger([0, 0, 0], 1) === nothing
    end

    @testset "thresh_cross" begin
        sig = [0, 1, 0, 1, 0]
        @test thresh_cross(sig, 1, <) == [2, 4]
    end

    @testset "filtermap" begin
        @test filtermap(x -> x > 0, x -> x^2, [-1, 2, -3, 4]) == [4, 16]
    end

    @testset "find_not_unique" begin
        @test sort(find_not_unique([1, 2, 1, 3, 2])) == [1, 2, 3, 5]
        @test find_not_unique([1, 2, 3]) == []
    end
end
