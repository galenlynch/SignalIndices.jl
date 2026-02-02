@testset "index and time conversion" begin
    @testset "ndx_to_t" begin
        @test ndx_to_t(1, 1, 0) == 0
        @test ndx_to_t(30000, 30000, 0) == 29999 / 30000
        @test ndx_to_t(1, 1, 20) == 20
        @test ndx_to_t(1:2, 1, 0) == 0.0:1.0:1.0
        @test ndx_to_t(1:1, 1, 20) == 20.0:1.0:20.0
        @test ndx_to_t([1, 2], 1, 0) == [0, 1]

        # Float32 array input
        @test eltype(ndx_to_t([1.0f0, 2.0f0], 1.0f0, 0.0f0)) == Float32

        # Type promotion: Int index + Float64 fs
        @test ndx_to_t(1, 1.0, 0.0) ≈ 0.0
    end

    @testset "t_to_ndx" begin
        @test t_to_ndx(1, 30000, 0) == 30001
        @test t_to_ndx(1, 30000, 1) == 1
        @test t_to_ndx(1:2, 1, 0) == [2, 3]
        @test t_to_ndx([1, 2], 1, 0) == [2, 3]

        # Integer time input
        @test t_to_ndx(0, 1, 0) == 1

        # Exact sample-aligned time
        @test t_to_ndx(1.0, 1.0, 0.0) == 2
    end

    @testset "t_to_last_ndx" begin
        # Scalar
        @test t_to_last_ndx(1.0, 1.0, 0.0) == 2

        # Array
        @test t_to_last_ndx([1.0, 2.0], 1.0, 0.0) == [2, 3]

        # For exact sample times, t_to_last_ndx == t_to_ndx
        @test t_to_last_ndx(1.0, 10.0, 0.0) == t_to_ndx(1.0, 10.0, 0.0)

        # For non-sample-aligned times, t_to_last_ndx <= t_to_ndx
        @test t_to_last_ndx(0.15, 10.0, 0.0) <= t_to_ndx(0.15, 10.0, 0.0)
    end

    @testset "t_sup_to_ndx" begin
        # Scalar: equivalent to t_to_ndx - 1
        @test t_sup_to_ndx(1.0, 1.0, 0.0) == t_to_ndx(1.0, 1.0, 0.0) - 1

        # Array
        @test t_sup_to_ndx([1.0, 2.0], 1.0, 0.0) == t_to_ndx([1.0, 2.0], 1.0, 0.0) .- 1
    end

    @testset "clip_ndx" begin
        @test clip_ndx(-1, 2) == 1
        @test clip_ndx(3, 2) == 2
        @test clip_ndx(1, 2) == 1
        @test clip_ndx(2, 2) == 2
        @test clip_ndx(Int32(1), Int64(2)) == 1
    end

    @testset "clip_ndx_deviance" begin
        # In range: deviance == 0
        @test clip_ndx_deviance(3, 5) == (3, 0)

        # Below range: clipped up, deviance > 0
        ndx, dev = clip_ndx_deviance(-2, 5)
        @test ndx == 1
        @test dev > 0

        # Above range: clipped down, deviance < 0
        ndx, dev = clip_ndx_deviance(10, 5)
        @test ndx == 5
        @test dev < 0
    end

    @testset "n_ndx" begin
        @test n_ndx(1, 1) == 1
        @test n_ndx(1, 2) == 2
    end

    @testset "duration" begin
        @test duration(2, 1) == 1.0
        @test duration(30001, 30000) == 1.0
    end

    @testset "time_interval" begin
        # Integer npoints
        @test time_interval(101, 100, 0.0) == (0.0, 1.0)

        # Float npoints with start_t offset
        @test time_interval(101, 100.0, 5.0) == (5.0, 6.0)

        # Vector input form
        v = rand(101)
        @test time_interval(v, 100, 0.0) == time_interval(101, 100, 0.0)
    end

    @testset "ndx_offset" begin
        @test ndx_offset(1, 1) == 1
        @test ndx_offset(1, 2) == 2
        @test ndx_offset(1, 0) == 0
        @test ndx_offset(3, -3) == 1

        # Negative npt: start_ndx is the end of range
        @test ndx_offset(10, -5) == 6
        @test ndx_offset(5, -1) == 5
    end

    @testset "ndx_wrap" begin
        @test ndx_wrap(1, 5) == 1
        @test ndx_wrap(6, 5) == 1
        @test ndx_wrap(5, 5) == 5

        # Wrap-around at boundaries
        @test ndx_wrap(0, 5) == 5
        @test ndx_wrap(10, 5) == 5
        @test ndx_wrap(11, 5) == 1

        # Values well beyond max
        @test ndx_wrap(106, 5) == 1
    end
end
