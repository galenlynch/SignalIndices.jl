@testset "signal processing" begin
    @testset "moving_sum" begin
        @test moving_sum([1, 2, 3, 4, 5], 3) == [6, 9, 12]
        @test moving_sum([1, 2, 3], 1) == [1, 2, 3]
        @test moving_sum([1, 2, 3], 0) == [1, 2, 3]

        # Empty input
        @test moving_sum(Int[], 3) == Int[]

        # Single element
        @test moving_sum([42], 1) == [42]
        @test moving_sum([42], 5) == [42]

        # Window == length
        @test moving_sum([1, 2, 3], 3) == [6]

        # Window > length
        @test moving_sum([1, 2, 3], 5) == [6]

        # Negative window throws
        @test_throws ArgumentError moving_sum([1, 2, 3], -1)
    end

    @testset "moving_sum!" begin
        s = [1, 2, 3, 4, 5]
        out = zeros(Int, 3)
        moving_sum!(out, s, 3)
        @test out == [6, 9, 12]

        # Verify matches moving_sum
        @test out == moving_sum(s, 3)
    end

    @testset "local_extrema" begin
        C = [2, 1, 2, 3, 2]
        @test local_extrema(C) == [4]
        @test local_extrema(C, <) == [2]
        @test local_extrema(C[1:4]) == []

        # Plateau: the >= comparison sees index 3 as extremum (2 >= 1 after 2 !> 2)
        @test local_extrema([1, 2, 2, 1]) == [3]

        # Monotonic
        @test local_extrema([1, 2, 3, 4]) == []

        # Length <= 2
        @test local_extrema([1, 2]) == []
        @test local_extrema([1]) == []
    end

    @testset "find_local_extrema" begin
        # Climb to max
        sig = [1, 3, 5, 4, 2]
        @test find_local_extrema(sig, 2; findmax=true) == 3

        # Climb to min
        sig = [5, 3, 1, 2, 4]
        @test find_local_extrema(sig, 4; findmax=false) == 3

        # right_on_ties: starting at a tied position, bias determines direction
        sig = [1, 3, 3, 1]
        @test find_local_extrema(sig, 2; findmax=true, right_on_ties=true) == 2
        @test find_local_extrema(sig, 3; findmax=true, right_on_ties=false) == 3

        # Start at extremum (returns immediately)
        sig = [1, 5, 1]
        @test find_local_extrema(sig, 2; findmax=true) == 2

        # Short signal (length 1)
        @test find_local_extrema([42], 1) == 1

        # Missing/NaN skipping
        sig_m = [missing, 1, 3, 2]
        @test find_local_extrema(sig_m, 1; findmax=true) == 3
    end

    @testset "find_all_edge_triggers" begin
        sig = [0, 0, 1, 1, 0, 0, 1, 1]
        @test find_all_edge_triggers(sig, 1) == [3, 7]

        # Custom comparator
        @test find_all_edge_triggers([0, 1, 0, 1], 1, >) == []

        # No crossings
        @test find_all_edge_triggers([0, 0, 0], 1) == []

        # Single-element input
        @test find_all_edge_triggers([1], 1) == []
    end

    @testset "find_first_edge_trigger" begin
        sig = [0, 0, 1, 1, 0, 0, 1, 1]
        @test find_first_edge_trigger(sig, 1) == 3
        @test find_first_edge_trigger([0, 0, 0], 1) === nothing

        # Custom comparator
        @test find_first_edge_trigger([1, 0, 1], 1, >) === nothing
    end

    @testset "thresh_cross" begin
        sig = [0, 1, 0, 1, 0]
        @test thresh_cross(sig, 1, <) == [2, 4]

        # Empty input
        @test thresh_cross(Int[], 1) == Int[]

        # No crossings
        @test thresh_cross([0, 0, 0], 1) == Int[]
    end

    @testset "indices_above_thresh" begin
        # Contiguous spans
        @test indices_above_thresh([0, 1, 2, 0, 3, 4, 0], 1) == [2:3, 5:6]

        # Single-element spans
        @test indices_above_thresh([0, 5, 0, 5, 0], 1) == [2:2, 4:4]

        # Entire array above
        @test indices_above_thresh([5, 5, 5], 1) == [1:3]

        # None above
        @test indices_above_thresh([0, 0, 0], 1) == []

        # Ends above threshold (range extends to end)
        @test indices_above_thresh([0, 1, 2], 1) == [2:3]
    end

    @testset "filter_no_collisions" begin
        # No collisions (all pass)
        @test filter_no_collisions([1, 5, 10], [20, 30], 2) == [1, 5, 10]

        # All collide
        @test filter_no_collisions([1, 2, 3], [1, 2, 3], 0) == Int[]

        # Partial
        @test filter_no_collisions([1, 5, 10], [5], 1) == [1, 10]

        # Empty inputs
        @test filter_no_collisions(Int[], [1, 2], 1) == Int[]
        @test filter_no_collisions([1, 2, 3], Int[], 1) == [1, 2, 3]
    end

    @testset "window_counts" begin
        # Basic counting
        ts = [1.0, 2.0, 3.0, 5.0, 10.0]
        cnts = window_counts(ts, 2.0)
        @test cnts[1] == 3  # [1,2,3] within 1.0+2.0
        @test cnts[5] == 1  # only 10.0

        # Single event
        @test window_counts([5.0], 1.0) == [1]

        # 4-arg form with range bounds
        cnts2, ib = window_counts(ts, 2.0, 2.0, 5.0)
        @test ib == 2
        @test length(cnts2) == 3  # elements 2.0, 3.0, 5.0
    end

    @testset "uniformhist" begin
        @test uniformhist([0.5, 1.5, 2.5], 0:1:4) == [1, 1, 1, 0]
        @test uniformhist([0.5, 0.6, 2.5], 0:1:3) == [2, 0, 1]

        # Out-of-range values ignored
        @test uniformhist([-10.0, 0.5, 100.0], 0:1:3) == [1, 0, 0]

        # Empty input
        @test uniformhist(Float64[], 0:1:3) == [0, 0, 0]

        # Float64 eltype variant
        @test uniformhist(Float64, [0.5, 1.5], 0:1:3) == [1.0, 1.0, 0.0]
        @test eltype(uniformhist(Float64, [0.5], 0:1:2)) == Float64
    end

    @testset "uniformhist!" begin
        # Pre-allocated
        cnts = zeros(Int, 3)
        uniformhist!(cnts, [0.5, 1.5, 2.5], 0:1:3)
        @test cnts == [1, 1, 1]

        # Accumulation (call twice without zeroing)
        uniformhist!(cnts, [0.5, 1.5, 2.5], 0:1:3)
        @test cnts == [2, 2, 2]

        # Wrong-size error
        @test_throws Exception uniformhist!(zeros(Int, 2), [0.5], 0:1:3)
    end

    @testset "centered_basis" begin
        # Odd n (symmetric)
        cb3 = centered_basis(3)
        @test collect(cb3) ≈ [-1.0, 0.0, 1.0]

        # Even n
        cb4 = centered_basis(4)
        @test length(cb4) == 4

        # n=1
        @test collect(centered_basis(1)) ≈ [0.0]
    end

    @testset "stepsize" begin
        # StepRangeLen
        @test stepsize(0.0:0.5:2.0) == 0.5

        # UnitRange
        @test stepsize(1:10) == 1

        # Plain vector
        @test stepsize([0.0, 0.25, 0.5, 0.75]) == 0.25
    end
end
