@testset "filtering and equality" begin
    @testset "filtermap" begin
        @test filtermap(x -> x > 0, x -> x^2, [-1, 2, -3, 4]) == [4, 16]

        # Empty input
        @test filtermap(x -> x > 0, x -> x^2, Int[]) == Int[]

        # Nothing passes predicate
        @test filtermap(x -> x > 100, x -> x^2, [1, 2, 3]) == Int[]
    end

    @testset "find_not_unique" begin
        @test sort(find_not_unique([1, 2, 1, 3, 2])) == [1, 2, 3, 5]
        @test find_not_unique([1, 2, 3]) == []

        # All duplicates
        @test sort(find_not_unique([1, 1, 1])) == [1, 2, 3]

        # Single element
        @test find_not_unique([1]) == []
    end

    @testset "skipoftype" begin
        # Skip Int from mixed
        result = collect(skipoftype(Int, Any[1, "a", 2, "b"]))
        @test result == ["a", "b"]

        # Skip Nothing (same as skipnothing)
        @test collect(skipoftype(Nothing, [1, nothing, 2])) == [1, 2]

        # Collect preserves order
        @test collect(skipoftype(Nothing, [3, nothing, 1, nothing, 2])) == [3, 1, 2]

        # sum/reduce
        @test sum(skipoftype(Nothing, [1, nothing, 2, nothing, 3])) == 6
    end

    @testset "skipnothing" begin
        @test sum(skipnothing([1, nothing, 2, nothing, 3])) == 6
        @test collect(skipnothing([1, nothing, 2])) == [1, 2]

        # All-nothing input (empty result)
        @test collect(skipnothing([nothing, nothing])) == []

        # Matrix input (returns vector)
        @test collect(skipnothing([1 nothing; 2 nothing])) == [1, 2]

        # mapreduce optimization path
        @test mapreduce(identity, +, skipnothing([1, nothing, 2, nothing, 3])) == 6
    end

    @testset "allsame" begin
        A = ones(3, 3)
        @test allsame(A)
        A[1] = 0
        @test !allsame(A)

        @test allsame(1, 1)
        @test !allsame(1, 2)
        @test allsame(1)
        @test allsame(length, (1, 2), (3, 4))

        # Empty array returns true
        @test allsame(Int[])

        # Single element
        @test allsame([42])
    end

    @testset "anyeq" begin
        # Found
        @test anyeq(2, [1, 2, 3])

        # Not found
        @test !anyeq(5, [1, 2, 3])

        # Empty collection
        @test !anyeq(1, Int[])
    end

    @testset "absdiff" begin
        # Unsigned integers (no overflow)
        @test absdiff(UInt(5), UInt(3)) == UInt(2)
        @test absdiff(UInt(3), UInt(5)) == UInt(2)

        # Signed integers
        @test absdiff(5, 3) == 2
        @test absdiff(-3, 5) == 8

        # Floats
        @test absdiff(1.5, 3.0) ≈ 1.5

        # Zero difference
        @test absdiff(7, 7) == 0
    end
end
