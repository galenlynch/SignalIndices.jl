@testset "types" begin
    @testset "div_type" begin
        # Scalar types
        @test div_type(Int) == Float64
        @test div_type(Int, Int) == Float64
        @test div_type(Float32) == Float32
        @test div_type(Float32, Int) == Float32
        @test div_type(Float32, Float64) == Float64
        @test div_type(2, 5) == Float64
        @test div_type(1.0f0, 1) == Float32

        # Array types
        @test div_type(Vector{Int}) == Vector{Float64}
        @test div_type(Vector{Float32}) == Vector{Float32}
        @test div_type(Matrix{Int}) == Matrix{Float64}
        @test div_type(Matrix{Float64}) == Matrix{Float64}
    end
end
