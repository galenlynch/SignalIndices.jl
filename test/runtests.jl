using SignalIndices
using Test

@testset "SignalIndices.jl" begin
    include("test_types.jl")
    include("test_index_time.jl")
    include("test_bin_window.jl")
    include("test_signal.jl")
    include("test_array.jl")
    include("test_pairwise.jl")
    include("test_filtering.jl")
end
