module SignalIndices

using IterTools: imap

using SharedArrays: SharedVector

using Statistics: mean, std

import Base:
    eltype,
    mapreduce,
    _mapreduce,
    mapreduce_empty,
    mapreduce_first,
    pairwise_blocksize,
    iterate,
    IteratorSize,
    IteratorEltype,
    SizeUnknown

export
    # Type utilities
    div_type,
    # Index/time conversion
    ndx_to_t,
    t_to_ndx,
    t_to_last_ndx,
    t_sup_to_ndx,
    clip_ndx,
    clip_ndx_deviance,
    n_ndx,
    ndx_offset,
    ndx_wrap,
    duration,
    time_interval,
    bin_bounds,
    bin_center,
    expand_selection,
    copy_length_check,
    copy_length_dest_check,
    make_slice_idx,
    make_expand_idx,
    view_trailing_slice,
    invert_perm,
    # Signal processing
    moving_sum!,
    moving_sum,
    local_extrema,
    find_local_extrema,
    find_all_edge_triggers,
    find_first_edge_trigger,
    thresh_cross,
    indices_above_thresh,
    filter_no_collisions,
    window_counts,
    centered_basis,
    stepsize,
    uniformhist!,
    uniformhist,
    # Array utilities
    weighted_mean,
    weighted_mean_dim,
    simple_summary_stats,
    find_closest,
    subselect,
    find_subseq,
    rev_view,
    pairwise_idxs,
    pairwise_idx,
    map_pairwise,
    imap_product,
    filtermap,
    trailing_zeros_idx,
    clipsize!,
    find_not_unique,
    skipoftype,
    skipnothing,
    allsame,
    anyeq,
    absdiff

"""
    div_type(::Type{N})
    div_type(::Type{N}, ::Type{D})
    div_type(num, den)

Return the result type of dividing values of the given numeric or array type(s).
Float types are preserved; integer types promote to `Float64`.
"""
div_type(::Type{N}) where {N<:AbstractFloat} = N
div_type(::Type{N}) where {N<:Integer} = Float64
div_type(::Type{A}) where {T<:AbstractFloat,N,A<:AbstractArray{T,N}} = Array{T,N}
div_type(::Type{A}) where {T<:Integer,N,A<:AbstractArray{T,N}} = Array{Float64,N}
function div_type(::Type{N}, ::Type{D}) where {N<:Number,D<:Number}
    div_type(promote_type(N, D))
end
div_type(num::N, den::D) where {N<:Number,D<:Number} = div_type(N, D)
div_type(num::N) where {N<:Number} = div_type(N)

# =============================================================================
# Index and time conversion utilities
# =============================================================================

"""
    ndx_to_t(i, fs, start_t=0)

Convert a 1-based index (or array of indices) to its time in a regularly sampled
time series with sample rate `fs` and optional `start_t`.

See also [`t_to_ndx`](@ref), [`t_to_last_ndx`](@ref), [`t_sup_to_ndx`](@ref).
"""
function ndx_to_t end
function ndx_to_t(i::AbstractUnitRange, fs::R, start_t::R = zero(fs)) where {R<:Real}
    (i .- 1) ./ fs .+ start_t
end
function ndx_to_t(i::AbstractUnitRange, fs::R, start_t::R = zero(fs)) where {R<:Integer}
    (i .- 1) ./ fs .+ start_t
end
function ndx_to_t(i::Real, fs::R, start_t::R = zero(fs)) where {R<:Real}
    (i - 1) / fs + start_t
end
function ndx_to_t!(dest::AbstractArray, A::AbstractArray, fs, args...)
    dest .= ndx_to_t.(A, fs, args...)
end
function ndx_to_t(A::AbstractArray, fs::R, start_t::R = zero(fs)) where {R<:Integer}
    ts = Vector{Float64}(undef, length(A))
    ndx_to_t!(ts, A, fs, start_t)
end
function ndx_to_t(
    A::AbstractArray{R},
    fs::R,
    start_t::R = zero(fs),
) where {R<:AbstractFloat}
    ts = Vector{R}(undef, length(A))
    ndx_to_t!(ts, A, fs, start_t)
end
function ndx_to_t(a::AbstractArray{T}, fs::R, start_t::S) where {T,R<:Real,S<:Real}
    P = promote_type(T, R, S)
    promoted = convert.(P, (fs, start_t))
    return ndx_to_t(a, promoted...)
end
function ndx_to_t(i::Real, fs::Real, start_t::Real)
    return ndx_to_t(promote(i, fs, start_t)...)
end


"""
    t_to_ndx(x, fs, start_t=0, T=Int)

Convert a time to its 1-based index in a regularly sampled time series. Returns
the index of the first sample at or after the specified time.

See also [`ndx_to_t`](@ref), [`t_to_last_ndx`](@ref), [`t_sup_to_ndx`](@ref).
"""
function t_to_ndx end
function t_to_ndx(x::AbstractFloat, fs::Real, start_t::Real = zero(fs), T::DataType = Int)
    ceil(T, (x - start_t) * fs - sqrt(eps(typeof(x)))) + one(T)
end
function t_to_ndx(x::Integer, fs::Real, start_t::Real = zero(fs), T::DataType = Int)
    ceil(T, (x - start_t) * fs) + one(T)
end
function t_to_ndx!(dest::AbstractArray, A::AbstractArray, args...)
    dest .= t_to_ndx.(A, args...)
end
function t_to_ndx(a::AbstractArray, fs, start_t = zero(fs), T::DataType = Int)
    dest = similar(a, T)
    t_to_ndx!(dest, a, fs, start_t)
end

"""
    t_to_last_ndx(x, fs, start_t=0, T=Int)

Like [`t_to_ndx`](@ref), but returns the index of the last sample at or before
the specified time.

See also [`ndx_to_t`](@ref), [`t_sup_to_ndx`](@ref).
"""
function t_to_last_ndx(x::Real, fs::Real, start_t::Real = zero(fs), T::DataType = Int)
    floor(T, (x - start_t) * fs) + one(T)
end
function t_to_last_ndx!(dest::AbstractArray, A::AbstractArray, args...)
    dest .= t_to_last_ndx.(A, args...)
end
function t_to_last_ndx(a::AbstractArray, fs, start_t = zero(fs), T::DataType = Int)
    dest = similar(a, T)
    t_to_last_ndx!(dest, a, fs, start_t)
end

"""
    t_sup_to_ndx(x, fs, start_t=0, T=Int)

Convert a time to the index of the last sample strictly before the specified time
in a regularly sampled time series. Equivalent to `t_to_ndx(...) - 1`.

See also [`t_to_ndx`](@ref), [`t_to_last_ndx`](@ref), [`ndx_to_t`](@ref).
"""
function t_sup_to_ndx(x::Real, fs::Real, start_t::Real = zero(fs), T::DataType = Int)
    t_to_ndx(x, fs, start_t, T) - one(T)
end
function t_sup_to_ndx(a::AbstractArray, fs, start_t = zero(fs), T::DataType = Int)
    t_to_ndx(a, fs, start_t, T) .- one(T)
end

"""
    clip_ndx(ndx, l)

Clamp index `ndx` to the valid range `1:l` for an array of length `l`.

See also [`clip_ndx_deviance`](@ref).
"""
function clip_ndx end
@inline clip_ndx(ndx::T, l::T) where {T<:Integer} = clamp(ndx, one(T), l)
@inline clip_ndx(ndx::Integer, l::Integer) = clip_ndx(promote(ndx, l)...)

"""
    clip_ndx_deviance(ndx, l)

Like [`clip_ndx`](@ref), but also returns the deviance (clipped - original) as a
second value.
"""
function clip_ndx_deviance(ndx, l)
    ndxout = clip_ndx(ndx, l)
    return ndxout, ndxout - ndx
end

"""
    n_ndx(start_idx::T, stop_idx::T) where {T<:Integer}

Find the number of indices between `start_idx` and `stop_idx`.
"""
n_ndx(start_idx::T, stop_idx::T) where {T<:Integer} = stop_idx - start_idx + one(T)

"""
    expand_selection(ib, ie, imax, expansion)

Expand the index range `[ib, ie]` by `expansion` in both directions, clamping to
`[1, imax]`. Returns a tuple `(ib_expanded, ie_expanded)`.
"""
function expand_selection(ib::T, ie::T, imax::T, expansion::T) where {T<:Integer}
    ib_exp = max(one(T), ib - expansion)
    ie_exp = min(ie + expansion, imax)
    return (ib_exp, ie_exp)
end

function expand_selection(ib::Integer, ie::Integer, imax::Integer, expansion::Integer)
    expand_selection(promote(ib, ie, imax, expansion)...)
end

"""
    ndx_offset(start_ndx, npt)

Finds the index that will select npt number of elements starting at start_ndx.

If npt is negative, then start_ndx is treated like the end of a range of elements,
and the index required to return npt number of elements is returned.
"""
function ndx_offset(start_ndx::T, npt::T) where {T<:Integer}
    adjust = ifelse(npt < zero(T), one(T), -one(T))
    return start_ndx + npt + adjust
end

"""
    duration(npoints, fs)

Return the duration of a regularly sampled time series with `npoints` samples at
sample rate `fs`. Assumes 1-based indexing, so duration is `(npoints - 1) / fs`.
"""
duration(npoints::Integer, fs::Real) = (npoints - 1) / fs

"""
    time_interval(npoints, fs, start_t=0)
    time_interval(a::AbstractVector, fs, start_t=0)

Return the `(start_t, end_t)` time interval of a regularly sampled time series
with `npoints` samples (or a vector `a`) at sample rate `fs`.
"""
function time_interval end
function time_interval(npoints::Integer, fs::T, start_t::T) where {T<:AbstractFloat}
    return (start_t, start_t + duration(npoints, fs))
end
function time_interval(n::Integer, fs::Real, start_t::Real = zero(fs))
    return time_interval(n, convert(Float64, fs), convert(Float64, start_t))
end
time_interval(a::AbstractVector, args...) = time_interval(length(a), args...)

"""
    bin_bounds(binno, binsize)
    bin_bounds(binno, binsize, max_ndx)

Return the `(start_idx, stop_idx)` 1-based index bounds for bin number `binno`
with `binsize` elements per bin. The three-argument form clamps the result to
`max_ndx`.

See also [`bin_center`](@ref).
"""
function bin_bounds end
# Intended to work with binno as an integer or ranges
# though I can't figure out how to express that
function bin_bounds(
    binno::Union{AbstractUnitRange{T},T},
    binsize::S,
) where {T<:Integer,S<:Integer}
    R = promote_type(T, S)
    idx_start = (binno .- one(R)) .* binsize .+ one(R)
    idx_stop = idx_start .+ binsize .- one(R)
    return (idx_start, idx_stop)
end
function bin_bounds(binno::Real, binsize::Real, max_ndx::Real)
    bounds = bin_bounds(binno, binsize)
    clipped_bounds = min.(bounds, max_ndx)
    return clipped_bounds
end

"""
    bin_center(idxs::NTuple{2})
    bin_center(binno, binsize, ...)

Return the center index of a bin, either from a `(start, stop)` tuple or by
computing [`bin_bounds`](@ref) first.

See also [`bin_bounds`](@ref).
"""
function bin_center end
bin_center(idxs::NTuple{2,<:Real}) = mean(idxs)
bin_center(i::Real, args...) = bin_center(bin_bounds(i, args...))
bin_center(rs::NTuple{2,R}) where {R<:AbstractRange} = (rs[1] + rs[2]) / 2
bin_center(r::AbstractRange, binsize::Real) = bin_center(bin_bounds(r, binsize))
function bin_center!(
    dest::AbstractArray{<:AbstractFloat},
    a::AbstractArray{<:NTuple{2,<:Real}},
)
    dest .= bin_center.(a)
end
function bin_center(a::AbstractArray{<:NTuple{2,F}}) where {F<:AbstractFloat}
    dest = similar(a, F)
    bin_center!(dest, a)
end
function bin_center(a::AbstractArray{<:NTuple{2,<:Integer}})
    dest = similar(a, Float64)
    bin_center!(dest, a)
end

"""
    view_trailing_slice(a::AbstractArray{<:Any,N}, idx)

Return a `view` of `a` sliced along its last dimension at `idx`, with all
leading dimensions selected via `:`. For a 3D array, this is equivalent to
`view(a, :, :, idx)`.
"""
@generated function view_trailing_slice(
    a::AbstractArray{<:Any,N},
    idx::T,
) where {N,T<:Union{Integer,OrdinalRange{<:Integer}}}
    view_trailing_slice_impl(a)
end

function view_trailing_slice_impl(a::Type{<:AbstractArray{<:Any,N}}) where {N}
    exprargs = Vector{Any}(undef, N + 2)
    exprargs[1] = :view
    exprargs[2] = :a
    exprargs[3:(end-1)] .= Ref(:(Colon()))
    exprargs[end] = :idx
    Expr(:call, exprargs...)
end

"""
    make_slice_idx(ndims, dimno, idx)

Construct a tuple of indices that selects `idx` along dimension `dimno` and `:` along
all other dimensions. Useful for programmatic slicing of N-dimensional arrays.
"""
function make_slice_idx(
    ndims::Integer,
    dimno::Integer,
    idx::T,
) where {T<:Union{Integer,OrdinalRange{<:Integer}}}
    idxes = Array{Union{Colon,T}}(undef, ndims)
    idxes .= Colon()
    idxes[dimno] = idx
    return (idxes...,)
end

"""
    make_expand_idx(ndims, dimno)

Construct an index tuple that selects `:` along dimension `dimno` and `1` along
all other dimensions. Useful for broadcasting a vector along one dimension of an
N-dimensional array.
"""
function make_expand_idx(ndims::Integer, dimno::Integer)
    idxes = Array{Union{Colon,Int}}(undef, ndims)
    idxes[:] .= 1
    idxes[dimno] = Colon()
    return (idxes...,)
end

"""
    copy_length_check(n_dest, n_source, d_off=1, s_off=1, n=n_ndx(s_off, n_source))
    copy_length_check(dest::AbstractArray, source::AbstractArray, ...)
    copy_length_check(dest::AbstractArray, d_off, source::AbstractArray, s_off, ...)

Return `true` if `dest` has enough room to accept `n` elements copied from
`source` at the given offsets. The three-argument array form interleaves each
offset with its array: `(dest, d_off, source, s_off, ...)`.
"""
function copy_length_check end

"""
    copy_length_dest_check(n_dest, d_off, n)

Return `true` if a destination of length `n_dest` can accept `n` elements
starting at offset `d_off`.
"""
function copy_length_dest_check(n_dest::Integer, d_off::Integer, n::Integer)
    n_dest >= d_off && n_ndx(d_off, n_dest) >= n
end
function copy_length_check(
    n_dest::Integer,
    n_source::Integer,
    d_off::Integer = 1,
    s_off::Integer = 1,
    n::Integer = n_ndx(s_off, n_source),
)
    source_ok = n_source >= s_off && n_ndx(s_off, n_source) >= n
    source_ok && copy_length_dest_check(n_dest, d_off, n)
end

function copy_length_check(dest::AbstractArray, source::AbstractArray, args...)
    copy_length_check(length(dest), length(source), args...)
end

function copy_length_check(
    dest::AbstractArray,
    d_off::Integer,
    source::AbstractArray,
    s_off::Integer,
    args...,
)
    copy_length_check(dest, source, d_off, s_off, args...)
end

"""
    ndx_wrap(i, max_ndx)

Wrap a 1-based index `i` into the range `1:max_ndx` using modular arithmetic.

# Examples
```jldoctest
julia> SignalIndices.ndx_wrap(5, 3)
2

julia> SignalIndices.ndx_wrap(3, 3)
3
```
"""
function ndx_wrap(i::T, max_ndx::Integer) where {T<:Integer}
    mod(i - one(T), T(max_ndx)) + one(T)
end

"""
    invert_perm(p)

Return the inverse of permutation vector `p`, such that `invert_perm(p)[p[i]] == i`.
"""
function invert_perm(p)
    ip = similar(p)
    for i = 1:length(p)
        ip[p[i]] = i
    end
    ip
end

# =============================================================================
# Signal processing and array utilities
# =============================================================================

"""
    weighted_mean_dim(summand, weights, dim=N, total_weight=sum(weights))

Compute the weighted mean of `summand` along dimension `dim` using `weights`
(a vector whose length matches `size(summand, dim)`). Returns an array with
dimension `dim` reduced to size 1.

See also [`weighted_mean`](@ref).
"""
function weighted_mean_dim(
    summand::AbstractArray{E,N},
    weights::AbstractVector{T},
    dim::Integer = N,
    total_weight::T = sum(weights),
) where {E<:Number,N,T<:Number}
    F = promote_type(E, T)

    dims = collect(size(summand))
    dims[dim] = 1

    outdims = (dims...,)
    reduced = zeros(F, outdims)
    for i = 1:size(summand, dim)
        slice = selectdim(summand, dim, i)
        reduced .+= reshape(slice, outdims) .* weights[i]
    end
    reduced .= reduced ./ total_weight
    return reduced
end

"""
    weighted_mean(summand, weights, total_weight=sum(weights))

Compute the weighted mean of `summand` using element-wise `weights`. Both arrays
must have the same size. Returns a scalar.

See also [`weighted_mean_dim`](@ref).
"""
function weighted_mean(
    summand::AbstractArray{E,N},
    weights::AbstractArray{T,N},
    total_weight::T = sum(weights),
) where {E<:Number,N,T<:Number}
    if size(summand) != size(weights)
        throw(ArgumentError("Sizes are not the same"))
    end
    F = promote_type(E, T)
    reduced = zero(F)
    for i in eachindex(summand)
        reduced += summand[i] * weights[i]
    end
    return reduced / total_weight
end

function weighted_mean(
    summand::AbstractArray{E,N},
    weights::AbstractArray{T,N},
    total_weight::T = sum(weights),
) where {G<:Number,E<:AbstractArray{G},N,T<:Number}
    # assumes elements of summand have the same size
    if size(summand) != size(weights)
        throw(ArgumentError("Sizes are not the same"))
    end
    if isempty(summand)
        throw(ArgumentError("summand is empty"))
    end
    F = promote_type(G, T)
    reduced = zeros(F, size(summand[1]))
    for i in eachindex(summand)
        reduced += summand[i] * weights[i]
    end
    return reduced / total_weight
end

"""
    local_extrema(s, comp = >)

Return indices of local extrema in vector `s`. An index `i` is a local extremum
when `comp(s[i], s[i+1])` is true and `comp(s[i-1], s[i])` is false (i.e., a
direction change). Use `>` for maxima (default), `<` for minima.

Note: plateau-to-descent transitions are reported as extrema. For example,
`local_extrema([1, 2, 2, 1])` returns `[3]` because `2 > 1` but `!(2 > 2)`.
"""
function local_extrema(s::AbstractVector, comp::Function = >)
    ns = length(s)
    idxes = Vector{Int}(undef, div(ns, 2))
    out_i = 0
    if ns > 2
        @inbounds last_comp = comp(s[1], s[2])
        for i = 2:(ns-1)
            @inbounds this_comp = comp(s[i], s[i+1])
            if this_comp && !last_comp
                out_i += 1
                @inbounds idxes[out_i] = i
            end
            last_comp = this_comp
        end
    end
    resize!(idxes, out_i)
    idxes
end

"""
    rev_view(a::AbstractVector)

Return a reversed view of vector `a` without copying.
"""
rev_view(a::AbstractVector) = @view a[end:-1:1]

"""
    pairwise_idxs(n::Integer) -> Vector{NTuple{2, Int}}

Returns the possible combinations of `1:n` indices excluding self-pairs.
The first index is always greater than the second, allowing for easy
subtraction of ordered lists.

# Examples
```julia-repl
julia> pairwise_idxs(3)
3-element Vector{Tuple{Int64,Int64}}:
 (2, 1)
 (3, 1)
 (3, 2)
```
"""
function pairwise_idxs(n::Integer)
    idxs = Vector{NTuple{2,Int}}(undef, convert(Int, n * (n - 1) / 2))
    offset = 0
    @inbounds for i = 1:(n-1)
        @simd for j = 1:(n-i)
            idxs[offset+j] = (i + j, i)
        end
        offset += n - i
    end
    idxs
end

"n is the number of elements in the INPUT of pairwise diff"
_pairwise_idx(i, j, n) = i - j + div((j - 1) * (2 * (n - 1) - (j - 2)), 2)

"""
    pairwise_idx(i, j, n)

Return the linear index into the output of [`pairwise_idxs`](@ref) for the pair
`(i, j)` where `n` is the number of input elements. Throws `ArgumentError` if
`i == j`.

See also [`pairwise_idxs`](@ref), [`map_pairwise`](@ref).
"""
function pairwise_idx(i::T, j::T, n) where {T}
    i == j && throw(ArgumentError("invalid indices"))
    lower, upper = ifelse(i < j, (i, j), (j, i))
    _pairwise_idx(upper, lower, n)
end

"""
    map_pairwise(f, a, R=eltype(a))
    map_pairwise(f, as, bs, R=eltype(as))

Apply `f` to all unique pairs from vector `a` (single-argument form) or to the
Cartesian product of `as` and `bs` (two-argument form). The single-argument form
returns a vector ordered consistently with [`pairwise_idxs`](@ref); the
two-argument form returns a matrix.

See also [`pairwise_idxs`](@ref), [`pairwise_idx`](@ref).
"""
function map_pairwise(f::F, a::AbstractVector{T}, ::Type{R} = T) where {F,T,R}
    n = length(a)
    out = Vector{R}(undef, convert(Int, n * (n - 1) / 2))
    offset = 0
    for i = 1:(n-1)
        @inbounds @simd for j = 1:(n-i)
            out[offset+j] = f(a[i+j], a[i])
        end
        offset += n - i
    end
    out
end

function map_pairwise(
    f::F,
    as::AbstractVector{T},
    bs::AbstractVector,
    ::Type{R} = T,
) where {F,T,R}
    na = length(as)
    nb = length(bs)
    out = Matrix{R}(undef, nb, na)
    @inbounds @simd for i = 1:na
        for j = 1:nb
            out[j, i] = f(as[i], bs[j])
        end
    end
    out
end

"""
    imap_product(f, as, bs)

Return a lazy iterator that applies `f(a, b)` to every element of the Cartesian
product of `as` and `bs`.
"""
imap_product(f, as, bs) = imap(x -> @inbounds(f(x[1], x[2])), Iterators.product(as, bs))

"""
    find_subseq(subseq, seq)

Return a vector of starting indices where `subseq` occurs in `seq`. Uses
`isequal` for element comparison.
"""
function find_subseq(subseq, seq)
    nsub = length(subseq)
    nsub > 0 || throw(ArgumentError("subseq cannot be empty"))
    p = isequal(subseq[1])
    if nsub == 1
        return findall(p, seq)
    end
    nseq = length(seq)
    nsub > nseq && return Vector{Int}()
    max_idx = nseq - nsub + 1
    imatch = Vector{Int}(undef, max_idx)
    nmatch = 0
    idx = 1
    while (idx = findnext(p, seq, idx)) !== nothing
        idx > max_idx && break
        ismatch = true
        @inbounds for i = 2:nsub
            if seq[idx+i-1] != subseq[i]
                ismatch = false
                break
            end
        end
        if ismatch
            nmatch += 1
            imatch[nmatch] = idx
        end
        idx += 1
    end
    resize!(imatch, nmatch)
    imatch
end

"""
    subselect(base_vec, idx_tup_vec, outtype=...)

Extract sub-vectors from `base_vec` using a vector of `(start, stop)` index
tuples. Returns a `Vector{outtype}` where each element is the slice
`base_vec[start:stop]`.
"""
function subselect(
    base_vec,
    idx_tup_vec::AbstractVector{<:NTuple{2}},
    outtype::Type{T} = ifelse(
        base_vec isa AbstractVector,
        typeof(base_vec),
        Vector{eltype(base_vec)},
    ),
) where {T<:AbstractVector}
    nout = length(idx_tup_vec)
    out = Vector{T}(undef, nout)
    for i = 1:nout
        ib, ie = idx_tup_vec[i]
        out[i] = convert(outtype, view(base_vec, ib:ie))
    end
    out
end

function subselect(
    base_vec,
    idx_tup_vec::AbstractVector{<:NTuple{2}},
    outtype::Type{T},
) where {T<:SharedVector}
    nout = length(idx_tup_vec)
    out = Vector{T}(undef, nout)
    for i = 1:nout
        ib, ie = idx_tup_vec[i]
        out[i] = outtype(base_vec[ib:ie])
    end
    out
end

"""
    simple_summary_stats(a)

Return `(mean, std, sem)` for array `a`, where `sem` is the standard error of
the mean.
"""
function simple_summary_stats(a::AbstractArray)
    m = mean(a)
    s = std(a)
    sem = s / sqrt(length(a))
    m, s, sem
end

"""
    find_closest(a, target, ...)

Find the index of the element in `a` that minimizes the absolute difference from
the target.
"""
function find_closest end

find_closest(arr::AbstractVector, target) = argmin(abs.(arr .- target))

function find_closest(arr::AbstractVector, target, eligible::AbstractVector)
    elig_ndxs = findall(eligible)
    rel_ndx = find_closest(arr[eligible], target)
    elig_ndxs[rel_ndx]
end

function find_closest(f::Function, arr::AbstractVector, target, args...)
    find_closest(map(f, arr), f(target), args...)
end

"""
    skipoftype(::Type{T}, itr)
Return an iterator over the elements in `itr` skipping values of type `T`.
Use `collect` to obtain an `Array` containing the non-`T` values in
`itr`. Note that even if `itr` is a multidimensional array, the result will always
be a `Vector` since it is not possible to remove nothings while preserving dimensions
of the input.
# Examples
```jldoctest
julia> sum(SignalIndices.skipoftype(nothing, [1, nothing, 2]))
3
julia> collect(SignalIndices.skipoftype(nothing, [1, nothing, 2]))
2-element Vector{Int64}:
 1
 2
julia> collect(SignalIndices.skipoftype(nothing, [1 nothing; 2 nothing]))
2-element Vector{Int64}:
 1
 2
```
"""
skipoftype(::Type{T}, itr::A) where {T,A} = SkipOfType{T,A}(itr)
skipoftype(::T, itr) where {T} = skipoftype(T, itr)

struct SkipOfType{T,A}
    x::A
end

IteratorSize(::Type{<:SkipOfType}) = SizeUnknown()
IteratorEltype(::Type{SkipOfType{T,A}}) where {T,A} = IteratorEltype(A)
eltype(::Type{SkipOfType{T,A}}) where {T,A} = union_poptype(T, eltype(A))

function iterate(itr::SkipOfType{T}, state...) where {T}
    y = iterate(itr.x, state...)
    isnothing(y) && return nothing
    item, state = y
    while item isa T
        y = iterate(itr.x, state)
        isnothing(y) && return nothing
        item, state = y
    end
    item, state
end

# Optimized mapreduce implementation
# The generic method is faster when !(eltype(A) >: Nothing) since it does not need
# additional loops to identify the two first non-nothing values of each block
function mapreduce(f, op, itr::SkipOfType{T,<:AbstractArray}) where {T}
    _mapreduce(f, op, IndexStyle(itr.x), eltype(itr.x) >: T ? itr : itr.x)
end

function _mapreduce(f, op, ::IndexLinear, itr::SkipOfType{T,<:AbstractArray}) where {T}
    A = itr.x
    local ai
    inds = LinearIndices(A)
    i = first(inds)
    ilast = last(inds)
    while i <= ilast
        @inbounds ai = A[i]
        ai isa T || break
        i += 1
    end
    i > ilast && return mapreduce_empty(f, op, eltype(itr))
    a1 = ai
    i += 1
    while i <= ilast
        @inbounds ai = A[i]
        ai isa T || break
        i += 1
    end
    i > ilast && return mapreduce_first(f, op, a1)
    # We know A contains at least two non-nothing entries: the result cannot be nothing
    something(mapreduce_impl(f, op, itr, first(inds), last(inds)))
end

_mapreduce(f, op, ::IndexCartesian, itr::SkipOfType) = mapfoldl(f, op, itr)

mapreduce_impl(f, op, A::SkipOfType, ifirst::Integer, ilast::Integer) =
    mapreduce_impl(f, op, A, ifirst, ilast, pairwise_blocksize(f, op))

# Returns nothing when the input contains only nothing values
@noinline function mapreduce_impl(
    f,
    op,
    itr::SkipOfType{T,<:AbstractArray},
    ifirst::Integer,
    ilast::Integer,
    blksize::Int,
) where {T}
    A = itr.x
    if ifirst == ilast
        @inbounds a1 = A[ifirst]
        if a1 isa T
            return nothing
        else
            return Some(mapreduce_first(f, op, a1))
        end
    elseif ifirst + blksize > ilast
        # sequential portion
        local ai
        i = ifirst
        while i <= ilast
            @inbounds ai = A[i]
            ai isa T || break
            i += 1
        end
        i > ilast && return nothing
        a1 = ai::eltype(itr)
        i += 1
        while i <= ilast
            @inbounds ai = A[i]
            ai isa T || break
            i += 1
        end
        i > ilast && return Some(mapreduce_first(f, op, a1))
        a2 = ai::eltype(itr)
        i += 1
        v = op(f(a1), f(a2))
        @simd for i = i:ilast
            @inbounds ai = A[i]
            if !(ai isa T)
                v = op(v, f(ai))
            end
        end
        return Some(v)
    else
        # pairwise portion
        imid = (ifirst + ilast) >> 1
        v1 = mapreduce_impl(f, op, itr, ifirst, imid, blksize)
        v2 = mapreduce_impl(f, op, itr, imid+1, ilast, blksize)
        if isnothing(v1) && isnothing(v2)
            return nothing
        elseif isnothing(v1)
            return v2
        elseif isnothing(v2)
            return v1
        else
            return Some(op(something(v1), something(v2)))
        end
    end
end

_union_poptype(::Type{T}, ::Type{Union{T,S}}) where {T,S} = S
_union_poptype(::Type{T}, ::Type{T}) where {T} = Union{}

# Necessary for cases like _union_poptype(Union{A,B}, Union{A,C})
union_poptype(::Type{T}, ::Type{S}) where {T,S} = _union_poptype(T, Union{T,S})

"""
    skipnothing(itr)

Return an iterator that skips `nothing` values in `itr`. Equivalent to
`skipoftype(Nothing, itr)`.

See also [`skipoftype`](@ref).
"""
skipnothing(itr) = skipoftype(Nothing, itr)

# Does not check input lengths
function _moving_sum!(out, s, nav, nout)
    if nav <= 1
        copyto!(out, 1, s, 1, nout)
    elseif nout > 0
        @inbounds out[1] = 0
        @inbounds @simd for i = 1:nav
            out[1] += s[i]
        end
        @inbounds for i = 2:nout
            # The only thing that changes is the first and last part of the window
            out[i] = out[i-1] + s[i+nav-1] - s[i-1]
        end
    end
    out
end

"""
    moving_sum!(out, s, nav)

Sum `s` in a sliding window of `nav` points, placing the result into `out`.
The length of `out` should be `max(length(s) - nav + 1, 0)` if `nav > 0`, or
`length(s)` otherwise. When `nav` is 0 or 1 the result is a copy of `s`.

Does not zero-pad.
"""
function moving_sum!(out, s, nav)
    nav < 0 && throw(ArgumentError("nav must be non-negative, got $nav"))
    nin = length(s)
    nout = length(out)
    if nout != ifelse(nav == 0, nin, max(nin - nav + 1, min(nin, 1)))
        throw(ArgumentError("out is not the right size"))
    end
    _moving_sum!(out, s, min(nav, nin), nout)
end

"""
    moving_sum(s, nav)

Same as [`moving_sum!`](@ref), but returns a new array. When `nav` is 0 or 1 the
result is a copy of `s`.
"""
function moving_sum(s::AbstractVector, nav::Integer)
    nav < 0 && throw(ArgumentError("nav must be non-negative, got $nav"))
    nin = length(s)
    nout = ifelse(nav == 0, nin, max(nin - nav + 1, min(nin, 1)))
    _moving_sum!(similar(s, nout), s, min(nav, nin), nout)
end

"""
    trailing_zeros_idx(arr)

return last index that is not zero
"""
function trailing_zeros_idx(arr)
    l = length(arr)
    last_idx = l
    while last_idx > 0 && arr[last_idx] == 0
        last_idx -= 1
    end
    last_idx
end

"""
    thresh_cross(arr, thresh, comp = <)

Return indices where `arr` crosses upward through `thresh` according to `comp`.
An index `i+1` is returned when `comp(arr[i], thresh)` is true and
`comp(arr[i+1], thresh)` is false (i.e., the signal leaves the `comp` region).
"""
function thresh_cross(arr, thresh, comp = <)
    l = length(arr)
    idx_cross = Vector{Int}(undef, div(l, 2))
    out_no = 0
    for i = 1:(l-1)
        if comp(arr[i], thresh) & (!comp(arr[i+1], thresh))
            out_no += 1
            idx_cross[out_no] = i + 1
        end
    end
    resize!(idx_cross, out_no)
    idx_cross
end

"""
    centered_basis(n_point)

Return a range of `n_point` values centered around zero (e.g., `[-1.0, 0.0, 1.0]`
for `n_point=3`).
"""
centered_basis(n_point) = (0:(n_point-1)) .- (n_point - 1) / 2

"""
    stepsize(r)

Return the step size of a range or regularly sampled vector. For a `UnitRange`,
returns `1`. For a generic vector, returns `a[2] - a[1]` (assumes regular spacing).
"""
stepsize(r::StepRangeLen) = Float64(r.step)
stepsize(::UnitRange) = 1
stepsize(a::AbstractVector) = a[2] - a[1]

@inline Base.@propagate_inbounds function _uniformhist_push!(cnts, x, nbin, m, offset)
    binndx = floor(Int, muladd(m, x, offset)) + 1
    inbounds = (binndx > 0) & (binndx <= nbin)
    trunc_ndx = ifelse(inbounds, binndx, 1)
    cnts[trunc_ndx] += inbounds
end

function _uniformhist!(cnts, xs, first, nbin::Integer, step)
    # Approximating division with multiplication of inverse is 20x faster
    m = 1 / step
    offset = -m * first
    for x in xs
        @inbounds _uniformhist_push!(cnts, x, nbin, m, offset)
    end
    cnts
end

_uniformhist!(cnts, xs, r) = _uniformhist!(cnts, xs, first(r), length(r) - 1, stepsize(r))

"""
    uniformhist!(cnts, xs, r)

Compute a histogram of `xs` into pre-allocated count vector `cnts`, using the bin edges
defined by range `r`. Bins are left-inclusive: a value `x` falls into bin `i` when
`r[i] <= x < r[i+1]`.

`r` must have uniform step size (e.g. a `StepRangeLen` or `UnitRange`) and at least 2
elements. The length of `cnts` must equal `length(r) - 1`. Values outside the range are
silently ignored.

`cnts` is not zeroed before accumulation, so it must be initialized (e.g. with `zeros`).
This also means `uniformhist!` can be called repeatedly to accumulate counts from multiple
datasets.

Uses multiplication by the reciprocal of the step size instead of division for a ~20x
speedup over naive binning.

See also [`uniformhist`](@ref).

# Examples
```jldoctest
julia> cnts = zeros(Int, 3);

julia> uniformhist!(cnts, [0.1, 0.5, 1.2, 2.9], 0.0:1.0:3.0)
3-element Vector{Int64}:
 2
 1
 1
```
"""
function uniformhist!(cnts, xs, r)
    length(cnts) == length(r) - 1 || error("cnts must be length length(r) - 1")
    _uniformhist!(cnts, xs, r)
end

"""
    uniformhist([::Type{T} = Int,] xs, r) where T

Compute a histogram of `xs` using the bin edges defined by range `r`, returning a new
vector of counts with element type `T` (default `Int`). Bins are left-inclusive: a value
`x` falls into bin `i` when `r[i] <= x < r[i+1]`.

`r` must have uniform step size (e.g. a `StepRangeLen` or `UnitRange`). The returned
vector has length `length(r) - 1`. Values outside the range are silently ignored.

See also [`uniformhist!`](@ref).

# Examples
```jldoctest
julia> uniformhist([0.1, 0.5, 1.2, 2.9], 0.0:1.0:3.0)
3-element Vector{Int64}:
 2
 1
 1

julia> uniformhist(Float64, [1, 1, 2, 3, 3, 3], 1:4)
3-element Vector{Float64}:
 2.0
 1.0
 3.0
```
"""
uniformhist(::Type{T}, xs, r) where {T} = _uniformhist!(zeros(T, length(r) - 1), xs, r)
uniformhist(xs, r) = uniformhist(Int, xs, r)

"""
    find_local_extrema(sig, start_ndx=div(length(sig), 2); findmax=true, right_on_ties=true)

Starting from `start_ndx`, hill-climb through `sig` to find a local maximum
(or minimum if `findmax=false`). When tied, moves right if `right_on_ties=true`.
Skips `missing` and `NaN` values. Returns the index of the found extremum.
"""
function find_local_extrema(
    sig::AbstractVector,
    start_ndx::Integer = div(length(sig), 2);
    findmax::Bool = true,
    right_on_ties::Bool = true,
)
    sigl = length(sig)
    sigl < 2 && return start_ndx
    checkbounds(sig, start_ndx)
    comp = ifelse(findmax, >=, <=)
    bias = ifelse(right_on_ties, 1, -1)
    if ismissing(sig[start_ndx]) || isnan(sig[start_ndx])
        newstart = findfirst(x -> !(ismissing(x) | isnan(x)), sig)
        isnothing(newstart) && error("None valid")
        start_ndx = newstart
    end
    search_ndx = start_ndx
    iterno = 0
    maxiter = 2 * sigl
    @inbounds while iterno < maxiter
        notleft =
            search_ndx == 1 ||
            ismissing(sig[search_ndx-1]) ||
            isnan(sig[search_ndx-1]) ||
            comp(sig[search_ndx], sig[search_ndx-1])
        notright =
            search_ndx == sigl ||
            ismissing(sig[search_ndx+1]) ||
            isnan(sig[search_ndx+1]) ||
            comp(sig[search_ndx], sig[search_ndx+1])
        if notleft & notright
            return search_ndx
        else
            search_ndx += ifelse(notleft, 1, ifelse(notright, -1, bias))
        end
        iterno += 1
    end
    error("Did not converge")
end

"""
    filter_no_collisions(as, bs, coll_rad)

Filter elements of `as` to keep elements that are not within `coll_rad` of any
element in `bs`. Assumes both are sorted.
"""
function filter_no_collisions(as, bs, coll_rad)
    out = similar(as)
    outno = 0
    nb = length(bs)
    na = length(as)
    ib = 1
    for (i, a) in enumerate(as)
        # Skip over bs that are too far back to matter
        while ib <= nb && a - bs[ib] > coll_rad
            ib += 1
        end
        if ib > nb
            # No times to avoid, push the rest of as into out
            nremainder = na - i + 1
            copyto!(out, outno + 1, as, i, nremainder)
            outno += nremainder
            break
        end
        # Only keep elements of as that do not collide with bs
        # If the next element of bs does not collide, none of the others will
        if abs(a - bs[ib]) > coll_rad
            outno += 1
            out[outno] = a
        end
    end
    resize!(out, outno)
    out
end

"""
    window_counts(ts, window_dur)

For each event in `ts`, count the number of events in `ts` that are in the
range of `ts[i]` and `ts[i] + window_dur`. Assumes `ts` is sorted, and
that elements of `ts` are unique.
"""
function window_counts(ts, window_dur)
    cnts = similar(ts, Int)
    for (i, t) in enumerate(ts)
        se = searchsortedlast(ts, t + window_dur)
        cnts[i] = se - i + 1
    end
    cnts
end

"`window_counts` in a certain range`"
function window_counts(ts, window_dur, tb, te)
    ib = searchsortedfirst(ts, tb)
    ie = searchsortedlast(ts, te)
    subset = view(ts, ib:ie)
    cnts = window_counts(subset, window_dur)
    cnts, ib
end

"""
    filtermap(p::Function, f::Function, xs)

Filters the input vector, then maps the remaining values. For each element
of `xs` which predicate function `p` returns true for, use mapping function `f`
to transform the result."""
function filtermap(p::Function, f::Function, xs::AbstractVector)
    if isempty(xs)
        return similar(xs, Base.promote_op(f, eltype(xs)))
    end
    # Find the first element passing the predicate to determine output type
    first_pass = findfirst(p, xs)
    if isnothing(first_pass)
        return similar(xs, Base.promote_op(f, eltype(xs)), 0)
    end
    @inbounds begin
        m_el1 = f(xs[first_pass])
        out = similar(xs, typeof(m_el1))
        nout = 1
        out[1] = m_el1
        for elno = (first_pass + 1):length(xs)
            if p(xs[elno])
                nout += 1
                out[nout] = f(xs[elno])
            end
        end
    end
    resize!(out, nout)
    out
end

"""
    find_not_unique(a::AbstractArray)

Returns the indices of all redundant elements in a. The time a value is
seen, it is not considered redundant
"""
function find_not_unique(a::AbstractArray{T}) where {T}
    na = length(a)

    # Stores the first seen index, and if the index is a known
    # duplicate
    seen_els = Dict{T,Tuple{Int,Bool}}()
    sizehint!(seen_els, na)
    redundant_ndxs = Vector{Int}(undef, na)
    outno = 0
    for (i, el) in enumerate(a)
        if haskey(seen_els, el)
            (first_ndx, duplicated) = seen_els[el]
            if !duplicated
                seen_els[el] = (first_ndx, true)
                outno += 1
                redundant_ndxs[outno] = first_ndx
            end
            outno += 1
            redundant_ndxs[outno] = i
        else
            seen_els[el] = (i, false)
        end
    end
    resize!(redundant_ndxs, outno)
    redundant_ndxs
end

"""
    clipsize!(a::AbstractVector, n)

Resize vector `a` to length `n` and release excess memory via `sizehint!`.
"""
clipsize!(a::AbstractVector, n::Integer) = sizehint!(resize!(a, n), n)

"""
    find_all_edge_triggers(arr, thr, comp = >=)

Return all indices where `arr` triggers an edge crossing of threshold `thr`.
An edge occurs at index `i` when `comp(arr[i], thr)` is true and
`comp(arr[i-1], thr)` is false.

See also [`find_first_edge_trigger`](@ref).
"""
function find_all_edge_triggers(arr, thr, comp = >=)
    indices = Int[]
    @inbounds for i = 2:length(arr)
        if comp(arr[i], thr) & !comp(arr[i-1], thr)
            push!(indices, i)
        end
    end
    indices
end

"""
    find_first_edge_trigger(arr, thr, comp = >=)

Return the index of the first edge trigger in `arr` crossing `thr`, or `nothing`
if none is found.

See also [`find_all_edge_triggers`](@ref).
"""
function find_first_edge_trigger(arr, thr, comp = >=)
    @inbounds for i = 2:length(arr)
        if comp(arr[i], thr) & !comp(arr[i-1], thr)
            return i
        end
    end
    nothing
end

"""
    indices_above_thresh(arr, thr)

Return a vector of `UnitRange{Int}` spans where `arr[i] >= thr`. Consecutive
above-threshold indices are merged into a single range.
"""
function indices_above_thresh(arr, thr)
    out = Vector{UnitRange{Int}}()
    curr_start = nothing
    lasti = 0
    for (i, el) in enumerate(arr)
        above_thr = el >= thr
        if isnothing(curr_start)
            if above_thr
                curr_start = i
            end
        else
            if !above_thr
                push!(out, curr_start:(i-1))
                curr_start = nothing
            end
        end
        lasti = i
    end
    if !isnothing(curr_start)
        push!(out, curr_start:lasti)
    end
    return out
end

"""
    allsame(a::AbstractArray)
    allsame(first, second, others...)
    allsame(f::Function, first, second, others...)

Return `true` if all arguments (or elements of `a`) are equal. The `f` form
compares `f(x)` values.
"""
allsame(f::Function, first) = true

function allsame(f::Function, first, second, others...)
    f(first) == f(second) && allsame(f, second, others...)
end

allsame(first, args...) = allsame(identity, first, args...)

function allsame(a::AbstractArray)
    isempty(a) && return true
    @inbounds first = a[1]
    for e in @view a[2:end]
        if !isequal(e, first)
            return false
        end
    end
    true
end

"""
    anyeq(el, iter)

Return `true` if any element of `iter` equals `el`.
"""
anyeq(el, iter) = any(a -> a == el, iter)

"""
    absdiff(a, b)

Return the absolute difference `|a - b|`. Uses branch-free subtraction for
unsigned integers to avoid overflow.
"""
@inline absdiff(a::Unsigned, b::Unsigned) = ifelse(a <= b, b - a, a - b)
@inline absdiff(a::Signed, b::Signed) = abs(a - b)
@inline absdiff(a::Real, b::Real) = abs(a - b)

end
