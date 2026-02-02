```@meta
CurrentModule = SignalIndices
```

# [Usage Guide](@id guide)

Functions are grouped by category. Each entry links to the full docstring in the [API Reference](@ref api).

## Index/Time Conversion

Convert between 1-based sample indices and time values for regularly sampled data.

[`ndx_to_t`](@ref) converts an index (or array of indices) to time:

```julia
fs = 30_000.0
ndx_to_t(15001, fs)        # 0.5
ndx_to_t(15001, fs, 1.0)   # 1.5 (with start_t offset)
```

[`t_to_ndx`](@ref) converts a time to the first sample index at or after that time.
[`t_to_last_ndx`](@ref) returns the last sample at or before the time.
[`t_sup_to_ndx`](@ref) returns the last sample strictly before the time.

```julia
t_to_ndx(0.5, fs)           # 15001
t_to_last_ndx(0.5, fs)      # 15001
t_sup_to_ndx(0.5, fs)       # 15000
```

[`n_ndx`](@ref) counts the number of indices in a range. [`ndx_offset`](@ref) computes
the end index for a given start and number of points. [`ndx_wrap`](@ref) wraps an index
using modular arithmetic.

[`clip_ndx`](@ref) clamps an index to `1:length`. [`clip_ndx_deviance`](@ref) also
returns how far the index was clipped.

```julia
clip_ndx(0, 100)              # 1
clip_ndx(150, 100)            # 100
clip_ndx_deviance(150, 100)   # (100, -50)
```

[`duration`](@ref) and [`time_interval`](@ref) compute timing information for a signal:

```julia
duration(30001, fs)            # 1.0 (seconds)
time_interval(30001, fs, 0.0)  # (0.0, 1.0)
```

## Bin and Window Operations

[`bin_bounds`](@ref) returns the `(start, stop)` indices for a bin number and bin size.
[`bin_center`](@ref) returns the center index.

```julia
bin_bounds(2, 10)     # (11, 20)
bin_center(2, 10)     # 15.5
```

[`expand_selection`](@ref) widens an index range by a given amount, clamped to `[1, imax]`:

```julia
expand_selection(5, 10, 100, 3)  # (2, 13)
```

[`make_slice_idx`](@ref) and [`make_expand_idx`](@ref) construct index tuples for
programmatic slicing and broadcasting along arbitrary dimensions.

[`view_trailing_slice`](@ref) returns a view sliced along the last dimension:

```julia
A = rand(3, 4, 5)
v = view_trailing_slice(A, 2)  # equivalent to view(A, :, :, 2)
```

## Signal Processing

[`moving_sum`](@ref) / [`moving_sum!`](@ref) compute a sliding window sum without
zero-padding:

```julia
moving_sum([1, 2, 3, 4, 5], 3)  # [6, 9, 12]
```

[`local_extrema`](@ref) finds indices of local maxima (default) or minima:

```julia
local_extrema([1, 3, 2, 4, 1])        # [2] (local max at index 2)
local_extrema([1, 3, 2, 4, 1], <)     # [3] (local min at index 3)
```

[`find_local_extrema`](@ref) hill-climbs from a starting index to the nearest extremum.

[`find_all_edge_triggers`](@ref) / [`find_first_edge_trigger`](@ref) detect threshold
crossings (rising edges):

```julia
find_all_edge_triggers([0, 0, 1, 1, 0, 1], 1)  # [3, 6]
find_first_edge_trigger([0, 0, 1, 1, 0, 1], 1)  # 3
```

[`thresh_cross`](@ref) finds where a signal crosses upward through a threshold.
[`indices_above_thresh`](@ref) returns contiguous ranges where the signal exceeds a
threshold.

[`filter_no_collisions`](@ref) removes elements from one sorted array that are within
a collision radius of elements in another:

```julia
filter_no_collisions([1, 5, 10, 15], [4, 14], 2)  # [1, 10]
```

[`window_counts`](@ref) counts the number of events within a sliding time window.

[`uniformhist`](@ref) / [`uniformhist!`](@ref) compute fast histograms for
regularly-spaced bins, using reciprocal multiplication for speed:

```julia
uniformhist([0.1, 0.5, 1.2, 2.9], 0.0:1.0:3.0)  # [2, 1, 1]
```

[`centered_basis`](@ref) returns a range of points centered around zero.
[`stepsize`](@ref) extracts the step size from a range or regularly sampled vector.

## Array Utilities

[`weighted_mean`](@ref) and [`weighted_mean_dim`](@ref) compute weighted averages:

```julia
weighted_mean([1.0, 2.0, 3.0], [0.5, 0.3, 0.2])  # 1.7
```

[`simple_summary_stats`](@ref) returns `(mean, std, sem)` in one call.

[`find_closest`](@ref) returns the index of the nearest value:

```julia
find_closest([1.0, 3.0, 5.0], 2.8)  # 2
```

[`find_subseq`](@ref) finds all starting indices where a subsequence occurs.

[`subselect`](@ref) extracts sub-vectors using `(start, stop)` index tuples.

[`rev_view`](@ref) returns a reversed view without copying.

[`trailing_zeros_idx`](@ref) finds the last non-zero index.

[`clipsize!`](@ref) resizes a vector and releases excess memory.

[`copy_length_check`](@ref) and [`copy_length_dest_check`](@ref) validate that a
destination has room for a copy operation.

[`invert_perm`](@ref) returns the inverse of a permutation vector.

## Pairwise Operations

[`pairwise_idxs`](@ref) generates all unique `(i, j)` pairs where `i > j`:

```julia
pairwise_idxs(3)  # [(2,1), (3,1), (3,2)]
```

[`pairwise_idx`](@ref) converts a pair back to its linear index.

[`map_pairwise`](@ref) applies a function to all unique pairs from a vector:

```julia
map_pairwise(-, [10, 20, 30])  # [10, 20, 10]
```

[`imap_product`](@ref) lazily maps a function over the Cartesian product of two
collections.

## Filtering and Comparison

[`filtermap`](@ref) combines `filter` and `map` in a single pass.

[`find_not_unique`](@ref) returns indices of all duplicate elements.

[`skipoftype`](@ref) / [`skipnothing`](@ref) return iterators that skip elements of
a given type:

```julia
sum(skipnothing([1, nothing, 2, nothing, 3]))  # 6
```

[`allsame`](@ref) checks whether all elements (or mapped values) are equal:

```julia
allsame([1, 1, 1])             # true
allsame(length, [1,2], [3,4])  # true
```

[`anyeq`](@ref) checks whether an element exists in a collection.

[`absdiff`](@ref) computes the absolute difference, with branch-free unsigned integer
support.

## Type Utilities

[`div_type`](@ref) determines the appropriate floating-point result type for division.
Float types are preserved; integer types promote to `Float64`.

```julia
div_type(Int)      # Float64
div_type(Float32)  # Float32
```
