# SignalIndices [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://galenlynch.github.io/SignalIndices.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://galenlynch.github.io/SignalIndices.jl/dev/) [![Build Status](https://github.com/galenlynch/SignalIndices.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/galenlynch/SignalIndices.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/galenlynch/SignalIndices.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/galenlynch/SignalIndices.jl)

Utilities for working with regularly sampled time series data, including index-time conversions, signal processing primitives, and array manipulation functions.

## Overview

This package provides low-level building blocks for signal processing workflows, particularly for neuroscience and electrophysiology applications where data is regularly sampled at a known frequency. The focus is on efficient, allocation-conscious implementations suitable for processing large datasets.

## Features

### Index/Time Conversion

Convert between sample indices and time values for regularly sampled data:

```julia
using SignalIndices

fs = 30000.0  # 30 kHz sampling rate
start_t = 0.0

# Convert time to index
idx = t_to_ndx(0.5, fs, start_t)  # First sample at or after 0.5s

# Convert index to time
t = ndx_to_t(15001, fs, start_t)  # Time of sample 15001

# Duration and time intervals
dur = duration(30001, fs)  # Duration spanned by 30001 samples
interval = time_interval(30001, fs, start_t)  # (start_t, end_t) tuple
```

### Signal Processing Primitives

```julia
# Sliding window sum (no zero-padding)
sums = moving_sum(signal, window_size)

# Find threshold crossings (edge triggers)
rising_edges = find_all_edge_triggers(signal, threshold)
first_edge = find_first_edge_trigger(signal, threshold)

# Find local extrema
maxima_indices = local_extrema(signal, >)
minima_indices = local_extrema(signal, <)

# Fast histogram for regularly-spaced bins
counts = glhist(values, bin_edges)
```

### Array Utilities

```julia
# Pairwise operations
idxs = pairwise_idxs(n)  # All (i,j) pairs where i > j
diffs = map_pairwise(-, sorted_times)  # Pairwise differences

# Filtering and searching
idx = find_closest(array, target)
non_colliding = filter_no_collisions(events_a, events_b, min_gap)
duplicates = find_not_unique(array)

# Iterator utilities
result = sum(skipnothing([1, nothing, 2, nothing, 3]))  # 6
```

## Relation to Other Packages

### vs DSP.jl

[DSP.jl](https://github.com/JuliaDSP/DSP.jl) provides comprehensive digital signal processing including filters, FFTs, and spectral analysis. SignalIndices.jl is complementary, focusing on:
- Index/time bookkeeping for sampled data
- Simple threshold and edge detection
- Utilities commonly needed before/after DSP operations

### vs SampledSignals.jl

[SampledSignals.jl](https://github.com/JuliaAudio/SampledSignals.jl) focuses on audio I/O with rich type hierarchies for audio buffers. This package has different goals:
- No audio-specific functionality
- Focus on index/time math and signal primitives
- Designed for neuroscience/electrophysiology workflows

### vs TimeSeries.jl

[TimeSeries.jl](https://github.com/JuliaStats/TimeSeries.jl) provides `TimeArray` types for financial/irregular time series. SignalIndices.jl assumes regular sampling and works with plain arrays, computing times from indices and sampling rates rather than storing timestamps.

## Installation

```julia
# Not yet registered - install from path or URL
using Pkg
Pkg.develop(path="/path/to/SignalIndices")
```

## API Reference

### Type Utilities
- `div_type(T)` - Determine appropriate floating-point type for division

### Index/Time Conversion
- `ndx_to_t(i, fs, start_t)` - Index to time
- `t_to_ndx(t, fs, start_t)` - Time to index (first sample at or after)
- `t_to_last_ndx(t, fs, start_t)` - Time to index (last sample at or before)
- `t_sup_to_ndx(t, fs, start_t)` - Time to index (last sample before)
- `clip_ndx(idx, len)` - Clamp index to valid range
- `n_ndx(start, stop)` - Count indices in range
- `duration(npoints, fs)` - Time span of samples
- `time_interval(npoints, fs, start_t)` - Start/end time tuple

### Bin/Window Operations
- `bin_bounds(binno, binsize)` - Index bounds for a bin
- `bin_center(binno, binsize)` - Center index of a bin
- `expand_selection(ib, ie, imax, expansion)` - Expand index range

### Signal Processing
- `moving_sum(signal, window)` / `moving_sum!(out, signal, window)` - Sliding sum
- `find_all_edge_triggers(arr, threshold, comp)` - Find all threshold crossings
- `find_first_edge_trigger(arr, threshold, comp)` - Find first crossing
- `thresh_cross(arr, threshold, comp)` - Threshold crossing indices
- `indices_above_thresh(arr, threshold)` - Contiguous regions above threshold
- `local_extrema(signal, comp)` - Find local maxima/minima
- `find_local_extrema(signal, start; findmax, right_on_ties)` - Hill-climb to extremum
- `glhist(values, bins)` / `glhist!(counts, values, bins)` - Fast histogram
- `filter_no_collisions(as, bs, radius)` - Remove elements too close to reference
- `window_counts(times, duration)` - Count events in sliding window

### Array Utilities
- `weighted_mean(values, weights)` / `weighted_mean_dim(arr, weights, dim)` - Weighted averages
- `simple_summary_stats(arr)` - Mean, std, SEM tuple
- `find_closest(arr, target)` - Index of nearest value
- `find_subseq(subseq, seq)` - Find subsequence occurrences
- `subselect(vec, index_tuples)` - Extract subvectors by index ranges
- `pairwise_idxs(n)` - Generate all unique pairs
- `pairwise_idx(i, j, n)` - Linear index for pair
- `map_pairwise(f, vec)` - Apply function to all pairs
- `imap_product(f, as, bs)` - Lazy Cartesian product map
- `filtermap(predicate, f, xs)` - Combined filter and map
- `find_not_unique(arr)` - Indices of duplicate elements
- `trailing_zeros_idx(arr)` - Last non-zero index
- `clipsize!(vec, n)` - Resize and shrink capacity
- `rev_view(vec)` - Reversed view
- `skipoftype(T, itr)` / `skipnothing(itr)` - Skip elements by type

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
