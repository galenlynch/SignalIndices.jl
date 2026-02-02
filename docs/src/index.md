```@meta
CurrentModule = SignalIndices
```

# SignalIndices.jl

Low-level building blocks for regularly sampled time series — index/time conversion, threshold detection, sliding windows, and array utilities. Designed for neuroscience and electrophysiology workflows where data is sampled at a known frequency and stored in plain Julia arrays.

## Design Principles

- **Plain arrays, not signal types** — works with standard Julia arrays, computing times from indices and sample rates rather than wrapping data in custom types.
- **1-based index convention** — all functions assume Julia's native 1-based indexing. Index 1 corresponds to `start_t`.
- **Allocation-conscious** — in-place variants ([`moving_sum!`](@ref), [`uniformhist!`](@ref)) and pre-sized outputs where possible.
- **Composable primitives** — small functions that chain together rather than monolithic signal processing pipelines.

## Example Workflow

A typical neuroscience analysis: load spike times, convert to indices, detect threshold crossings, and build a histogram.

```julia
using SignalIndices

fs = 30_000.0   # 30 kHz sampling rate
start_t = 0.0

# Convert event times (seconds) to sample indices
spike_indices = t_to_ndx([0.103, 0.207, 0.515], fs, start_t)

# Convert back to verify
spike_times = ndx_to_t(spike_indices, fs, start_t)

# Detect where a filtered signal crosses a threshold
edges = find_all_edge_triggers(filtered_signal, 3.5)

# Smooth with a sliding window, then find local maxima
smoothed = moving_sum(raw_signal, 50)
peaks = local_extrema(smoothed)

# Build a peri-stimulus time histogram
counts = uniformhist(spike_times, 0.0:0.001:1.0)
```

## Next Steps

- [Usage Guide](@ref guide) — categorized function reference with examples
- [API Reference](@ref api) — complete docstrings for every exported function
