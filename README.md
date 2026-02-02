# SignalIndices [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://galenlynch.github.io/SignalIndices.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://galenlynch.github.io/SignalIndices.jl/dev/) [![Build Status](https://github.com/galenlynch/SignalIndices.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/galenlynch/SignalIndices.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/galenlynch/SignalIndices.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/galenlynch/SignalIndices.jl)

Utilities for working with regularly sampled time series data, including index-time conversions, signal processing primitives, and array manipulation functions.

## Quick Example

```julia
using SignalIndices

fs = 30_000.0  # 30 kHz sampling rate

# Convert between sample indices and time
t = ndx_to_t(15001, fs)         # 0.5 s
idx = t_to_ndx(0.5, fs)         # 15001

# Detect threshold crossings in a signal
edges = find_all_edge_triggers(signal, threshold)

# Sliding window sum (no zero-padding) and fast histogram
sums = moving_sum(signal, 100)
counts = uniformhist(spike_times, 0.0:0.001:1.0)
```

## When to Use What

| Package | Focus |
|---------|-------|
| **SignalIndices.jl** | Index/time math, edge detection, and array primitives for regularly sampled data |
| [DSP.jl](https://github.com/JuliaDSP/DSP.jl) | Filters, FFTs, spectral analysis — use *with* SignalIndices for bookkeeping |
| [SampledSignals.jl](https://github.com/JuliaAudio/SampledSignals.jl) | Audio I/O with rich buffer types — audio-specific |
| [TimeSeries.jl](https://github.com/JuliaStats/TimeSeries.jl) | `TimeArray` types for financial/irregular time series — stores timestamps per sample |

## Installation

```julia
# Not yet registered - install from path or URL
using Pkg
Pkg.develop(path="/path/to/SignalIndices")
```

## Documentation

See the [full documentation](https://galenlynch.github.io/SignalIndices.jl/dev/) for a usage guide and API reference.

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
