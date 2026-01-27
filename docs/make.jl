using SignalIndices
using Documenter

DocMeta.setdocmeta!(SignalIndices, :DocTestSetup, :(using SignalIndices); recursive=true)

makedocs(;
    modules=[SignalIndices],
    authors="Galen Lynch <galen@galenlynch.com>",
    sitename="SignalIndices.jl",
    format=Documenter.HTML(;
        canonical="https://galenlynch.github.io/SignalIndices.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/galenlynch/SignalIndices.jl",
    devbranch="main",
)
