# see documentation at https://juliadocs.github.io/Documenter.jl/stable/

using Documenter, SensitivityAnalysis

makedocs(
    modules = [SensitivityAnalysis],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Tamás K. Papp",
    sitename = "SensitivityAnalysis.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

# Some setup is needed for documentation deployment, see “Hosting Documentation” and
# deploydocs() in the Documenter manual for more information.
deploydocs(
    repo = "github.com/tpapp/SensitivityAnalysis.jl.git",
    push_preview = true
)
