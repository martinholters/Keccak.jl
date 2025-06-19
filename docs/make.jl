using Documenter, Keccak

makedocs(
    sitename="Keccak.jl",
    modules=[Keccak],
    #remotes=nothing,
)

deploydocs(
    repo = "github.com/martinholters/Keccak.jl.git",
)
