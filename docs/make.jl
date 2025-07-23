using Documenter, Keccak

makedocs(
    sitename="Keccak.jl",
    modules=[Keccak],
    pages=[
        "Home" => "index.md",
        "User Guide" => [
            "fips-202.md",
            "sponge.md",
            "sp800-185.md",
            "keccak.md",
        ],
        "reference.md",
        "internals.md",
    ],
)

deploydocs(
    repo = "github.com/martinholters/Keccak.jl.git",
)
