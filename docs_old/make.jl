using Documenter
using StellarChem, StellarConstants, StellarEOS, StellarEvolution, StellarOpacity

makedocs(
    sitename = "Jems.jl",
    format=Documenter.HTML(),
    modules = [StellarChem, StellarConstants, StellarEOS, StellarEvolution, StellarOpacity]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/orlox/medusa.jl.git",
    deploy_config = Documenter.GitHubActions()
)
