using Jems
using Documenter

DocMeta.setdocmeta!(Jems, :DocTestSetup, :(using Jems); recursive=true)

makedocs(;
    modules=[Jems],
    authors="Pablo Marchant <pablo.marchant@kuleuven.be>, Matthias Fabry <matthias.fabry@kuleuven.be>",
    repo="https://github.com/Pablo Marchant/Jems.jl/blob/{commit}{path}#{line}",
    sitename="Jems.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Pablo Marchant.github.io/Jems.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Pablo Marchant/Jems.jl",
    devbranch="main",
)
