using Jems
using Documenter
using Literate

# Parse examples using Literate
pkg_path = pkgdir(Jems)

function ignore_code_blocks(content)
    content = replace(content, "##\n" => "\n")  # remove code blocks 
    content = replace(content, "###" => "##")  # make level 3 headers level 2
end

Literate.markdown(pkg_path * "/examples/NuclearBurning.jl", pkg_path * "/docs/src/", preprocess=ignore_code_blocks)

DocMeta.setdocmeta!(Jems, :DocTestSetup, :(using Jems); recursive=true)

makedocs(modules=[Jems],
         authors="Pablo Marchant <pablo.marchant@kuleuven.be>, Matthias Fabry <matthias.fabry@kuleuven.be>",
         repo="https://github.com/orlox/Jems.jl/blob/{commit}{path}#{line}",
         sitename="Jems.jl",
         format=Documenter.HTML(prettyurls=get(ENV, "CI", "false") == "true",
                                canonical="https://orlox.github.io/Jems.jl",
                                edit_link="main",
                                assets=String[]),
         pages=["Home" => "index.md",
                "Examples" => ["NuclearBurning.md"],
                "Modules" => ["Chem.md", "Constants.md", "EOS.md", "Evolution.md", "Opacity.md"],
                "Developing" => ["style.md"]])

deploydocs(repo="github.com/orlox/Jems.jl", devbranch="main")