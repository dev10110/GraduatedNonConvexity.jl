using GraduatedNonConvexity
using Documenter

DocMeta.setdocmeta!(GraduatedNonConvexity, :DocTestSetup, :(using GraduatedNonConvexity); recursive=true)

makedocs(;
    modules=[GraduatedNonConvexity],
    authors="Devansh Ramgopal Agrawal <devansh@umich.edu> and contributors",
    repo="https://github.com/dev10110/GraduatedNonConvexity.jl/blob/{commit}{path}#{line}",
    sitename="GraduatedNonConvexity.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dev10110.github.io/GraduatedNonConvexity.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/dev10110/GraduatedNonConvexity.jl",
    devbranch="main",
)
