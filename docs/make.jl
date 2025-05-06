cd(@__DIR__)

using Documenter
using DocumenterCitations
using GeochemistryTools

pages = [
    "Introduction" => "index.md"
]

# bib = CitationBibliography(joinpath(@__DIR__, "src", "GeochemistryTools.jl.bib");
# style = :authoryear)
makedocs(
    sitename = "GeochemistryTools.jl",
    authors = "Jarred Lloyd",
    format = Documenter.HTML(),
    plugins = [bib],
    checkdocs= :export,
    pages = [
        "Home" => "index.md",
        "Geochronology" => [
            "Isochrons" => "Geochronology/isochrons.md",
        ],
        "Contributing" => "contributing.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
