using Documenter
using GeochemistryTools

makedocs(
    sitename = "GeochemistryTools.jl",
    format = Documenter.HTML(),
    modules = [GeochemistryTools]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
