using GeochemistryTools
using Test

include.(filter(contains(r".jl$"), readdir(@__DIR__; join=true)))
