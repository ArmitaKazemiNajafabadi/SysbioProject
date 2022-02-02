# importing packages
using Pkg
Pkg.add("SBML")
Pkg.add("JuMP")
Pkg.add("GLPK")

using SBML
using LinearAlgebra, SparseArrays
# model = readSBML("resources/e-coli-core.xml")

using JuMP, GLPK

S = Matrix{Float64}(undef, 4, 6)

S[1] = [-1,0,0,0,0,1]
S[2] = [2,-1,0,-1,0,0]
S[3] = [0,1,-1,0,0,0]
S[4] = [0,0,1,1,-1,1]
rev = [true false false]
#find_BI(S,rev)

