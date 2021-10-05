# importing packages
using Pkg
Pkg.add("SBML")
Pkg.add("JuMP")
Pkg.add("GLPK")
Pkg.add("DataStructures")
Pkg.add("LightGraphs")

using SBML
# model = readSBML("resources/e-coli-core.xml")

using JuMP, GLPK, LightGraphs

sparse()
