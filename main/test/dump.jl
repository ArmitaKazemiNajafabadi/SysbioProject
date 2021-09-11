# importing packages
using Pkg
Pkg.add("SBML")
Pkg.add("JuMP")
Pkg.add("GLPK")

using SBML
# model = readSBML("resources/e-coli-core.xml")

using JuMP, GLPK

sparse()
