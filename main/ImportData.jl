# importing packages
using Pkg
Pkg.add("SBML")

using SBML
model = readSBML("resources/e-coli-core.xml")

