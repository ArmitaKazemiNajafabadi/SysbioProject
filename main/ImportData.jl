module ImportData
export getParsedData

using SBML

function getParsedData(path=0)
    #=
        returns S, rev
    =#

    model = readSBML("Ec_core_flux1.xml")
    metabolites, reactions, S = getS(model)

end

end
