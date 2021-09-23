module ImportData
export getParsedData

using SBML

function getParsedData(path)
    #=
        returns S, rev
    =#

    model = readSBML(path)

    # metabolites, reactions, S = getS(model)
    metabolites, reactions, S = stoichiometry_matrix(model)

    reversibles = [v.reversible for (x,v) in model.reactions]
    return metabolites, reactions, S, reversible
end

end
