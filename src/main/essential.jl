module essential
export remove_essential_reactions

using 
    LinearAlgebra, 
    SparseArrays, 
    JuMP, 
    GLPK       

function remove_essential_reactions(S, rev)
    #=
        findes, then removes all essential reactions
        and returns new S, rev
    =#
    inessential = find_essential_reactions(S, rev)
    return S[:, inessential], rev[:, inessential]
end

function find_essential_reactions(S, rev)
#=
     findes all blocked essential reactions
        and returns an array of indices corresponding to reactions which were inessential.
=#

#= Ye hamchin LP bayad hal she:

min v_t
    such that:
            Sigma_{j=1}^{M} S_{ij}v_{j} = 0         i= 1, 2, ..., N
            v_{j}^{min} <= v_j <= v_{j}^{max}       j =1, 2, ..., M
            v_{biomass} >= v_{biomass}^{max}

baraye t = 1, 2, ..., M, yedune az in LP ha hal mikonim.
age min v_t shod == 0 va (r_t irreversible bud??) yani v_t inessential e!
age r_t reversible bud?? momkene majbur shim ""min -v_t"" ro hal konim (ba hamin natije)  

=#
    return inessential
end

end