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
    m, n = size(S)
    model = Model(GLPK.Optimizer)

    # set variables and constraints based on LP (10) in the article
    @variable(model, -Inf <= v[j=1:n] <= Inf)
    @variable(model, 0.0 <= u[j=1:n] <= 1.0)        # 0 <= u <= 1

    for j in 1:n
        if rev[j]
            @constraint(model, u[j] == 1.0)         # u[j] == 1 if rev[j], variable u just stands for irreversible reactions
        else
            @constraint(model, u[j] <= v[j])        # u <= v_I
        end
    end
    for j in 1:m
        @constraint(model, sum(S[j,k]*v[k] for k in 1:n) == 0.0)   # Sv == 0
    end
    # set objective function : sum u_i for i in R_I
    @objective(model, Max, sum([u[j] for j in 1:n]))
    optimize!(model)
    result = [value(u[j]) for j in 1:n]

    ############################ IN CHERA INO RETURN MIKONE?
    return result .≈ 1
end

end