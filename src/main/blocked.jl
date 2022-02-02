module new
export remove_blocked

using 
    LinearAlgebra, 
    SparseArrays, 
    JuMP, 
    GLPK       

function find_blocked_irrev(S, rev)
#=
    
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
    
    return result .≈ 1
end

function remove_blocked_irrev(S, rev)
#=
    removes all blocked irreversible reactions
    and returns new S, rev
=#
    unblocked = find_blocked_irrev(S, rev)
    return S[:, unblocked], rev[:, unblocked]
end 

function find_blocked_rev(S, rev)
#=
    S: dense stoich matrix (n*m)
=#
    K = nullspace(S)

    n, b = size(K)  #   n: number of reactions
                    #   b: basis cardinal
    
    is_blocked = falses(n)
    all_zero = false
    for i in 1:n
        all_zero = true
        for j in 1:b 
            if K[i,j] != 0
                all_zero = false
                break
            end
        end
        if all_zero
            is_blocked[i] = true
        end
    end

    return blocked
end

function remove_blocked_rev(S, rev)
    unblocked = .! find_blocked_rev(S, rev)
    return S[:, unblocked], rev[:, unblocked]
end

function remove_blocked(S, rev)
    S, rev = remove_blocked_irrev(S, rev)
    S, rev = remove_blocked_rev(S, rev)
    return S, rev
end

end