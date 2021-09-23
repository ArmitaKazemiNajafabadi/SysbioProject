module new
export find_blocked_irrev

using 
    LinearAlgebra, 
    SparseArrays, 
    JuMP, 
    GLPK       

function find_blocked_irrev(S, rev)
#=
    removes blocked irreversible functions
=#
    m, n = size(S)
    model = Model(GLPK.Optimizer)

    # set variables and constraints based on LP (10) in the article
    ub = [fill(10, n); fill(1.0, n)]
    lb = [fill(-10, n); fill(0.0, n)]
    @variable(model, lb[j] <= vu[j=1:n+n] <= ub[j])     # 0 <= u <= 1

    @variable(model, -Inf <= v[j=1:n] <= Inf)
    @variable(model, 0.0 <= u[j=1:n] <= 1.0)    # 0 <= u <= 1

    for j in 1:n
        if rev[j]
            @constraint(model, u[j] == 1.0)          # u[j] == 1 if rev[j], variable u just stands for irreversible reactions
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
    return result
end

end