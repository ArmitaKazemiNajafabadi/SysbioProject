module MRSFinder   
export findAllMRS
export find_BI
export is_rev_blocked

using LinearAlgebra, SparseArrays, JuMP, GLPK

function remove_blocked_reactions!(S, rev)
#=
    finds all blocked reactions
    finds irreversible blocked with one LP
    finds reversible blocked with many linear equations
=#
    irr_blocked = find_BI(S, rev)
    @remove_reactions!(S, rev, irr_blocked)
    n = size(S, 2) # n is number of columns of S
    
    rev_blocked = fill(1, n)
    blocked = irr_blocked .≈ -1

    for i in 1:n
        if rev[i]
            rev_blocked[i] = is_rev_blocked_LP(S, i, blocked)
        end
    end
    @remove_reactions!(S, rev, rev_blocked)
    return S, rev
end

function is_rev_blocked_EQ(S, i)
#=
    returns 1 whether system of equations Sx = 0 and e_i x = 1 has a solution. returns 0 otherwise.
=#

    # [Q, R, ~] = qr(transpose(S));
    # Z = Q(rev(blocked == 0) == 1, sum(abs(diag(R)) > tol)+1:end);
    # blocked(blocked == 0 & rev == 1) = diag(Z*Z.') < tol^2;
    # # removing the blocked reactions from the metabolic network
    # S(:, blocked == 1) = [];
    # rev(blocked == 1) = [];
    # rxns(blocked == 1) = [];

    blocked = result .≈ -1
    finalBlocked = copy(blocked)
    Z = nullspace(Matrix(S[:, .!blocked]))
    blocked = [norm(Z[j, :]) < norm(S, 2)*eps(Float64) for j in 1:size(Z, 1)]
    
    finalBlocked[.!finalBlocked] = blocked

end



function find_BI(S,rev)
    m, n = size(S)
    model = Model(GLPK.Optimizer)

    # set variables and constraints based on LP (10) in the article
    ub = [fill(10, n); fill(1.0, n)]
    lb = [fill(-10, n); fill(0.0, n)]
    @variable(model, lb[j] <= vu[j=1:n+n] <= ub[j])     # 0 <= u <= 1
    for j in 1:n
        if rev[j]
            @constraint(model, vu[n+j] == 1.0)          # u[j] == 1 if rev[j], variable u just stands for irreversible reactions
        else
            @constraint(model, vu[n+j] <= vu[j])        # u <= v_I
        end
    end
    for j in 1:m
        @constraint(model, sum(S[j,k]*vu[k] for k in 1:n) == 0.0)   # Sv == 0
    end

    # set bojective function : sum u_i for i in R_I
    @objective(model, Max, sum(vu[j] for j in [n + i for i in 1:n if !rev[i]]))

    optimize!(model)

    result = [value(vu[j]) for j in n+1:n+n]
    return result
end

function is_fully_coupled(S, i , j)
#=
    returns c, 1 when system (S^(j))^T x = (e_i)^(j) has solution. S^(j) denotes the matrix S without j-th row
    otherwise returns 0
=#
    # S2 = remove j-th row of S
    # e2 = remove j-th row of e_i
    # determine whether system S2 x = e2 has a solution

end

function reduce_variables_based_on_fully_coupled(model)
#=
    adds constraint v_i = c v_j for any i,j that Ri and Rj are fully coupled
    adds constraint z_i = z_j too (z is integer varialbes of MILP)
=#

end

function any_other_idea_for_reducing_variables()

end

function essential_identifier(S, rev)
    println(":|")
end

macro remove_reactions!(S, rev, reaction_set)
#=
    remove i-th column of S and i-th component of rev for each i in reactions_set
=#
    S = S[:, .!reaction_set]
    rev = rev[.!reaction_set]
end
function is_rev_blocked_LP(S, i, blocked)
#=
    returns 1 whether system of equations Sx = 0 and e_i x = 1 has a solution. returns 0 otherwise.
=#

    b = zeros(1,n)
    m, n = size(S)
    model = Model(GLPK.Optimizer)

    # set variables and constraints based on LP (10) in the article
    @variable(model, x[j=1:n] )   
    for j in 1:m
        @constraint(model, sum(S[j,k]*x[k] for k in 1:n) == 0.0)   # Sx == 0
    end
    @constraint(model, x[i] <= 1.0)   # e_i x = 1

    # set bojective function : sum u_i for i in R_I
    @objective(model, Max, x[i])
    optimize!(model)
    
    return x[i]
    
end

function is_rev_blocked_System(S, i)
    #=
        returns 1 whether system of equations Sx = 0 and e_i x = 1 has a solution. returns 0 otherwise.
    =#
        m, n = size(S)
        e_i = zeros(n)
        e_i[i] = 1
        S2 = [S;n]
        b = zeros(1,m)
        b[m+1] = 1
        res = S2 \ b
        return res
end

function find_BI(S,rev)
    m, n = size(S)
    model = Model(GLPK.Optimizer)

    # set variables and constraints based on LP (10) in the article
    ub = [fill(10, n); fill(1.0, n)]
    lb = [fill(-10, n); fill(0.0, n)]
    @variable(model, lb[j] <= vu[j=1:n+n] <= ub[j])     # 0 <= u <= 1
    for j in 1:n
        if rev[j]
            @constraint(model, vu[n+j] == 1.0)          # u[j] == 1 if rev[j], variable u just stands for irreversible reactions
        else
            @constraint(model, vu[n+j] <= vu[j])        # u <= v_I
        end
    end
    for j in 1:m
        @constraint(model, sum(S[j,k]*vu[k] for k in 1:n) == 0.0)   # Sv == 0
    end

    # set bojective function : sum u_i for i in R_I
    @objective(model, Max, sum(vu[j] for j in [n + i for i in 1:n if !rev[i]]))

    optimize!(model)

    result = [value(vu[j]) for j in n+1:n+n]
    return result
end

function is_fully_coupled(S, i , j)
#=
    returns c, 1 when system (S^(j))^T x = (e_i)^(j) has solution. S^(j) denotes the matrix S without j-th row
    otherwise returns 0
=#
    # S2 = remove j-th row of S
    # e2 = remove j-th row of e_i
    # determine whether system S2 x = e2 has a solution

end

function reduce_variables_based_on_fully_coupled(model)
#=
    adds constraint v_i = c v_j for any i,j that Ri and Rj are fully coupled
    adds constraint z_i = z_j too (z is integer varialbes of MILP)
=#

end

function any_other_idea_for_reducing_variables()

end

function essential_identifier(S, rev)
    println(":|")
end

macro remove_reactions!(S, rev, reaction_set)
#=
    remove i-th column of S and i-th component of rev for each i in reactions_set
=#
    S = S[:, .!reaction_set]
    rev = rev[.!reaction_set]
end

function findMRS()
   # create MILP 
   # add some constraints based on fully coupled reactions (or remove some variables)
   # remove more variables with other ideas
   # solve MILP
   
end

function findAllMRS()
    
end


#------------------------------------------------------------------------------------------------------------------------------
function is_rev_blocked_LP(S, i)
    #=
        returns 1 whether system of equations Sx = 0 and e_i x = 1 has a solution. returns 0 otherwise.
    =#
        b = zeros(1,n)
        m, n = size(S)
        model = Model(GLPK.Optimizer)
        @variable(model, x[j=1:n] )   
        for j in 1:m
            @constraint(model, sum(S[j,k]*x[k] for k in 1:n) == 0.0)   # Sx == 0
        end
        @constraint(model, x[i] <= 1.0)   # e_i x <= 1
        @objective(model, Max, x[i])
        optimize!(model)
        return x[i] ≈ 1 ? 1 : 0
    end

end