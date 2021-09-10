module MRSFinder   
export findAllIRGs
export find_BI

using LinearAlgebra, SparseArrays, JuMP, GLPK

function remove_blocked_reactions(S, rev)
#=
    finds all blocked reactions
        finds irreversible blocked with one LP
        finds reversible blocked with many linear equations
=#
    irr_blocked = find_BI(S,rev)
    S2, rev2 = remove_reactions(S,rev,irr_blocked)
    x,n = size(S2, 2) # n is number of columns of S
    rev_blocked = fill(1, n)
    for i in 1:n
        if rev2[i]
            rev_blocked[i] = is_rev_blocked(S2,i)
        end
    end
    S3, rev3 = remove_reactions(S2,rev2,rev_blocked)
    return S3, rev3
end

#copied from https://github.com/mtefagh/swiftcore/blob/master/src/swiftcc.m
function remove_blocked_reactions(S, rev, varargin)            
    [m, n] = size(S);
    consistent = true(n, 1)
    # setting up the LP solver
    if ~isempty(varargin)
        solver = varargin{1};
    else
        solver = 'linprog';
    end
                  
    # identifying the blocked irreversible reactions
    result = remove_irrev_blocked_reactions(S, rev, solver);
    consistent(result.x(m+1:end) < -0.5) = false;
    
    # setting up the zero-tolerance parameter
    tol = norm(S(:, consistent), 'fro')*eps(class(S));
    
    # identifying the blocked reversible reactions
    [Q, R, ~] = qr(transpose(S(:, consistent)));
    Z = Q(rev(consistent) == 1, sum(abs(diag(R)) > tol)+1:end);
    
    # finding the consistent reactions of the original metabolic network
    consistent(consistent & rev == 1) = diag(Z*Z.') > tol^2;
    consistent = find(consistent);

end

#copied from https://github.com/mtefagh/swiftcore/blob/master/src/blocked.m
function  remove_irrev_blocked_reactions(S, rev, solver)
    [m, n] = size(S);
    irev = m + find(rev == 0);
    model.obj = zeros(m+n, 1);
    model.obj(irev) = 1;
    model.A = [S.', -speye(n)];
    model.sense = repmat('=', n, 1);
    model.sense(rev == 0) = '<';
    model.rhs = zeros(n, 1);
    model.lb = [-Inf(m, 1); zeros(n, 1)];
    model.lb(irev) = -1;
    model.ub = [Inf(m, 1); zeros(n, 1)];
    if strcmp(solver, 'gurobi') % gurobi
        params.outputflag = 0;
        result = gurobi(model, params);
        if ~strcmp(result.status, 'OPTIMAL')
            warning('Optimization is unstable!');
            fprintf('Optimization returned status: %s\n', result.status);
        end
    elseif strcmp(solver, 'linprog') % linprog
        problem.f = model.obj;
        problem.Aineq = model.A(rev == 0, :);
        problem.bineq = model.rhs(rev == 0);
        problem.Aeq = model.A(rev ~= 0, :);
        problem.beq = model.rhs(rev ~= 0);
        problem.lb = model.lb;
        problem.ub = model.ub;
        problem.solver = 'linprog';
        problem.options = optimset('Display', 'off');
        [result.x, result.objval, result.status, ~] = linprog(problem);
        if result.status ~= 1
            warning('Optimization is unstable!');
            fprintf('Optimization returned status: %d\n', result.status);
        end
    elseif strcmp(solver, 'cplex') % cplex
        problem.f = model.obj;
        problem.Aineq = model.A(rev == 0, :);
        problem.bineq = model.rhs(rev == 0);
        problem.Aeq = model.A(rev ~= 0, :);
        problem.beq = model.rhs(rev ~= 0);
        problem.lb = model.lb;
        problem.ub = model.ub;
        [result.x, result.objval, result.status] = cplexlp(problem);
        if result.status ~= 1
            warning('Optimization is unstable!');
            fprintf('Optimization returned status: %d\n', result.status);
        end
    else % COBRA
        model.b = model.rhs;
        model.c = model.obj;
        model.osense = 1;
        model.sense(model.sense == '=') = 'E';
        model.sense(model.sense == '<') = 'L';
        model.csense = model.sense;
        solution = solveCobraLP(model, 'solver', solver);
        result.x = solution.full;
        result.objval = solution.obj;
        result.status = solution.stat;
        if result.status ~= 1
            warning('Optimization is unstable!');
            fprintf('Optimization returned status: %d\n', result.status);
        end
    end
    
end

function is_rev_blocked(S, i)
#=
    returns 1 whether system of equations Sx = 0 and e^i x = 1 has a solution. returns 0 otherwise.
=#
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

function remove_reactions(S, rev, reaction_set)
#=
    remove i-th column of S and i-th component of rev for each i in reactions_set
=#
end

function generate_substitutable_solutions(new_ans)
    candidate_answers = []
    possible_answers = []

    push!(candidate_answers, new_ans)
    for #each existed reaction in new_ans
        if(#there were any substitutable reactions, consider that){
            push!(candidate_answers,candidate_ans)

        }
        end
    end
    for #each candidate answer:
        #solve an LP 
        if(#the candidate was appropriate){
            push!(possible_answers,candidate_ans)

        }
        end
    end
end

function findIRG(S,V_biomass,additional_conditions=0)
   # create MILP 
   # add some constraints based on fully coupled reactions (or remove some variables)
   # remove more variables with other ideas
   # solve MILP
    
    return IRG
end

function findAllIRGs(S,V_biomass)
    IRGs = []
    while true
        new_ans = findIRG(S,V_biomass,additional_conditions);         
        new_ans_size = size(new_ans)                                 #number of reastions
        
        if(new_ans_size > prev_ans_size){
            break                                                    #the new_ans is not an optimal solution 
        } else { 
            #=
                Here some substitutable solutions are generated and checked if they really hold the required additional conditions.
                If there were any answers, some additional conditions are generated to avoid finding previously-founded solutions.
            =#
            generated_solutions, additional_conditions = generate_substitutable_solutions(new_ans)
            push!(IRGs, new_ans)                                     #the new_ans is a new IRG (new solution)
        }
        end
    end
return IRGs
end

end