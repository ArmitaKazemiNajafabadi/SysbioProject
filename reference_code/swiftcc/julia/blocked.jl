using LinearAlgebra, SparseArrays, JuMP, GLPK

function blocked(S, rev)
    m,n = size(S)
    irev = ?

    model = Model(GLPK.Optimizer)
    A = [S' -sparse(1:n, 1:n, 1)]
    lb = [-fill(Inf, m); fill(0.0, n)]
    ub = [fill(Inf, m); fill(0.0, n)]
    
end