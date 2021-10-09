module new
export calc_fctable

using LinearAlgebra

function calc_fctable(S, rev)
    
end

function is_fully_coupled(S, rev, i, j)
    #= 
    calculates whether reaction i is fully coupled to j
    returns c where v_i = c * v_j or 0 if not coupled
    =#
    
    m, n = size(S)
    
    jth, Sj = S[:, j], S[:, 1:end.!=j]

    # using nullspace
    ith, Sji = Sj[:, i], Sj[:, 1:end.!=i]
    K = nullspace(Sji')
    row, col = size(K)

    for j in 1:col
        inp = ith' * K[:, j]
        if isapprox(inp, 0, atol=1e-2)  # tod: change atol
            x = K[:,j] / inp
            c = jth' * x
            return c
        end
    end

    return 0
end

end