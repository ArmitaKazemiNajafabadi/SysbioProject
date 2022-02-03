module graph
export startSearchingNeighbourhood
export searchNeighbourhood
export graphSearchRequirements

using 
    LinearAlgebra, 
    SparseArrays, 
    JuMP, 
    GLPK     

function startSearchingNeighbourhood(
    eachMetaboliteProducers, eachMetaboliteConsumers, setOfReactions, 
    S, rows, vals, rev, depth
    )
    
    m, n = Size(S)
    coupled_reaction_rates = zeros(n,n)
    for i = 1:n
        searchNeighbourhood(eachMetaboliteProducers, eachMetaboliteConsumers, setOfReactions, S, rows, vals, rev, depth,i , coupled_reaction_rates)
    end
end

function searchNeighbourhood(
    eachMetaboliteProducers, eachMetaboliteConsumers, setOfReactions, 
    S, rows, vals, rev, depth-1, j, coupled_reaction_rates
    )
    if depth==0 
        return
    end
    for i  = nzrange(S, j)
        whichRow = rows[i]
        val = vals[i]
        for ract_rate in eachMetaboliteConsumers[whichRow]
            if(coupled_reaction_rates[j,ract_rate[1]]!=0)
                coupleness_coefficient = checkCoupleness(j,ract_rate[1])
                searchNeighbourhood(eachMetaboliteProducers, eachMetaboliteConsumers, setOfReactions, S, rows, vals, rev, depth-1,ract_rate[1], coupled_reaction_rates)
            end
        end
    end
end

function graphSearchRequirements(m,n,rows,vals, S,rev)
    eachMetaboliteProducers = Array{MutableLinkedList{Any}}(undef, 1,m)
    eachMetaboliteConsumers = Array{MutableLinkedList{Any}}(undef, 1,m)
    
    for i=1:m
        eachMetaboliteProducers[i] = MutableLinkedList{Any}()
        eachMetaboliteConsumers[i] = MutableLinkedList{Any}()
    end
    
    for j = 1:n
        for i in nzrange(S, j)
            whichRow = rows[i]
            val = vals[i]
            if(rev[j]) 
                push!(eachMetaboliteProducers[whichRow],(j, abs(val)))
                push!(eachMetaboliteConsumers[whichRow],(j, -abs(val)))
            elseif (val>0)
                push!(eachMetaboliteProducers[whichRow],(j, val))
            else
                push!(eachMetaboliteConsumers[whichRow],(j, val))
            end     
        end
    end
    return eachMetaboliteProducers, eachMetaboliteConsumers
end




end