include("MRSFinder.jl")
include("ImportData.jl")
include("DataStructures")
include("LightGraphs.jl")
using .MRSFinder
using .ImportData
using .DataStructures
using .LightGraphs

#=
    group_identifier: is an array, representing the number of each group members. Each group is assosiated with a binary variable.
    
=#

S, rev = getParsedData() #obtain S, rev, ... from metabolic network
blocked_reactions, S1, rev1 = remove_blocked_reactions(S, rev)
essential_reactions, S2, rev2 = identify_essential_reactions(S1, rev)
m2, n2 = Size(S2)
setOfReactions = IntDisjointSet(n2)
group_identifier = group_reactions(setOfReactions, S2, rev2)
all_IRGs = findAllIRGs(S,V_biomass)
MRSs = findAllMRSs(blocked_reactions,essential_reactions,all_IRGs,group_identifier)







function group_reactions(setOfReactions, S, rev)
    # identify incoming/outgoing boundary reactions
    incoming_reactions = []
    m, n = Size(S)
    rows = rowvals(S)
    vals = nonzeros(S)
    eachMetaboliteProducers, eachMetaboliteConsumers = graphSearchRequirements(m,n,rows,vals,S,rev)
    startSearchingNeighbourhood(eachMetaboliteProducers, eachMetaboliteConsumers, setOfReactions, S, rows, vals, rev, depth)
   
    flag = false
    break_flag = false
    is_positive = true
    
    for i in 1:n
        flag = false
        break_flag = false

        for j in 1:m
            if !flag
                if S[i,j] != 0
                    flag = true
                    is_positive = S[i,j] > 0
                    continue
                end
            else
                if (S[i,j] > 0) != is_positive
                    break_flag = true
                    break
                end
            end
            if !break_flag
                if is_positive == true || rev[i] == true    # either incoming or reversible exiting
                    push!(incoming_reactions, i)
                end
            end
        end
    end

    
    


    m, n = Size(S)
    edgeAttributes = zeros(3,n)
    FirstListOfCurrentEdges = extract_boundary_reactions(S)
    while true
        for edg in listOfFirstSeenEdges
        CurrentList = []
            for neigh in edg
                if were_equivalent(edg,neigh) && 
                    push!(CurrentList, neigh)
        end
        
    end
    
end

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

function checkCoupleness(j,whichRow)
end
