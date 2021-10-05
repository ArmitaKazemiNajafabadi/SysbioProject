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
    for j = 1:n
        for i in nzrange(A, j)
            row = rows[i]
            val = vals[i]
        end
    end
        



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