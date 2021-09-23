include("MRSFinder.jl")
include("ImportData.jl")
include("DataStructures")
using .MRSFinder
using .ImportData
using .DataStructures

#=
    group_identifier: is an array, representing the number of each group members. Each group is assosiated with a binary variable.
    
=#

S, rev = getParsedData() #obtain S, rev, ... from metabolic network
blocked_reactions, S1 = remove_blocked_reactions(S, rev)
essential_reactions, S2 = identify_essential_reactions(S1, rev)
m2, n2 = Size(S2)
setOfReactions = IntDisjointSet(n2)
group_identifier = group_reactions(setOfReactions, S2)
all_IRGs = findAllIRGs(S,V_biomass)
MRSs = findAllMRSs(blocked_reactions,essential_reactions,all_IRGs,group_identifier)



function group_reactions(setOfReactions, S2)
    m, n = Size(S2)
    edgeAttributes = zeros(3,n)
    FirstListOfCurrentEdges = extract_boundary_reactions(S2)
    while true
        for edg in listOfFirstSeenEdges
        
        end

    end
    
end