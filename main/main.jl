include("MRSFinder.jl")
include("ImportData.jl")
using .MRSFinder
using .ImportData


#=
obtain S, rev, ... from metabolic network
=#
S, rev = getParsedData()
blocked_reactions = remove_blocked_reactions(S,rev)
essential_reactions = identify_essential_reactions(S,rev)
group_identifier = group_reactions(S,rev)
all_IRGs = findAllIRGs(S,V_biomass)
MRSs = findAllMRSs(blocked_reactions,essential_reactions,all_IRGs,group_identifier)

