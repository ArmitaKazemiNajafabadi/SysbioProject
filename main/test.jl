include("MRSFinder.jl")


using .MRSFinder
#=
println("hello world2")
blocked_identifier(1,2)
essential_identifier(1,2)

=#

#= 
    he he he
=#


S = [-1 0 0 0 0 0 0;
     1 -1 1 0 0 0 0;
     0 1 -1 0 0 0 0;
     0 -1 0 1 0 0 0;
     0 1 0 0 0 -1 0;
     0 0 1 0 -1 0 0;
     0 0 -1 0 0 0 1]
rev = [false false false false false false false]
result = find_BI(S,rev)
println(result)
