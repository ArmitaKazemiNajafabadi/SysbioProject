# function do_change(S)
#     S[1] = 2;
#     S = [3,2,3]
# end

# S = [1,2,3]
# do_change(S)
# println(S)

# a = [false, true, false, true]
# .!a

# a = fill(1, 3)
# b = fill(true, 3)
# c = [a' fill(2, 3)']
# d = rand(10)

# d[c]

# println("TEEEEEEEST")
# A = [-1 1; 2 1; -1 1; -1 1]
# b = [0, 3, 1, 2]
# println(A)
# A\b

# function f(x, y=2) 
#     print(x+y)
# end

# f(1,1)
using SBML

model = readSBML("resources/e-coli-core.xml")
metabolites, reactions, S = getS(model)
    
reactions.
# flux_bounds