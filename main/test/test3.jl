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



# using SBML

# model = readSBML("resources/e-coli-core.xml")
# metabolites, reactions, S = getS(model)

# reactions.


# flux_bounds


using DataStructures
a = DataStructures::IntDisjointSets(10)  # creates a forest comprised of 10 singletons
union!(a, 3, 5)          # merges the sets that contain 3 and 5 into one and returns the root of the new set
root_union!(a, x, y)     # merges the sets that have root x and y into one and returns the root of the new set
find_root!(a, 3)          # finds the root element of the subset that contains 3
in_same_set(a, x, y)     # determines whether x and y are in the same set
elem = push!(a)          # adds a single element in a new set; returns the new element
                         # (this operation is often called MakeSet)
num_groups(a)