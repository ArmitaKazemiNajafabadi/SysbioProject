using 
    SBML, 
    LinearAlgebra, 
    SparseArrays, 
    JuMP, 
    GLPK

model = Model(GLPK.Optimizer)

@variable(model, x<=10)
@variable(model, u<=10)
# @constraint(model, x<=10)
# @constraint(model, u<=10)
@objective(model, Max, x+u)
optimize!(model)

result = value(x), value(u)
println(result)