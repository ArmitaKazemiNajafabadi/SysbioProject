S = [-1 0 0 0 0 0 0;
1 -1 1 0 0 0 0;
0 1 -1 0 0 0 0;
0 -1 0 1 0 0 0;
0 1 0 0 0 -1 0;
0 0 1 0 -1 0 0;
0 0 -1 0 0 0 1]
rev = [false false false false false false false]

m, n = size(S)
model = Model(GLPK.Optimizer)

# set variables and constraints based on LP (10) in the article
ub = [fill(10, n); fill(1.0, n)]
lb = [fill(-10, n); fill(0.0, n)]
@variable(model, lb[j] <= vu[j=1:n+n] <= ub[j])     # 0 <= u <= 1
for j in 1:n
    if rev[j]
        @constraint(model, vu[n+j] == 1.0)          # u[j] == 1 if rev[j], variable u just stands for irreversible reactions
    else
        @constraint(model, vu[n+j] <= vu[j])        # u <= v_I
    end
end
for j in 1:m
    @constraint(model, sum(S[j,k]*vu[k] for k in 1:n) == 0.0)   # Sv == 0
end

# set bojective function : sum u_i for i in R_I
@objective(model, Max, sum(vu[j] for j in [n + i for i in 1:n if !rev[i]]))

optimize!(model)

result = [value(vu[j]) for j in n+1:n+n]