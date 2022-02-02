include("old_blocked.jl")
include("blocked.jl")

using .old
using .new

S = [-1 0 0 0 0 0 0;
1 -1 1 0 0 0 0;
0 1 -1 0 0 0 0;
0 -1 0 1 0 0 0;
0 1 0 0 0 -1 0;
0 0 1 0 -1 0 0;
0 0 -1 0 0 0 1]
rev = [false false false false false false false]

a = find_BI_old(S, rev)
b = find_blocked_irrev(S, rev)

a == b