# Script to determine the range of safe controls to keep a state in some safe set.
# We calculate the minimal control needed so that the lower bound is invariant,
# and the maximal control making the upper bound invariant.
# These 2 range of controls might not necessarily intersect,
# yielding no feasible constant control making the safe set invariant.

using JLD, MosekTools, SCS, CSDP, SeDuMi, SDPA

# list_solvers = ["Mosek", "SCS", "CSDP", "SeDuMi", "SDPA"]
solver_name = "Mosek"
# solver_name = "SeDuMi"

# Safe_controls_ub = zeros(4,2)
# Safe_controls_lb = zeros(4,2)
U_min_sys1_v = zeros(1,10)
U_max_sys1_v = zeros(1,10)
i = 1

for np_ratio = 0.1:0.1:1
nq_ratio = 1


# choose an inverter node/bus for simulation
for sys_idx = 2   # [1], [2], [3] or [4]
    # choose the variable to study
    for var_idx = 1  # [1] frequency or [2] voltage

        flush(stdout)
        println("Using $solver_name.")
        println("Droop coefficients are $(100*np_ratio)% for the frequency and $(100*nq_ratio)% for the voltage.")
        eval(Meta.parse("solver = optimizer_with_attributes($solver_name.Optimizer, MOI.Silent() => true)"))
        model = SOSModel(solver)

        @time begin
        local out = sos_constraint(var_idx, sys_idx, model, np_ratio, nq_ratio)
        end # @time begin

        U_max_sys1_v[i] = out[1]
        U_min_sys1_v[i] = out[2]
        global i = i + 1

        # Safe_controls_ub[sys_idx, var_idx] = out[1]
        # Safe_controls_lb[sys_idx, var_idx] = out[2]
        # save("/Users/bouvier3/Documents/PNNL internship/Julia/invariance_data.jld", "Safe_controls_ub", Safe_controls_ub, "Safe_controls_lb", Safe_controls_lb)
        println("\n")
    end
end


end # np_ratio

# end#solvers

println(U_min_sys1_v)
println(U_max_sys1_v)
