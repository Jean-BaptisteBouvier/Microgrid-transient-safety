# Finding the maximum of active power without using SOS
using SumOfSquares, MosekTools, SCS, CSDP, SeDuMi, SDPA, JLD, Suppressor, MAT


sys_idx = 2 # 1, 2, 3, 4
var_idx = 1 # [1] frequency, [2] voltage


# @suppress begin # suppress some irritating warning
# global dictionary = load("/Users/bouvier3/Documents/PNNL internship/Julia/droop_data.jld")
# end
# lambda_saved = dictionary["Droop_ratios"][sys_idx, var_idx]

# Lower and upper bounds on the voltage
lb_v = -0.4
ub_v = 0.2
# Lower and upper bounds on the frequency
lb_w = -3
ub_w = 3
# Lower and upper bounds on the angle
lb_theta = convert(Float64, -pi/6)
ub_theta = convert(Float64, pi/6)
# v_nei is in the range [v - coupling, v + coupling].
coupling_v = 0.02

flush(stdout)
if var_idx == 1 # frequency case
    println("Safe set is w in [$lb_w  $ub_w] with a voltage coupling of $coupling_v at node $sys_idx")
elseif var_idx == 2 # voltage case
    println("Safe set is v in [$lb_v $ub_v] with a coupling of $coupling_v at node $sys_idx")
else
    error("Need to select either voltage or frequency.")
end



cd("/Users/bouvier3/Documents/PNNL internship/Matlab/results")
vars = matread("SysModel.mat")
SysDesc = vars["SysDesc"]
neighbors = convert(Vector{Int}, SysDesc["nei"][sys_idx][2:end])

solver = optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)
model = SOSModel(solver)


# x_id in [[1,12]] refers to the x variable of interest
x_id = 3*(sys_idx-1) + var_idx + 1

@polyvar x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 np nq
global (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, np, nq)
f = tayDynamics(np, nq, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12)
x_dot = f[x_id] # dynamics of the chosen variable
inv_tau = eval(Meta.parse("subs(x_dot, np=>0, nq=>0, x$x_id=>-1)"))
dP = eval(Meta.parse("x_dot + inv_tau*x$x_id"))
inv_tau = convert(Float64, inv_tau)
dP = dP/inv_tau
dP = subs(dP, np=>1, nq=>1)

if mod(x_id,3) == 0  # voltage study, frequency does not matter

    # always true initialization since x1 = 0
    global ineq_ub = @set x1 >= 0
    global ineq_lb = @set x1 >= 0
    global ineq_theta = @set x1 >= 0
    # inequality constraints for upper and lower bound
    if sys_idx > 1 # theta of the other subsystems is not constant at 0, but bounded
        local theta_id = x_id-2
        global ineq_theta = @set ineq_theta && eval(Meta.parse("x$theta_id")) >= lb_theta && ub_theta >= eval(Meta.parse("x$theta_id"))
    end
    for nei = neighbors
        if nei > 1
            local theta_nei = 3*(nei-1) + 1; # id of the angle of the neighbor
            global ineq_theta = @set ineq_theta && eval(Meta.parse("x$theta_nei")) >= lb_theta && ub_theta >= eval(Meta.parse("x$theta_nei"))
        end
        local v_nei = nei*3; # id of the voltage of the neighbor
        global ineq_ub = @set ineq_ub && eval(Meta.parse("x$v_nei")) >= max(ub_v-coupling_v, lb_v) && ub_v >= eval(Meta.parse("x$v_nei"))
        global ineq_lb = @set ineq_lb && eval(Meta.parse("x$v_nei")) >= lb_v && min(lb_v+coupling_v, ub_v) >= eval(Meta.parse("x$v_nei"))
    end

    ineq_ub = @set ineq_ub && ineq_theta
    ineq_lb = @set ineq_lb && ineq_theta

    # upper bound of v
    dP_ub = subs(dP, eval(Meta.parse("x$x_id"))=>ub_v) # evaluate v at its upper bound
    @variable(model, t)
    @constraint(model, dP_ub <= t, domain = ineq_ub, maxdegree = 4)
    @objective(model, Min, t)
    optimize!(model)
    P_ub = value(t)

    # lower bound of v
    dP_lb = subs(dP, eval(Meta.parse("x$x_id"))=>lb_v) # evaluate v at its lower bound
    @variable(model, s)
    @constraint(model, dP_lb >= s, domain = ineq_lb, maxdegree = 4)
    @objective(model, Max, s)
    optimize!(model)
    P_lb = value(s)

    flush(stdout)
    println("Q_min = $P_lb  and Q_max = $P_ub")

    lambda_star = (ub_v - lb_v)/(P_ub - P_lb)


else # frequency study with 4 cases, because both v and w are coupled
    # Initialization of the inequality constraints
    global ineq_theta = @set x1 >= 0

    v_id = x_id+1
    global ineq_v_ub = @set eval(Meta.parse("x$v_id")) >= (ub_v-coupling_v) && ub_v >= eval(Meta.parse("x$v_id"))
    global ineq_v_lb = @set eval(Meta.parse("x$v_id")) >= lb_v && (lb_v+coupling_v) >= eval(Meta.parse("x$v_id"))

    if sys_idx > 1
        local theta_id = x_id-1
        global ineq_theta = @set ineq_theta && eval(Meta.parse("x$theta_id")) >= lb_theta && ub_theta >= eval(Meta.parse("x$theta_id"))
    end

    for nei = neighbors
        if nei > 1
            local theta_nei = 3*(nei-1) + 1; # id of the angle of the neighbor
            global ineq_theta = @set ineq_theta && eval(Meta.parse("x$theta_nei")) >= lb_theta && ub_theta >= eval(Meta.parse("x$theta_nei"))
        end
        local v_nei = nei*3; # id of the voltage of the neighbor
        global ineq_v_ub = @set ineq_v_ub && eval(Meta.parse("x$v_nei")) >= (ub_v-coupling_v) && ub_v >= eval(Meta.parse("x$v_nei"))
        global ineq_v_lb = @set ineq_v_lb && eval(Meta.parse("x$v_nei")) >= lb_v && (lb_v+coupling_v) >= eval(Meta.parse("x$v_nei"))
    end

    # upper bound of w
    dP_ub = subs(dP, eval(Meta.parse("x$x_id"))=>ub_w) # evaluate w at its upper bound
    @variable(model, ub_ub)
    ub_ub_domain = @set ineq_v_ub && ineq_theta
    @constraint(model, dP_ub <= ub_ub, domain = ub_ub_domain, maxdegree = 4)
    @objective(model, Min, ub_ub)
    optimize!(model)
    P_ub = value(ub_ub)

    @variable(model, ub_lb)
    ub_lb_domain = @set ineq_v_lb && ineq_theta
    @constraint(model, dP_ub <= ub_lb, domain = ub_lb_domain, maxdegree = 4)
    @objective(model, Min, ub_lb)
    optimize!(model)
    P_ub = max(value(ub_lb), P_ub)

    # lower bound of w
    dP_lb = subs(dP, eval(Meta.parse("x$x_id"))=>lb_w) # evaluate w at its lower bound
    @variable(model, lb_ub)
    lb_ub_domain = @set ineq_v_ub && ineq_theta
    @constraint(model, dP_lb >= lb_ub, domain = lb_ub_domain, maxdegree = 4)
    @objective(model, Max, lb_ub)
    optimize!(model)
    P_lb = value(lb_ub)

    @variable(model, lb_lb)
    lb_lb_domain = @set ineq_v_lb && ineq_theta
    @constraint(model, dP_lb >= lb_lb, domain = lb_lb_domain, maxdegree = 4)
    @objective(model, Max, lb_lb)
    optimize!(model)
    P_lb = min(value(lb_lb), P_lb)

    flush(stdout)
    println("P_min = $P_lb  and P_max = $P_ub")

    lambda_star = (ub_w - lb_w)/(P_ub - P_lb)
end

flush(stdout)
println("lambda_star = $lambda_star")
# println("lambda_star = $lambda_star while saved value is $lambda_saved")
