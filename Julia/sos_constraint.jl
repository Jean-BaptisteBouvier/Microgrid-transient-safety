using SumOfSquares, DynamicPolynomials, MAT

# var_idx is the id of the variable of sys_idx to be kept in safe set,
# var_idx is [1]frequency or [2]voltage
# Determine the upper bound and lower of the controls needed to keep state in safe set.
function sos_constraint(var_idx, sys_idx, model, np_ratio, nq_ratio)



# Lower and upper bounds on the voltage
lb_v = -0.4
ub_v = 0.2
# Lower and upper bounds on the frequency
lb_w = -3
ub_w = 3
# Lower and upper bounds on the angle
lb_theta = -pi/6
ub_theta = pi/6
# v_nei is in the range [v - coupling, v + coupling].
coupling_v = 0.02
# w_nei is in the range [w - coupling, w + coupling].
coupling_w = 0.12

flush(stdout)
if var_idx == 1 # frequency case
    println("Safe set is w in [$lb_w  $ub_w] with a coupling of $coupling_w for node $sys_idx")
elseif var_idx == 2 # voltage case
    println("Safe set is v in [$lb_v $ub_v] with a coupling of $coupling_v for node $sys_idx")
else
    error("Need to select either voltage or frequency.")
end


cd("/Users/bouvier3/Documents/PNNL internship/Matlab/results")
vars = matread("SysModel.mat")
SysDesc = vars["SysDesc"]
neighbors = convert(Vector{Int}, SysDesc["nei"][sys_idx][2:end])

# x_id in [[1,12]] refers to the x variable of interest
x_id = 3*(sys_idx-1) + var_idx + 1

@polyvar x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12
global (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12)
f = tayDynamics(np_ratio, nq_ratio, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12)
x_dot = f[x_id] # dynamics of the chosen variable


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
        global ineq_ub = @set ineq_ub && eval(Meta.parse("x$v_nei")) >= (ub_v-coupling_v) && ub_v >= eval(Meta.parse("x$v_nei"))
        global ineq_lb = @set ineq_lb && eval(Meta.parse("x$v_nei")) >= lb_v && (lb_v+coupling_v) >= eval(Meta.parse("x$v_nei"))
    end

    ineq_ub = @set ineq_ub && ineq_theta
    ineq_lb = @set ineq_lb && ineq_theta

    # upper bound of v
    v_dot_ub = subs(x_dot, eval(Meta.parse("x$x_id"))=>ub_v) # evaluate v at its upper bound
    @variable(model, t)
    @constraint(model, v_dot_ub <= t, domain = ineq_ub, maxdegree = 4)
    @objective(model, Min, t)
    optimize!(model)
    max_u_ub = -value(t)

    # lower bound of v
    v_dot_lb = subs(x_dot, eval(Meta.parse("x$x_id"))=>lb_v) # evaluate v at its lower bound
    @variable(model, s)
    @constraint(model, v_dot_lb >= s, domain = ineq_lb, maxdegree = 4)
    @objective(model, Max, s)
    optimize!(model)
    min_u_lb = -value(s)

else # frequency study with 4 cases, because both v and w are coupled
    # Initialization of the inequality constraints
    global ineq_w_ub = @set x1 >= 0
    global ineq_w_lb = @set x1 >= 0
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

        local w_nei = 3*(nei-1) + 2; # id of the frequency of the neighbor
        global ineq_w_ub = @set ineq_w_ub && eval(Meta.parse("x$w_nei")) >= (ub_w-coupling_w) && ub_w >= eval(Meta.parse("x$w_nei"))
        global ineq_w_lb = @set ineq_w_lb && eval(Meta.parse("x$w_nei")) >= lb_w && (lb_w+coupling_w) >= eval(Meta.parse("x$w_nei"))
    end

    # upper bound of w
    w_dot_ub = subs(x_dot, eval(Meta.parse("x$x_id"))=>ub_w) # evaluate w at its upper bound
    @variable(model, ub_ub)
    ub_ub_domain = @set ineq_w_ub && ineq_v_ub && ineq_theta
    @constraint(model, w_dot_ub <= ub_ub, domain = ub_ub_domain, maxdegree = 4)
    @objective(model, Min, ub_ub)
    optimize!(model)
    max_u_ub = -value(ub_ub)

    @variable(model, ub_lb)
    ub_lb_domain = @set ineq_w_ub && ineq_v_lb && ineq_theta
    @constraint(model, w_dot_ub <= ub_lb, domain = ub_lb_domain, maxdegree = 4)
    @objective(model, Min, ub_lb)
    optimize!(model)
    max_u_ub = min(-value(ub_lb), max_u_ub)

    # lower bound of w
    w_dot_lb = subs(x_dot, eval(Meta.parse("x$x_id"))=>lb_w) # evaluate w at its lower bound
    @variable(model, lb_ub)
    lb_ub_domain = @set ineq_w_lb && ineq_v_ub && ineq_theta
    @constraint(model, w_dot_lb >= lb_ub, domain = lb_ub_domain, maxdegree = 4)
    @objective(model, Max, lb_ub)
    optimize!(model)
    min_u_lb = -value(lb_ub)

    @variable(model, lb_lb)
    lb_lb_domain = @set ineq_w_lb && ineq_v_lb && ineq_theta
    @constraint(model, w_dot_lb >= lb_lb, domain = lb_lb_domain, maxdegree = 4)
    @objective(model, Max, lb_lb)
    optimize!(model)
    min_u_lb = max(-value(lb_lb), min_u_lb)
end

flush(stdout)
println("max_u_ub = $max_u_ub    min_u_lb = $min_u_lb")
return [max_u_ub, min_u_lb]

end #function sos_constraint
