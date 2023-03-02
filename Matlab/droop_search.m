%%%% Script to calculate the maximal droop coefficients at which there
%%%% exists a range of feasible constant controls to keep the state in the
%%%% safe set. Uses the dynamics from tayDynamics.m


clear variables; 
clc;
close all

% directory where the model, and results are stored
resul_dir = sprintf('%s/results',pwd);
load(sprintf('%s/SysModel.mat',resul_dir),'SysDesc')

% Lower and upper bounds on the voltage
options.lb_v = -0.4;
options.ub_v = 0.2; 

% Lower and upper bounds on the frequency
options.lb_w = -0.5;
options.ub_w = 0.5; 

% number of points to test for theta in [-pi/2, pi/2]
options.lb_theta = -pi/2;
options.ub_theta = pi/2;

% v_nei is in the range [v - coupling, v + coupling].
options.coupling_v = 0.05;
% w_nei is in the range [w - coupling, w + coupling].
options.coupling_w = 0.05;


% Droop_ratios = zeros(4,2);
% Safe_u = zeros(4,2);

Droop_ratios = [0.0361481000000000,0.448120000000000;0.0683898925781250,0.857177734375000;0.0288162231445313,0.100250244140625;0.0295791625976563,0.114715576171875];
Safe_u = [0.315068000000000,0.191356000000000;0.317248071497488,0.165475341252886;0.0875985993876611,0.316342129937311;0.0576412848785139,0.380674246649505];

% choose an inverter node/bus for simulation
for sys_idx = 2:4    % [1], [2], [3] or [4]
    % choose the variable to study
    for var_idx = 1:2 % [1]frequency or [2]voltage

        neighbors = SysDesc(sys_idx).nei(2:end);

        if var_idx == 1 % frequency case
            fprintf('Safe set is w in [%0.2f %0.2f] with a coupling of %0.2f for node %i\n',  options.lb_w, options.ub_w, options.coupling_w, sys_idx);   
        elseif var_idx == 2 % voltage case
            fprintf('Safe set is v in [%0.2f %0.2f] with a coupling of %0.2f for node %i\n',  options.lb_v, options.ub_v, options.coupling_v, sys_idx);
        else
            error('Need to select either [1]frequency or [2]voltage for var_idx.')
        end

        %%%%%%%%%%%%% COMPUTATIONS %%%%%%%%%%%%%%%

        tic();
        [u, droop_ratio] = sos_droop(neighbors, sys_idx, options, var_idx);
        toc()
        fprintf('For the droop ratio %0.4f the safe set is invariant with u = %0.4f\n', droop_ratio, u)
        Droop_ratios(sys_idx, var_idx) = droop_ratio;
        Safe_u(sys_idx, var_idx) = u;
%         save('Droop.mat', 'Droop_ratios', 'Safe_u')
        
    end
end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% var_idx is the id of the variable of sys_idx to be kept in safe set,
% var_idx is 1:frequency or 2:voltage
% Determine the single constant control u needed to keep state
% in the safe set at the maximal droop_ratio
function [u, droop_ratio] = sos_droop(neighbors, sys_idx, options, var_idx)

droop_ratio = 1;

% x_id in [[1,12]] refers to the x variable of interest
x_id = 3*(sys_idx-1) + var_idx + 1;


syms np_ratio nq_ratio
f = tayDynamics('sym');
x_dot = f(x_id-1);


%%% Voltage inputs
lb_v = options.lb_v; % lower bound on the voltage
ub_v = options.ub_v; % upper bound on the voltage
coupling_v = options.coupling_v; % coupling of voltages with neighbors

%%% Frequency inputs
lb_w = options.lb_w; % lower bound on the frequency
ub_w = options.ub_w; % upper bound on the frequency
coupling_w = options.coupling_w; % coupling of frequency with neighbors

%%% Angle inputs
lb_theta = options.lb_theta; % lower bound on the angle
ub_theta = options.ub_theta; % upper bound on the angle


syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12
sos_options.params.fid = 0; % no print

if mod(x_id,3) == 0  % voltage study, frequency does not matter
    vars = symvar(x_dot);
    if any([vars == x2, vars == x5, vars == x8, vars == x11])
        disp('Frequency matters, set to 0.\n')
        x_dot = subs(x_dot, [x2, x5, x8, x11], [0, 0, 0, 0]);
    end
    equality_constraints = [];
    
    %%% inequality constraints for upper and lower bound
    if sys_idx == 1
        [ineq_ub, ineq_lb] = deal([]); % theta_1 = 0
    else % theta of the other subsystems is not constant at 0, but bounded
        [ineq_ub, ineq_lb] = deal([ eval(sprintf('x%i-lb_theta', x_id-2)), eval(sprintf('ub_theta-x%i', x_id-2)) ]);
    end
    for nei = neighbors
        if nei == 1
            ineq_theta = []; % no theta_1
        else
            theta_nei = 3*(nei-1) + 1; % id of the angle of the neighbor
            ineq_theta = [ eval(sprintf('x%i-lb_theta', theta_nei)), eval(sprintf('ub_theta-x%i', theta_nei)) ];
        end
        v_nei = nei*3; % id of the voltage of the neighbor
        ineq_ub = [ineq_ub, ineq_theta, eval(sprintf('x%i - (ub_v-coupling_v)', v_nei)), eval(sprintf('ub_v - x%i', v_nei)) ];
        ineq_lb = [ineq_lb, ineq_theta, eval(sprintf('x%i - lb_v', v_nei)), eval(sprintf('(lb_v+coupling_v) - x%i', v_nei)) ];
    end
    
    v_dot_ub_droop = subs(x_dot, eval(sprintf('x%i', x_id)), ub_v); % evaluate v at its upper bound
    v_dot_lb_droop = subs(x_dot, eval(sprintf('x%i', x_id)), lb_v); % evaluate v at its lower bound
    
    % Determine an upper bound to the droop ratio
    while ~exist('max_ratio')
        droop_ratio = 2*droop_ratio;
        
        v_dot_ub = subs(v_dot_ub_droop, [np_ratio, nq_ratio], [droop_ratio, droop_ratio]);
        v_dot_lb = subs(v_dot_lb_droop, [np_ratio, nq_ratio], [droop_ratio, droop_ratio]);

        [max_u_ub, ~, ~] = findbound(-v_dot_ub, ineq_ub, equality_constraints, 4, sos_options);  
        [min_u_lb, ~, ~] = findbound(v_dot_lb, ineq_lb, equality_constraints, 4, sos_options);
        min_u_lb = -min_u_lb;
        if min_u_lb < max_u_ub
            min_ratio = droop_ratio;
        else
            max_ratio = droop_ratio;
        end
    end   
    % Determine a lower bound to the droop ratio
    while ~exist('min_ratio')
        droop_ratio = droop_ratio/2;
        
        v_dot_ub = subs(v_dot_ub_droop, [np_ratio, nq_ratio], [droop_ratio, droop_ratio]);
        v_dot_lb = subs(v_dot_lb_droop, [np_ratio, nq_ratio], [droop_ratio, droop_ratio]);

        [max_u_ub, ~, ~] = findbound(-v_dot_ub, ineq_ub, equality_constraints, 4, sos_options);  
        [min_u_lb, ~, ~] = findbound(v_dot_lb, ineq_lb, equality_constraints, 4, sos_options);
        min_u_lb = -min_u_lb;
        if min_u_lb < max_u_ub
            min_ratio = droop_ratio;
        else
            max_ratio = droop_ratio;
        end
    end   
    % Now that we have both bound, proceed to a bisection
    for iter = 1:10
        
        droop_ratio = (max_ratio + min_ratio)/2;
        
        v_dot_ub = subs(v_dot_ub_droop, [np_ratio, nq_ratio], [droop_ratio, droop_ratio]);
        v_dot_lb = subs(v_dot_lb_droop, [np_ratio, nq_ratio], [droop_ratio, droop_ratio]);

        [max_u_ub, ~, ~] = findbound(-v_dot_ub, ineq_ub, equality_constraints, 4, sos_options);  
        [min_u_lb, ~, ~] = findbound(v_dot_lb, ineq_lb, equality_constraints, 4, sos_options);
        min_u_lb = -min_u_lb;
        if min_u_lb < max_u_ub
            min_ratio = droop_ratio;
        else
            max_ratio = droop_ratio;
        end
    end   
    u = (max_u_ub + min_u_lb)/2;
    droop_ratio = (max_ratio + min_ratio)/2;
    
else % frequency study
    % 4 cases here, because both v and w are coupled
    equality_constraints = [];
    
    %%% inequality constraints
    [ineq_theta, ineq_w_ub, ineq_w_lb, ineq_v_ub, ineq_v_lb] = deal([]);
    
    if sys_idx >= 2
        ineq_theta = [ eval(sprintf('x%i-lb_theta', x_id-1)), eval(sprintf('ub_theta-x%i', x_id-1)) ];
    end
    ineq_v_ub = [ineq_v_ub, eval(sprintf('x%i - (ub_v-coupling_v)', x_id+1)), eval(sprintf('ub_v - x%i', x_id+1)) ];
    ineq_v_lb = [ineq_v_lb, eval(sprintf('x%i - lb_v', x_id+1)), eval(sprintf('(lb_v+coupling_v) - x%i', x_id+1)) ];
    
    
    for nei = neighbors
        
        if nei >= 2
            theta_nei = 3*(nei-1) + 1; % id of the angle of the neighbor
            ineq_theta = [ineq_theta, eval(sprintf('x%i-lb_theta', theta_nei)), eval(sprintf('ub_theta-x%i', theta_nei)) ];
        end
        v_nei = nei*3; % id of the voltage of the neighbor
        ineq_v_ub = [ineq_v_ub, eval(sprintf('x%i - (ub_v-coupling_v)', v_nei)), eval(sprintf('ub_v - x%i', v_nei)) ];
        ineq_v_lb = [ineq_v_lb, eval(sprintf('x%i - lb_v', v_nei)), eval(sprintf('(lb_v+coupling_v) - x%i', v_nei)) ];
        
        w_nei = 3*(nei-1) + 2; % id of the frequency of the neighbor
        ineq_w_ub = [ineq_w_ub, eval(sprintf('x%i - (ub_w-coupling_w)', w_nei)), eval(sprintf('ub_w - x%i', w_nei)) ];
        ineq_w_lb = [ineq_w_lb, eval(sprintf('x%i - lb_w', w_nei)), eval(sprintf('(lb_w+coupling_w) - x%i', w_nei)) ];
        
    end
    
    w_dot_ub_droop = subs(x_dot, eval(sprintf('x%i', x_id)), ub_w); % evaluate w at its upper bound
    w_dot_lb_droop = subs(x_dot, eval(sprintf('x%i', x_id)), lb_w); % evaluate w at its lower bound
    
    % Determine an upper bound to the droop ratio
    while ~exist('max_ratio')
        droop_ratio = 2*droop_ratio;
        
        w_dot_ub = subs(w_dot_ub_droop, [np_ratio, nq_ratio], [droop_ratio, droop_ratio]);
        w_dot_lb = subs(w_dot_lb_droop, [np_ratio, nq_ratio], [droop_ratio, droop_ratio]);

        [max_u_ub_v_lb, ~, ~] = findbound(-w_dot_ub, [ineq_theta, ineq_w_ub, ineq_v_lb], equality_constraints, 4, sos_options);
        [max_u_ub_v_ub, ~, ~] = findbound(-w_dot_ub, [ineq_theta, ineq_w_ub, ineq_v_ub], equality_constraints, 4, sos_options);
        max_u_ub = min(max_u_ub_v_lb, max_u_ub_v_ub);

        [min_u_lb_v_lb, ~, ~] = findbound(w_dot_lb, [ineq_theta, ineq_w_lb, ineq_v_lb], equality_constraints, 4, sos_options);
        [min_u_lb_v_ub, ~, ~] = findbound(w_dot_lb, [ineq_theta, ineq_w_lb, ineq_v_ub], equality_constraints, 4, sos_options);
        min_u_lb = max(-min_u_lb_v_lb, -min_u_lb_v_ub);
        
        if min_u_lb < max_u_ub
            min_ratio = droop_ratio;
        else
            max_ratio = droop_ratio;
        end
    end   
    % Determine a lower bound to the droop ratio
    while ~exist('min_ratio')
        droop_ratio = droop_ratio/2;
        
        w_dot_ub = subs(w_dot_ub_droop, [np_ratio, nq_ratio], [droop_ratio, droop_ratio]);
        w_dot_lb = subs(w_dot_lb_droop, [np_ratio, nq_ratio], [droop_ratio, droop_ratio]);

        [max_u_ub_v_lb, ~, ~] = findbound(-w_dot_ub, [ineq_theta, ineq_w_ub, ineq_v_lb], equality_constraints, 4, sos_options);
        [max_u_ub_v_ub, ~, ~] = findbound(-w_dot_ub, [ineq_theta, ineq_w_ub, ineq_v_ub], equality_constraints, 4, sos_options);
        max_u_ub = min(max_u_ub_v_lb, max_u_ub_v_ub);

        [min_u_lb_v_lb, ~, ~] = findbound(w_dot_lb, [ineq_theta, ineq_w_lb, ineq_v_lb], equality_constraints, 4, sos_options);
        [min_u_lb_v_ub, ~, ~] = findbound(w_dot_lb, [ineq_theta, ineq_w_lb, ineq_v_ub], equality_constraints, 4, sos_options);
        min_u_lb = max(-min_u_lb_v_lb, -min_u_lb_v_ub);
        
        if min_u_lb < max_u_ub
            min_ratio = droop_ratio;
        else
            max_ratio = droop_ratio;
        end
    end   
    % Now that we have both bound, proceed to a bisection
    for iter = 1:10
        
        droop_ratio = (max_ratio + min_ratio)/2;
            w_dot_ub = subs(w_dot_ub_droop, [np_ratio, nq_ratio], [droop_ratio, droop_ratio]);
        w_dot_lb = subs(w_dot_lb_droop, [np_ratio, nq_ratio], [droop_ratio, droop_ratio]);

        [max_u_ub_v_lb, ~, ~] = findbound(-w_dot_ub, [ineq_theta, ineq_w_ub, ineq_v_lb], equality_constraints, 4, sos_options);
        [max_u_ub_v_ub, ~, ~] = findbound(-w_dot_ub, [ineq_theta, ineq_w_ub, ineq_v_ub], equality_constraints, 4, sos_options);
        max_u_ub = min(max_u_ub_v_lb, max_u_ub_v_ub);

        [min_u_lb_v_lb, ~, ~] = findbound(w_dot_lb, [ineq_theta, ineq_w_lb, ineq_v_lb], equality_constraints, 4, sos_options);
        [min_u_lb_v_ub, ~, ~] = findbound(w_dot_lb, [ineq_theta, ineq_w_lb, ineq_v_ub], equality_constraints, 4, sos_options);
        min_u_lb = max(-min_u_lb_v_lb, -min_u_lb_v_ub);
        
        if min_u_lb < max_u_ub
            min_ratio = droop_ratio;
        else
            max_ratio = droop_ratio;
        end
    end   
    u = (max_u_ub + min_u_lb)/2;
    droop_ratio = (max_ratio + min_ratio)/2;
end


end