%%% Whole network of the microgrid time simulation
%%% with neighbor interactions and control choices.

clear variables
clc
close all

% directory where the model, and results are stored
resul_dir = sprintf('%s/results',pwd);

%%% LOAD/READ the necessary model files, and control barrier functions
load(sprintf('%s/SysModel.mat',resul_dir),'x','f','SysDesc')

%%%%%%%%%%%%%%%% SETTINGS %%%%%%%%%%%%
% A change in settings modifies the upper and lower invariance controls.
% They must be re-calculated for each setting using the JULIA codes as
% explained below.

% Lower and upper bounds on the voltage
options.lb_v = -0.4;
options.ub_v = 0.2; 
% Lower and upper bounds on the frequency
options.lb_w = -3;
options.ub_w = 3; 
% Lower and upper bounds on the phase angle
options.lb_theta = -pi/2;
options.ub_theta = pi/2;
% v_nei is in the range [v - coupling, v + coupling].
options.coupling_v = 0.1;
% w_nei is in the range [w - coupling, w + coupling].
options.coupling_w = 0.2;
% Base value of the droop coefficients
options.np_value = 0.2*2*pi/5; % = 0.2513 rad/s
options.nq_value = 0.1/5; % = 0.02 p.u.
% Ratio of the droop coefficients considered
options.np_ratio = 0.17;
options.nq_ratio = 0.09;


%%%%%%%%%% UPPER AND LOWER INVARIANT CONTROLS %%%%%%%%%%%%%%

% If the settings above are modified and the error message below is
% displayed, follow the steps:
% 1) Verify that both np_ratio and nq_ratio (written above) are smaller than the lambda_star value obtained
% with  smart_droop.jl (JULIA code) for each inverter with the updated settings.
% 2) Update the settings in sos_constraint.jl (JULIA code).
% 3) Update the settings in invariance.jl (JULIA code) and run it to obtain
% the values of the invariant controls below.
if options.ub_w == 3 && options.ub_theta == pi/6 && options.np_ratio == 0.3 && options.nq_ratio == 0.8 && options.coupling_v == 0.02 && options.coupling_w == 0.12
    options.safe_u_lb = [-2.12, -0.36, -3.85, -0.57, -0.068, 0.0021, -2.21, -0.19];
    options.safe_u_ub = [2.76, 0.0095, 4.19, 0.18, 1.43, 0.0095, 3.28, 0.21];
    
elseif options.ub_w == 3 && options.ub_theta == pi/2 && options.np_ratio == 0.17 && options.nq_ratio == 0.09 && options.coupling_v == 0.1 && options.coupling_w == 0.2
    options.safe_u_lb = [0.11, -0.61, -2.69, -0.68, 0.43, 0.28, 0.074, 0.21];
    options.safe_u_ub = [2.85, 0.37, 4.27, 0.37, 0.6, 0.32, 0.58, 0.38];

else
    error('No precalculated inputs for safe set for this setup.\nUse invariance.jl on JULIA to calculate the control bounds for your specific setting.')
end


%%%%%%%%%%% INITIAL STATES %%%%%%%%%%

%%%% Random state in safe set within coupling distance of each other
coupling_theta = min(options.coupling_w, options.ub_theta);
theta = [0, coupling_theta*(rand(1,3)-0.5)];
v = options.lb_v + options.coupling_v/2 + rand*(options.ub_v - options.lb_v - options.coupling_v) + options.coupling_v*(rand(1,4)-0.5);
% w = options.lb_w + options.coupling_w/2 + rand*(options.ub_w - options.lb_w - options.coupling_w) + options.coupling_w*(rand(1,4)-0.5);
% If frequencies are too large, the angles quickly exit their assumed set, and the assumptions
% for the invariance controls are not verified anymore.
w = options.coupling_w*(rand(1,4)-0.5);

%%%%% Random state in safe set
% theta = [0, options.lb_theta + rand(1,3)*(options.ub_theta - options.lb_theta)];
% v = options.lb_v + rand(1,4)*(options.ub_v - options.lb_v);
% w = options.lb_w + rand(1,4)*(options.ub_w - options.lb_w);

%%%%% Fixed initial states for repeatability
% theta = [0, 0.3, -0.3, 0.1];
% w = [-0.1, -0.5, 0.5, 0.1];
% v = [-0.3, -0.2, -0.1, 0.1];

options.initial_state = [w(1), v(1), theta(2), w(2), v(2), theta(3), w(3), v(3), theta(4), w(4), v(4)];



%%%%%%%%%%%%%% OTHER SIMULATION PARAMETERS %%%%%%%%%%%%%%

% select control option
options.cflag = 2; % [0] no control, [1] random control, [2] invariance control

% simulation time-steps in seconds
options.tFinal = 2;     % total simulation duration in [s]
options.tSteps = 0.001;  % simulation time-step in [s]
options.tContr = 0.5;  % control changes at this time-step [s]
 

%%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% T: an array of time-stamps at the specified tSteps interval
% X: an array of the time-stamped values of the state-variables
% uval: an array of the time-stamped values of the applied control inputs
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[T,X,uval] = runSim(SysDesc, resul_dir, options);

 
 
%%%%%%%%%%%%%%%%%% PLOT RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%

% list of the inverters, whose controls need to be plotted
plot_controls_inverters = [1,2,4]; % any sublist of [1,2,3,4]


fontsize = 14;
nb_nodes = length(SysDesc);

%%%%% phase angle plot
figure
hold on
grid on
plot([T(1), T(end)], [0, 0], 'LineWidth',1.5)
angle_var_id = 3;
for sys_idx = 2:nb_nodes
    nb_states = length(SysDesc(sys_idx).x);
    plot(T, X(:,angle_var_id),'LineWidth',1.5)
    angle_var_id = angle_var_id + nb_states;
end
plot(T, options.lb_theta*ones(size(T)),'r--','LineWidth',1.5)
plot(T, options.ub_theta*ones(size(T)),'r--','LineWidth',1.5)
xlabel('time [s]')
ylabel('phase angle [rad]')
set(gca,'FontSize',fontsize)
legend('inv_1', 'inv_2', 'inv_3', 'inv_4')


%%%%% frequency plot
figure
hold on
grid on
frequency_var_id = -1;
for sys_idx = 1:nb_nodes
    nb_states = length(SysDesc(sys_idx).x);
    frequency_var_id = frequency_var_id + nb_states;
    plot(T, X(:,frequency_var_id),'LineWidth',1.5)
end
plot(T, options.lb_w*ones(size(T)),'r--','LineWidth',1.5)
plot(T, options.ub_w*ones(size(T)),'r--','LineWidth',1.5)
xlabel('time [s]')
ylabel('frequency [Hz]')
set(gca,'FontSize',fontsize)
legend('inv_1', 'inv_2', 'inv_3', 'inv_4')


%%%% voltage plot
figure
hold on
grid on
voltage_var_id = 0;
for sys_idx = 1:nb_nodes
    nb_states = length(SysDesc(sys_idx).x);
    voltage_var_id = voltage_var_id + nb_states;
    plot(T, X(:,voltage_var_id),'LineWidth',1.5)
end
plot(T, options.lb_v*ones(size(T)),'r--','LineWidth',1.5)
plot(T, options.ub_v*ones(size(T)),'r--','LineWidth',1.5)
xlabel('time [s]')
ylabel('voltage [p.u.]')
set(gca,'FontSize',fontsize)
legend('inv_1', 'inv_2', 'inv_3', 'inv_4')

 

%%%% control plots for the inverters given in the list plot_controls_inverters
for sys_idx = plot_controls_inverters
    for u_id = 2*sys_idx-1:2*sys_idx
    
        figure
        hold on
        grid on
        plot(T, uval(:, u_id),'LineWidth',1.5)
        plot(T, options.safe_u_ub(u_id)*ones(size(T)),'r--','LineWidth',1.5)
        plot(T, options.safe_u_lb(u_id)*ones(size(T)),'r--','LineWidth',1.5)
        xlabel('time [s]')
        set(gca,'FontSize',fontsize)
        
        if rem(u_id, 2) == 0
            ylabel(sprintf('control input, u_%i^q [p.u.]', sys_idx))
        else
            ylabel(sprintf('control input, u_%i^p [p.u.]', sys_idx))
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [tPoint, X, uval] = runSim(SysDesc, resul_dir, options)
 
tFinal = options.tFinal;
tSteps = options.tSteps;
tContr = options.tContr;
 
if rem(tContr,tSteps)~=0
    error('control timestep has to be in multiple of the simulation timestep')
end
 
tPoint = (0:tSteps:tFinal)';
 
% initialize the time-series vectors
nb_nodes = length(SysDesc);
nb_variables = 0;
for node_id = 1:nb_nodes
    nb_variables = nb_variables + length(SysDesc(node_id).x);
end
X = NaN(length(tPoint), nb_variables);
uval = NaN(length(tPoint),2*nb_nodes);
 
% Writing differential equation file
fname_ode = genSimModel_np_nq(SysDesc, resul_dir, options);
 
oldDir = cd(resul_dir); % move to the location of the ODE-file

X(1,:) = genInitial(SysDesc, options); % initalize X

% compute control
u0 = zeros(1,2*nb_nodes);
uval(1,:) = updateControl(1, SysDesc, X, u0, options);

for iT = 2:length(tPoint)
    eval(sprintf('[tTemp,xTemp] = ode45(@(t,x) %s(t,x,uval(iT-1,:)),tPoint(iT-1:iT),X(iT-1,:));',fname_ode))
    X(iT,:) = xTemp(end,:);
    if isnan(sum(X(iT,:)))
        error('reduce simulation time-step, or decrease the control effort')
    end

    % compute/update control

    if rem(iT-1,tContr/tSteps) == 0   % update uval
        uval(iT,:) = updateControl(iT, SysDesc, X, uval(iT-1,:), options);
    else    % re-use from previous time-step
        uval(iT,:) = uval(iT-1,:);
    end
end

cd(oldDir)
 
end
 




% Generate function network_ode in results folder, differential equation of
% the state variables
function [fname_ode] = genSimModel_np_nq(SysDesc, resul_dir, options)

fname_ode = sprintf('network_ode_np_%i_percent_nq_%i_percent', round(100*options.np_ratio), round(100*options.nq_ratio));
path = strcat(resul_dir,'/',fname_ode,'.m');

%%% Obtaining the dynamics
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 np_ratio nq_ratio
f = tayDynamics('val', options.np_ratio, options.nq_ratio);

%%% Writing ODE file
fid = fopen(path,'w');
fprintf(fid,'function f = %s(t,x,u)\n', fname_ode); 
fprintf(fid,'\n');
nb_nodes = length(SysDesc);
i_var_net = 1;
for sys_idx = 1:nb_nodes
    for i_var_sys = 1:length(SysDesc(sys_idx).x)
        c = char(SysDesc(sys_idx).x(i_var_sys));
        fprintf(fid, strcat(c,' = x(',num2str(i_var_net),');\n'));
        i_var_net = i_var_net + 1;
    end
end

fprintf(fid,'\n\nf = zeros(size(x));\n\n');
sys_idx = 1; i_var_sys = 1;
for i_var_net = 1:length(f)
    fexpr = string( f(i_var_net) );
    
    if i_var_sys == length(SysDesc(sys_idx).x)-1
        fprintf(fid,strcat('f(',string(i_var_net),') = ',fexpr,' + u(', num2str(2*(sys_idx-1)+1),');\n'));
    elseif i_var_sys == length(SysDesc(sys_idx).x)
        fprintf(fid,strcat('f(',string(i_var_net),') = ',fexpr,' + u(', num2str(2*sys_idx),');\n'));
        sys_idx = sys_idx + 1; i_var_sys = 0;
    else
        fprintf(fid,strcat('f(',string(i_var_net),') = ',fexpr,';\n'));
    end
    i_var_sys = i_var_sys + 1;
end
end


 


% update control at each time step 
function uval = updateControl(iT, SysDesc, X, u0, options)
 
nb_nodes = length(SysDesc);
switch options.cflag
    case 0      % no control
        uval = zeros(size(u0));
    case 1      % any other control input
        uval = -1 + 2*rand(size(u0));
    case 2      % invariance control from Bony-Brezis
        uval = zeros(size(u0));
        var_id = 1;
        for sys_idx = 1:nb_nodes % iterates on nodes
            nb_var = length(SysDesc(sys_idx).x); % number of variables for this node
            uval(2*(sys_idx-1)+1:2*sys_idx) = invariance_control(sys_idx, X(iT, var_id:var_id+nb_var-1), options);
            var_id = var_id + nb_var;
        end
    otherwise
        error('Control option not available')
end
end

% Takes id of the system and current state
% compute the minimum control need for the safe set to be invariant, based
% on precalculated data
function uval = invariance_control(sys_idx, x, options)

% lower and upper bounds of the safe sets for voltage (v) and frequency (w)
lb_v = options.lb_v;
ub_v = options.ub_v;
lb_w = options.lb_w;
ub_w = options.ub_w;

% voltage and frequency of the node
if sys_idx == 1
    w = x(1);
    v = x(2);
else
    w = x(2);
    v = x(3);
end

% upper and lower bounds controls for v and w
u_w_lb = options.safe_u_lb(2*(sys_idx-1)+1);
u_w_ub = options.safe_u_ub(2*(sys_idx-1)+1);
u_v_lb = options.safe_u_lb(2*sys_idx);
u_v_ub = options.safe_u_ub(2*sys_idx);

%%%% Random control inside safety admissible range
% if u_w_lb < u_w_ub && u_v_lb < u_v_ub
%     uval(1) = u_w_lb + rand*(u_w_ub - u_w_lb);
%     uval(2) = u_v_lb + rand*(u_v_ub - u_v_lb);
% else
%     error('No safety admissible range of controls.')
% end
% Choosing a random control can lead to safety violations as the coupling
% assumption between states is not always verified, and that theta is also not guaranteed
% to stay inside its given set, since there is no direct control over theta.
% A smoother control signal will be better: use double smooth step instead.

%%%% double smooth step control 
ratio_v = (v - lb_v)/(ub_v - lb_v);
ratio_w = (w - lb_w)/(ub_w - lb_w);
y_w = double_smooth_step(ratio_w);
uval(1) = u_w_lb + (u_w_ub - u_w_lb)*y_w;
y_v = double_smooth_step(ratio_v);
uval(2) = u_v_lb + (u_v_ub - u_v_lb)*y_v;

end

% Double smooth step function taking x in [0, 1] and giving him a value y in [0, 1]
function y = double_smooth_step(x)

%%%%% data from doubleSmoothStep.m
start_step_1 = 0.1;
end_step_1 = 0.4;
start_step_2 = 0.6;
end_step_2 = 0.9;
p1 = 74.074074074074074074074074074074*x^3 - 55.555555555555555555555555555556*x^2 + 8.8888888888888888888888888888889*x + 0.59259259259259259259259259259259;
p2 = 74.074074074074074074074074074074*x^3 - 166.66666666666666666666666666667*x^2 + 120.0*x - 28.0;

if x < start_step_1
    y = 1;
elseif x >= start_step_1 && x < end_step_1
    y = p1;
elseif x >= end_step_1 && x < start_step_2
    y = 0;
elseif x >= start_step_2 && x < end_step_2
    y = p2;
elseif x >= end_step_2
    y = -1;
else
    error('Problem with step function')
end
y = (-y+1)/2; % rescale from [-1,1] to [0,1]
end

% Initalize the state X of the whole network
function X = genInitial(SysDesc, options)

if isfield(options, 'initial_state') % initial state specified by user
    X = options.initial_state;
    
else % random initial state in the safe set
    nb_nodes = length(SysDesc);
    X = [];
    for sys_idx = 1:nb_nodes
        if sys_idx > 1
            X(end+1) = options.lb_theta + rand*(options.ub_theta - options.lb_theta); % theta
        end
        X(end+1) = options.lb_w + rand*(options.ub_w - options.lb_w); % omega
        X(end+1) = options.lb_v + rand*(options.ub_v - options.lb_v); % v
    end
end

end

