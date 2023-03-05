%%% Microgrid time simulation with neighbor interactions and control
%%% choices


clear variables
clc
close all

% directory where the model, and results are stored
resul_dir = sprintf('%s/results',pwd);

%%% LOAD/READ the necessary model files, and control barrier functions
load(sprintf('%s/SysModel.mat',resul_dir),'x','f','SysDesc')

% Lower and upper bounds on the voltage
options.lb_v = -0.4;
options.ub_v = 0.2; 

% Lower and upper bounds on the frequency
options.lb_w = -3;
options.ub_w = 3; 

% Lower and upper bounds on the phase angle
options.lb_theta = -pi/6;
options.ub_theta = pi/6;

% v_nei is in the range [v - coupling, v + coupling].
options.coupling_v = 0.02;
% w_nei is in the range [w - coupling, w + coupling].
options.coupling_w = 0.12;

% droop coefficients
options.np_value = 0.2*2*pi/5; % = 0.2513 rad/s
options.nq_value = 0.1/5; % = 0.02 p.u.

%%% Stochastic time-varying inputs or constant worst case inputs
options.varying_inputs = false; options.worst_inputs = true;
% options.varying_inputs = true; options.worst_inputs = false;


sys_idx = 1;
options.np_ratio = 0.4;
options.nq_ratio = 1;


if options.lb_w == -0.1 && options.ub_w == 0.1 && options.np_ratio * options.nq_ratio == 1
    % Safe controls upper and lower bounds for the whole system for w in [-0.1 0.1]
    Safe_controls_ub = [-181.878947904333,-5.31442471242917;-86.1469143007416,-2.51587816051324;-319.431909354404,-4.04772784064699;-275.843945916742,-1.66274765365961];
    Safe_controls_lb = [336.392376461998,24.8743677933817;167.288937195204,11.2651646246086;388.318909179419,103.813042751048;308.883668094873,88.2148208972067];
    
elseif options.lb_w == -0.5 && options.ub_w == 0.5 && options.np_ratio * options.nq_ratio == 1
    % Safe controls upper and lower bounds for the whole system for w in [-0.5 0.5]
    Safe_controls_ub = [-181.024472783785,-5.31442471242917;-85.3770072478801,-2.51587816051324;-318.660837336443,-4.04772784064699;-275.112240413204,-1.66274765365961];
    Safe_controls_lb = [335.537901422940,24.8743677933817;166.519030130003,11.2651646246086;387.547836768441,103.813042751048;308.151962625794,88.2148208972067];
    
elseif sys_idx == 1 && options.lb_w == -0.5 && options.ub_w == 0.5 && options.np_ratio == 0.02 && options.nq_ratio == 0.2
    Safe_controls_ub = [0.319, 0.647];
    Safe_controls_lb = [-0.383, -0.298];
%     Safe_controls_ub = [0.647, 0.319];
%     Safe_controls_lb = [-0.298, -0.383];
    Q_min = -1.0948908219309998; Q_max = 0.24433658560683955;
    P_min = -17.948679254058195; P_max = 9.705731855571944;
    
elseif sys_idx == 1 && options.lb_w == -3 && options.ub_w == 3 && options.np_ratio == 0.4 && options.nq_ratio == 1 && options.coupling_v == 0.02 && options.coupling_w == 0.12
    Safe_controls_ub = [1.571, -0.0937];
    Safe_controls_lb = [-0.724, -0.25];

else
    error('No precalculated inputs for safe set for these bounds on the frequency.')
end
options.safe_u_ub = Safe_controls_ub(sys_idx,:);
options.safe_u_lb = Safe_controls_lb(sys_idx,:);



% select control option
options.cflag = 3; % [0] no control, [1] barrier control, [2] random control, [3] invariance control  
 
% simulation time-steps in seconds
options.tFinal = 10;     % total simulation duration in [s]
options.tSteps = 0.01;  % simulation time-step in [s]
options.tContr = 1;  % control changes at this time-step [s]
options.tDistu = 0.01;  % neighboring states change at this time-step [s]
 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NOW ready to run the time-domain simulations
%
% T: an array of time-stamps at the specified tSteps interval
% X: an array of the time-stamped values of the subsystem state-variables
% yval: an array of the time-stamped values of the neighbor states
% uval: an array of the time-stamped values of the applied control inputs
% PQval: an array of the time-stamped values of the active and reactive power
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[T,X,yval,uval,PQval] = runSim(SysDesc, sys_idx, resul_dir, options);

 
 
%%%%%%%%%%%%%%%%%% SAVE AND PLOT results %%%%%%%%%%%%%%%%%%%%%%%%%%
fontsize = 14;

if 1    % plot time-domain state trajectories   
    % voltage plot
    figure
    hold on
    grid on
    plot(T,X(:,end),'LineWidth',1.5)
    plot(T, options.lb_v*ones(size(T)),'r--','LineWidth',1.5)
    plot(T, options.ub_v*ones(size(T)),'r--','LineWidth',1.5)
    xlabel('time [s]')
    ylabel(sprintf('(shifted) voltage v_%i [p.u.]', sys_idx))
    set(gca,'FontSize',fontsize)
    
    % frequency plot
    figure
    hold on
    grid on
    plot(T, X(:,end-1),'LineWidth',1.5)
    plot(T, options.lb_w*ones(size(T)),'r--','LineWidth',1.5)
    plot(T, options.ub_w*ones(size(T)),'r--','LineWidth',1.5)
    xlabel('time [s]')
    ylabel(sprintf('frequency \x03C9_%i [Hz]', sys_idx))
    set(gca,'FontSize',fontsize)
    
    % reactive power control input
    figure
    hold on
    grid on
    stairs(T(1:end-1), uval(1:end-1, 2),'k-','LineWidth',1.5)
    plot(T, options.safe_u_lb(2)*ones(size(T)),'r--','LineWidth',1.5)
    plot(T, options.safe_u_ub(2)*ones(size(T)),'r--','LineWidth',1.5)
    xlabel('time [s]')
    ylabel(sprintf('control input, u_%i^q [p.u.]', sys_idx))
    set(gca,'FontSize',fontsize)
   
    % active power control input
    figure
    hold on
    grid on
    stairs(T(1:end-1), uval(1:end-1, 1),'k-','LineWidth',1.5)
    plot(T, options.safe_u_lb(1)*ones(size(T)),'r--','LineWidth',1.5)
    plot(T, options.safe_u_ub(1)*ones(size(T)),'r--','LineWidth',1.5)
    xlabel('time [s]')
    ylabel(sprintf('control input, u_%i^p [p.u.]', sys_idx))
    set(gca,'FontSize',fontsize)
    
%     % active power
%     figure
%     hold on
%     grid on
%     stairs(T, PQval(:,1),'k-','LineWidth',1.5)
%     plot(T, P_min*ones(size(T)),'r--','LineWidth',1.5)
%     plot(T, P_max*ones(size(T)),'r--','LineWidth',1.5)
%     xlabel('time [s]')
%     ylabel(sprintf('active power, P_%i [p.u.]', sys_idx))
%     set(gca,'FontSize',fontsize)
%     
%     % reactive power
%     figure
%     hold on
%     grid on
%     stairs(T, PQval(:,2),'k-','LineWidth',1.5)
%     plot(T, Q_min*ones(size(T)),'r--','LineWidth',1.5)
%     plot(T, Q_max*ones(size(T)),'r--','LineWidth',1.5)
%     xlabel('time [s]')
%     ylabel(sprintf('reactive power, Q_%i [p.u.]', sys_idx))
%     set(gca,'FontSize',fontsize)    
    
    
    
    neighbors = SysDesc(sys_idx).nei(2:end);
        
    
    for nei_id = 1:length(yval)
        
        nei = neighbors(nei_id);
        
        if nei > 1.1
            % phase angle plot
            figure
            hold on
            grid on
            plot(T, yval{nei_id}(:,end-2),'LineWidth',1.5)
            plot(T, options.lb_theta*ones(size(T)),'r--','LineWidth',1.5)
            plot(T, options.ub_theta*ones(size(T)),'r--','LineWidth',1.5)
            xlabel('time [s]')
            ylabel(sprintf('neighbor phase angle \x03B8_%i [rad]', nei))
            yticks([options.lb_theta 0 options.ub_theta])
            yticklabels({'-\pi/6','0','\pi/6'})
            set(gca,'FontSize',fontsize)
        end
        
        
        % voltage plot
        figure
        hold on
        grid on
        plot(T, yval{nei_id}(:,end),'LineWidth',1.5)
        plot(T, options.lb_v*ones(size(T)),'r--','LineWidth',1.5)
        plot(T, options.ub_v*ones(size(T)),'r--','LineWidth',1.5)
        xlabel('time [s]')
        ylabel(sprintf('(shifted) neighbor stochastic voltage v_%i [p.u.]', nei))
        set(gca,'FontSize',fontsize)
        
        % frequency plot
        figure
        hold on
        grid on
        plot(T, yval{nei_id}(:,end-1),'LineWidth',1.5)
        plot(T, options.lb_w*ones(size(T)),'r--','LineWidth',1.5)
        plot(T, options.ub_w*ones(size(T)),'r--','LineWidth',1.5)
        xlabel('time [s]')
        ylabel(sprintf('neighbor stochastic frequency \x03C9_%i [Hz]', nei))
        set(gca,'FontSize',fontsize)
    end
end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [tPoint,X,yval,uval,PQval] = runSim(SysDesc, sys_idx, resul_dir, options)
 

fprintf('\n*** Generating new conditions and disturbance for simulation.\n')   
tFinal = options.tFinal;
tSteps = options.tSteps;
tContr = options.tContr;
tDistu = options.tDistu;

 
if rem(tContr,tSteps)~=0
    error('control timestep has to be in multiple of the simulation timestep')
elseif rem(tDistu,tSteps)~=0
    error('disturbance timestep has to be in multiple of the simulation timestep')
end
 
tPoint = (0:tSteps:tFinal)';
 
% initialize the time-series vectors
X = NaN(length(tPoint),length(SysDesc(sys_idx).x));
uval = NaN(length(tPoint),2);
PQval = NaN(length(tPoint),2);
 
% Writing differential equation file
[fname_ode, fname_PQ] = genSimModel_np_nq(SysDesc, sys_idx, resul_dir, options);


%%%% neighbor nodes
neighbors = SysDesc(sys_idx).nei(2:end);
nb_nei = length(neighbors);

yval = cell(nb_nei,1);
 
varstr = '';
valstr = '';
for i_sys = 1:length(SysDesc(sys_idx).nei)
    if i_sys == 1
        for i_var = 1:length(SysDesc(sys_idx).x)
            varstr = strcat(varstr,'SysDesc(',string(sys_idx),').x(',string(i_var),'),');
            valstr = strcat(valstr,'X(iT,',string(i_var),'),');
        end
    else    % its neighbors
        nei_idx = SysDesc(sys_idx).nei(i_sys);
        for i_var = 1:length(SysDesc(nei_idx).x)
            varstr = strcat(varstr,'SysDesc(',string(nei_idx),').x(',string(i_var),'),');
            valstr = strcat(valstr,'yval{',string(i_sys-1),'}(iT,',string(i_var),'),');
        end
    end
end
 
oldDir = cd(resul_dir); % move to the location of the ODE-file

X(1,:) = genInitial(sys_idx, options); % initalize X

yTemp = genSamples(SysDesc, sys_idx, X(1,:), options);  % initialize Y
    
for i_nei = 1:nb_nei
    yval{i_nei} = NaN(length(tPoint),length(SysDesc(SysDesc(sys_idx).nei(i_nei+1)).x));
    yval{i_nei}(1,:) = yTemp{i_nei};
end

% compute control
if options.cflag == 3
    u0 = 0;
else
    error('Choose control = 3')
end
uval(1,:) = updateControl(1,SysDesc,X,yval,varstr,valstr,u0,options);

for iT = 2:length(tPoint)
    eval(sprintf('[tTemp,xTemp] = ode45(@(t,x) %s(t,x,yTemp,uval(iT-1,:)),tPoint(iT-1:iT),X(iT-1,:));',fname_ode))
    X(iT,:) = xTemp(end,:);
    if isnan(sum(X(iT,:)))
        error('reduce simulation time-step, or decrease the control effort')
    end
    
    [PQval(iT,1), PQval(iT,2)] = eval(sprintf('%s(X(iT,:), yTemp)', fname_PQ));
    
    if rem(iT-1,tDistu/tSteps) == 0    % update yTemp
        yTemp = genSamples(SysDesc, sys_idx, X(iT,:), options, yTemp);    % initialize Y
        for i_nei = 2:length(SysDesc(sys_idx).nei)
            yval{i_nei-1}(iT,:) = yTemp{i_nei-1};
        end
    else            % re-use yTemp from pevious time
        for i_nei = 2:length(SysDesc(sys_idx).nei)
            yval{i_nei-1}(iT,:) = yTemp{i_nei-1};
        end
    end

    % compute/update control

    if rem(iT-1,tContr/tSteps) == 0   % update uval
        uval(iT,:) = updateControl(iT,SysDesc,X,yval,varstr,valstr,u0,options);
    else    % re-use from previous time-step
        uval(iT,:) = uval(iT-1,:);
    end
end

cd(oldDir)
 
end
 




% Generate function sysodek with name of the system if not already existing
function [fname_ode, fname_PQ] = genSimModel_np_nq(SysDesc, sys_idx, resul_dir, options)

fname_ode = sprintf('sysode%i_np_%i_percent_nq_%i_percent', sys_idx, round(100*options.np_ratio), round(100*options.nq_ratio));
fname_PQ = sprintf('sys%i_PQ_np_%i_percent_nq_%i_percent', sys_idx, round(100*options.np_ratio), round(100*options.nq_ratio));
path = strcat(resul_dir,'/',fname_ode,'.m');

% if exist(path) == 2
%     return
% end

%%% Obtaining the dynamics
fprintf('Using originalDynamics and not sysode for the simulation\n')
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 np_ratio nq_ratio
% f = originalDynamics('sym');
f = tayDynamics('sym');
if sys_idx == 1
    x_dot = f(1:2);
else
    x_dot = f( 3*(sys_idx-1):3*sys_idx -1 );
end

%%% Obtaining the active power P
omega_id = 3*sys_idx - 1;
inv_tau_omega = eval(sprintf('subs(x_dot(end-1), [np_ratio, nq_ratio, x%i], [0, 0, -1])', omega_id));
dP = eval(sprintf('x_dot(end-1) + inv_tau_omega*x%i', omega_id));
dP = dP/inv_tau_omega;
dP = subs(dP, [np_ratio, nq_ratio], [1, 1]); % set lambda to 1 to remove it from x_dot as it does not belong to P
P = vpa(dP, 5);

%%% Obtaining the reactive power Q
v_id = omega_id + 1;
inv_tau_v = eval(sprintf('subs(x_dot(end), [np_ratio, nq_ratio, x%i], [0, 0, -1])', v_id));
dQ = eval(sprintf('x_dot(end) + inv_tau_v*x%i', v_id));
dQ = dQ/inv_tau_v;
dQ = subs(dQ, [np_ratio, nq_ratio], [1, 1]); % set lambda to 1 to remove it from x_dot as it does not belong to Q
Q = vpa(dQ, 5);

%%% State dynamics
x_dot = subs(x_dot, [np_ratio, nq_ratio], [options.np_ratio, options.nq_ratio]);
 
 
%%% Writing ODE file
fid = fopen(path,'w');
 
fprintf(fid,'function f = %s(t,x,y,u)\n', fname_ode);
 
fprintf(fid,'\n');
for i_var = 1:length(SysDesc(sys_idx).x)
    fprintf(fid,strcat(string(SysDesc(sys_idx).x(i_var)),' = x(',string(i_var),');\n'));
end
 
fprintf(fid,'\n');
for i_nei_idx = 2:length(SysDesc(sys_idx).nei)    % excluding self
    i_nei = SysDesc(sys_idx).nei(i_nei_idx);
    for i_var = 1:length(SysDesc(i_nei).x)
        fprintf(fid,strcat(string(SysDesc(i_nei).x(i_var)),' = y{',string(i_nei_idx-1),'}(',string(i_var),');\n'));
    end
end
 
fprintf(fid,'\nf = zeros(size(x));\n');
 
fprintf(fid,'\n');
for i_var = 1:length(SysDesc(sys_idx).x)
    fexpr = string( x_dot(i_var) );
    if i_var == length(SysDesc(sys_idx).x)-1
        fprintf(fid,strcat('f(',string(i_var),') = ',fexpr,' + u(1);\n'));
    elseif i_var == length(SysDesc(sys_idx).x)
        fprintf(fid,strcat('f(',string(i_var),') = ',fexpr,' + u(2);\n'));
    else
        fprintf(fid,strcat('f(',string(i_var),') = ',fexpr,';\n'));
    end
end




 
%%% Writing P and Q file
path = strcat(resul_dir,'/',fname_PQ,'.m');
fid = fopen(path,'w');
 
fprintf(fid,'function [P, Q] = %s(x,y)\n\n', fname_PQ);
 
for i_var = 1:length(SysDesc(sys_idx).x)
    fprintf(fid,strcat(string(SysDesc(sys_idx).x(i_var)),' = x(',string(i_var),');\n'));
end
 
fprintf(fid,'\n');
for i_nei_idx = 2:length(SysDesc(sys_idx).nei)    % excluding self
    i_nei = SysDesc(sys_idx).nei(i_nei_idx);
    for i_var = 1:length(SysDesc(i_nei).x)
        fprintf(fid,strcat(string(SysDesc(i_nei).x(i_var)),' = y{',string(i_nei_idx-1),'}(',string(i_var),');\n'));
    end
end

fprintf(fid,'\nP = %s;\nQ = %s;\n', string(P),string(Q));

end


 


% update control at each time step 
function uval = updateControl(iT,SysDesc,X,yval,varstr,valstr,u0,options)
 
switch options.cflag
    case 0      % no control
        uval = zeros(size(u0val));
    case 1      % computed control from CBF
        eval(strcat('uval = double(subs(u0,{',varstr,'},{',valstr,'}));'))
    case 2      % any other control input
        uval = -1 + 2*rand;
    case 3      % invariance control from Bony-Brezis
        c = char(varstr);
        sys_idx = str2num(c(9)); % 9th character of varstr is the id of the nod
        uval = invariance_control(sys_idx, X(iT, :), options);
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

% ratio in between bounds for the frequency and voltage of the node
ratio_v = (v - lb_v)/(ub_v - lb_v);
ratio_w = (w - lb_w)/(ub_w - lb_w);
% if min(ratio_v, ratio_w) < 0 || max(ratio_v, ratio_w) > 1
%     error('v or w is out of the safe set')
% end

% upper and lower bounds controls for v and w
u_w_lb = options.safe_u_lb(1);
u_w_ub = options.safe_u_ub(1);
u_v_lb = options.safe_u_lb(2);
u_v_ub = options.safe_u_ub(2);

%%%% Random control inside safety admissible range
if u_w_lb < u_w_ub && u_v_lb < u_v_ub
    uval(1) = u_w_lb + rand*(u_w_ub - u_w_lb);
    uval(2) = u_v_lb + rand*(u_v_ub - u_v_lb);
else
    error('No safety admissible range of controls.')
end

%%%% proportional control output % doesn't work great
% uval(1) = u_w_lb + ratio_w*(u_w_ub - u_w_lb);
% uval(2) = u_v_lb + ratio_v*(u_v_ub - u_v_lb);

%%%% single step control % chattering
% uval(1) = (ratio_w > 0.5)*u_w_ub + (ratio_w <= 0.5)*u_w_lb;
% uval(2) = (ratio_v > 0.5)*u_v_ub + (ratio_v <= 0.5)*u_v_lb;

%%%% double smooth step control % too complex
% y_w = double_smooth_step(ratio_w);
% uval(1) = (y_w > 0)* u_w_lb*y_w - (y_w <= 0)* u_w_ub*y_w;
% y_v = double_smooth_step(ratio_v);
% uval(2) = (y_v > 0)* u_v_lb*y_v - (y_v <= 0)* u_v_ub*y_v;

end

% Instead of y = 1 - 2x, here is a double smooth step function taking x in
% [0, 1] and giving him a value y in [-1, 1]
% see doubleSmoothStep.m
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

end



% Initalize the state x of the node
function x = genInitial(sys_idx, options)

% Random voltage in the bounds
v = options.lb_v + rand*(options.ub_v - options.lb_v);

% Random frequency in the bounds
w = 0;
% w = options.lb_w + rand*(options.ub_w - options.lb_w);

if sys_idx == 1
    x = [w, v];
else
    % Random angle in the bounds
    theta = options.lb_theta + rand*(options.ub_theta - options.lb_theta);
    x = [theta, w, v];
end

end



% Generate random samples for the states within some bounds
function y = genSamples(SysDesc, sys_idx, x, options, varargin)

if sys_idx == 1
    w = x(1); v = x(2);
else
    w = x(2); v = x(3);
end
neighbors = SysDesc(sys_idx).nei(2:end);
nb_nei = length(neighbors);
y = cell(nb_nei, 1);

% Angle inputs
lb_theta = options.lb_theta; % lower bound on the angle
ub_theta = options.ub_theta; % upper bound on the angle

% Coupling of voltage
coupling_v = options.coupling_v; % coupling of voltages with neighbors
% Coupling of frequency
coupling_w = options.coupling_w; % coupling of frequencies with neighbors

% Voltage inputs
lb_v = max(options.lb_v, v - coupling_v); % lower bound on the voltage
ub_v = min(options.ub_v, v + coupling_v); % upper bound on the voltage

% Frequency inputs
% lb_w = options.lb_w; % lower bound on the frequency
% ub_w = options.ub_w; % upper bound on the frequency
lb_w = max(options.lb_w, w - coupling_w); % lower bound on the frequency
ub_w = min(options.ub_w, w + coupling_w); % upper bound on the frequency


if options.varying_inputs
    
    for i = 1:nb_nei

        w_nei = lb_w + rand*(ub_w - lb_w);
        v_nei = lb_v + rand*(ub_v - lb_v);

        if neighbors(i) == 1
            y{i} = [w_nei, v_nei];
        else
            if isempty(varargin) % initialize neighbors
                theta_nei = 0.1*rand;
%                 theta_nei = lb_theta + rand*(ub_theta - lb_theta);
                
            else % update theta as integral of omega
                yTemp = varargin{1};
                old_theta = yTemp{i}(1);
                theta_nei = old_theta + w_nei*options.tSteps;
                if theta_nei > ub_theta
                    theta_nei = ub_theta;
                elseif theta_nei < lb_theta
                    theta_nei = lb_theta;
                end
            end
            y{i} = [theta_nei, w_nei, v_nei];
        end
    end
elseif options.worst_inputs
  
    worst_w = 0;
    worst_v = lb_v;
    worst_theta = lb_theta;
    
    for i = 1:nb_nei
        if neighbors(i) == 1
            y{i} = [worst_w, worst_v];
        else
            y{i} = [worst_theta, worst_w, worst_v];
        end
    end
end

end