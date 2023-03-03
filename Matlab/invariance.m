%%% Script to determine the range of safe controls to keep a state in some
%%% safe set. We calculate the minimal control needed so that the lower
%%% bound is invariant, and the maximal control making the upper bound
%%% invariant. These 2 range of controls might not necessarily intersect,
%%% yielding no feasible constant control making the safe set invariant.


% clear variables; 
% clc;
% close all

% directory where the model, and results are stored
resul_dir = sprintf('%s/results',pwd);

%%% LOAD/READ the necessary model files, and control barrier functions

% SysDesc: a structure array containing the interconnected network model,
% with each entry in SysDesc representing one subsystem
%
% SysDesc(i).x: set of symbolic variables associated with subsystem-i
% SysDesc(i).f: set of symbolic ODEs associated with subsystem-i
% SysDesc(i).nei: index set of the neighbors to subsystem-i (INCL. i)

droop_percentage = 100;
fprintf('Droop coefficients are at %i%% of the originals.\n', droop_percentage)
load(sprintf('%s/SysModel_%i_percent.mat',resul_dir,droop_percentage),'x','f','SysDesc')
  
 
%%% PERFORM simulations to test different control policies

% Lower and Upper bounds on the voltage
options.lb_v = -0.4;
options.ub_v = 0.2; 

% Lower and Upper bounds on the voltage of the neighbors
options.lb_v_nei = options.lb_v;
options.ub_v_nei = options.ub_v;
options.nb_pt_v = 31; % number of points to test for v_nei
% Lower and Upper bounds on the frequency
options.lb_w = -0.5;
options.ub_w = 0.5; 
% Lower and Upper bounds on the frequency of neighbors
options.lb_w_nei = options.lb_w;
options.ub_w_nei = options.ub_w; 
options.nb_pt_w = 10;
% number of points to test for theta in [-pi/2, pi/2]
options.lb_theta = -pi/2;
options.ub_theta = pi/2;
options.nb_pt_theta = 20; 

% v_nei is in the range [v - coupling, v + coupling].
options.coupling_v = 0.05;
% w_nei is in the range [w - coupling, w + coupling].
options.coupling_w = 0.05;


%%%% Safe controls upper and lower bounds for the whole system for w in [-0.1 0.1], v in [-0.4, 0.2]
% Safe_controls_ub = [-181.878947904333,-5.31442471242917,0;0,-86.1469143007416,-2.51587816051324;0,-319.431909354404,-4.04772784064699;0,-275.843945916742,-1.66274765365961];
% Safe_controls_lb = [336.392376461998,24.8743677933817,0;0,167.288937195204,11.2651646246086;0,388.318909179419,103.813042751048;0,308.883668094873,88.2148208972067];

%%%% Safe controls upper and lower bounds for the whole system for w in [-0.5 0.5], v in [-0.4, 0.2]
% Safe_controls_ub = [-181.024472783785,-5.31442471242917,0;0,-85.3770072478801,-2.51587816051324;0,-318.660837336443,-4.04772784064699;0,-275.112240413204,-1.66274765365961];
% Safe_controls_lb = [335.537901422940,24.8743677933817,0;0,166.519030130003,11.2651646246086;0,387.547836768441,103.813042751048;0,308.151962625794,88.2148208972067];

%%%% Safe controls upper and lower bounds for the whole system for w in [-0.1 0.1], v in [-0.1, 0.1]
% Safe_controls_ub = [-143.073727238458,-4.56735035090242;-67.9072685674764,-2.30345345211390;-159.440387082286,-3.42912558474200;-118.218881034758,-1.43063953227465];
% Safe_controls_lb = [267.603929264471,108.973033480113;133.198852762580,52.8184712365037;298.141765632667,219.054467986338;261.601674552076,192.650958425831];

%%%% Different dynamics: tayDynamics with np = nq = 1, and w in [-0.5 0.5], v in [-0.4, 0.2]
% Safe_controls_ub = [-19.4599581848015,-0.0937229651169673;-8.87530147968209,0.127143273620732;-27.3016570758773,-0.0781385527543534;-31.5317962629488,0.170167335600127];
% Safe_controls_lb = [36.8846904061502,1.46892850240329;18.1487937362981,0.325318770238482;33.3798592414241,9.61801457643742;35.4282163357129,9.61934458979522];


% choose an inverter node/bus for simulation
for sys_idx = 1    % OPTIONS: [1], [2], [3] or [4]
    % choose the variable to study
    for var_idx = 2 % for node 1: 1 = frequency, 2 = voltage
        % for other nodes: 1 = angle, 2 = frequency, 3 = voltage



        if (sys_idx == 1 && var_idx == 1) || (sys_idx >= 2 && var_idx == 2) % frequency case
            fprintf('Safe set is w in [%0.2f %0.2f] with a coupling of %0.2f for node %i\n',  options.lb_w, options.ub_w, options.coupling_w, sys_idx);   
        elseif (sys_idx == 1 && var_idx == 2) || (sys_idx >= 2 && var_idx == 3) % voltage case
            fprintf('Safe set is v in [%0.2f %0.2f] with a coupling of %0.2f for node %i\n',  options.lb_v, options.ub_v, options.coupling_v, sys_idx);
        else
            error('Need to select either voltage or frequency.')
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic();
% [f_min, f_max] = is_invariant(SysDesc, sys_idx, resul_dir, options, droop_percentage, var_idx);
% toc()

%%% takes two control signals u_lb and u_ub, with u_ub applied when v1 = ub_v
%%% while u_lb is applied when v1 = lb_v, whatever the states of the neighbor
% u_lb = -f_min + 1;
% u_ub = -f_max - 1;
% verif(SysDesc, sys_idx, resul_dir, options, droop_percentage, u_lb, u_ub)

        tic();
        [min_u_lb, max_u_ub] = sos_constraint(SysDesc, sys_idx, options, var_idx);
        toc()
        fprintf('SOS approach\n')
        fprintf('The set is invariant from below for all u > %f.\n', min_u_lb)
        fprintf('The set is invariant from above for all u < %f.\n\n', max_u_ub)
        
        
%         Safe_controls_ub(sys_idx, var_idx) = max_u_ub;
%         Safe_controls_lb(sys_idx, var_idx) = min_u_lb;
%         save('Safe_controls.mat', 'Safe_controls_ub', 'Safe_controls_lb');
        
    end
end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Generates the set of feasible controls for which the voltage safe set is
% invariant from below and above whatever the values of the states
function [f_min, f_max] = is_invariant(SysDesc, sys_idx, resul_dir, options, droop_percentage, var_idx)

% Inputs
if (sys_idx == 1 && var_idx == 1) || (sys_idx >= 2 && var_idx == 2) % frequency case
    lb_var = options.lb_w; % lower bound on the frequency
    ub_var = options.ub_w; % upper bound on the frequency
    
    variable_offset = 2;
elseif (sys_idx == 1 && var_idx == 2) || (sys_idx >= 2 && var_idx == 3) % voltage case
    lb_var = options.lb_v; % lower bound on the voltage
    ub_var = options.ub_v; % upper bound on the voltage
    
    variable_offset = 3;
end

coupling_w = 1.01*options.coupling_w; % coupling of frequency with neighbors
coupling_v = 1.01*options.coupling_v; % coupling of voltages with neighbors

neighbors = SysDesc(sys_idx).nei(2:end);
nb_nei = length(neighbors);

% Generates states for the subsystem studied
[~, X] = states(sys_idx, var_idx, options);
% Generates states for the neighbor subsystems
[Nb_nei_states, Y_nei] = neighbor_states(SysDesc, sys_idx, options);

% name of the ode function
fname = genSimModel(SysDesc, sys_idx, resul_dir, droop_percentage);
 
oldDir = cd(resul_dir); % move to the location of the ODE-file


%%%% Determine set of controls making safe set invariant

f_min = Inf;
f_max = -Inf;
for x = X
    var = x(var_idx);
    
    for i = 1:Nb_nei_states
        
% in case of coupling, only consider the neighbor states at +- coupling distance
        if ~isempty(coupling_v)
            nei_v = NaN(1, nb_nei);
            for j = 1:nb_nei
                nei_v(j) = Y_nei( 3*j, i);
            end
            if sum(abs(nei_v - x(end)) > coupling_v) > 0.5
                continue
            end
        end
        if ~isempty(coupling_w)
            nei_w = NaN(1, nb_nei);
            for j = 1:nb_nei
                nei_w(j) = Y_nei( 3*(j-1) + 2, i);
            end
            if sum(abs(nei_w - x(end-1)) > coupling_w) > 0.5
                continue
            end
        end
        
        y = cell(1, nb_nei);
        j = 1;
        for nei = neighbors
            if nei == 1
                y{j} = Y_nei( 3*(j-1) + 2:3*j, i);
            else
                y{j} = Y_nei( 3*(j-1) + 1:3*j, i);
            end
            j = j + 1;
        end
        
        
        % Calculate f = sysode...(0,x,Y_neighbor, [0,0]) to obtain dv/dt = f(end)
        eval(sprintf('f = %s(0,x,y,[0;0]);',fname))
        
        if var == lb_var && f(var_idx) < f_min
            f_min = f(var_idx);
%             min_nei = Y_nei(:,i);
                
        elseif var == ub_var && f(var_idx) > f_max
            f_max = f(var_idx);
%             max_nei = Y_nei(:,i);
        end
    end
end
% If f_min < 0, then we need u > -f_min to make the set invariant by below.
% If f_min > 0, then the set is invariant for u = 0, and we can
% even take u > -f_min and keep it invariant from below.

% If f_max > 0, then we need u < -f_max to make the set invariant by above.
% If f_max < 0, then the set is invariant for u = 0, and we can
% even take u < -f_max and keep it invariant from above.

if -f_min < -f_max
    fprintf('Safe set is invariant for u in [%f, %f].\n', -f_min, -f_max)
else
    fprintf('The set is invariant from below for all u > %f.\nThe set is invariant from above for all u < %f.\n', -f_min, -f_max)
end

cd(oldDir)
% min_nei
% max_nei
end
 

% Given two control signals u_lb and u_ub, with u_ub applied when v1 = ub_v
% while u_lb is applied when v1 = lb_v, whatever the states of the neighbor
% verify if the safe set is invariant
function verif(SysDesc, sys_idx, resul_dir, options, droop_percentage, u_lb, u_ub)

% Inputs
lb_v = options.lb_v; % lower bound on the voltage
ub_v = options.ub_v; % upper bound on the voltage
coupling = 1.01*options.coupling_v; % coupling of volatages with neighbors

nb_nei = length(SysDesc(sys_idx).nei)-1;

% Generates states for the subsystem studied
[~, X] = states(sys_idx, options);
% Generates states for the neighbor subsystems
[Nb_nei_states, Y_nei] = neighbor_states(SysDesc, sys_idx, options);

invariant = true;


% Some necessary stuff
fname = genSimModel(SysDesc, sys_idx, resul_dir, droop_percentage);
 
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


for x = X
    v = x(end);
    
    for i = 1:Nb_nei_states
        
        % in case of coupling, only consider the neighbor states where v_nei is at +- coupling from v
        if ~isempty(coupling)
            v_nei = NaN(1, nb_nei);
            for j = 1:nb_nei
                v_nei(j) = Y_nei(3*j, i);
            end
            if sum(abs(v_nei - v) > coupling) > 0.5
                continue
            end
        end
        
        y = cell(1, nb_nei);
        
        for j = 1:nb_nei
            y{j} = Y_nei( 3*(j-1) + 1:3*j, i);
        end   
        
        if v == lb_v 
            
            eval(sprintf('f = %s(0,x,y,[0;u_lb]);',fname))
            dv_dt = f(end);
            if dv_dt <= 0
                invariant = false;
            end
            
        elseif v == ub_v 
            
            eval(sprintf('f = %s(0,x,y,[0;u_ub]);',fname))
            dv_dt = f(end);
            if dv_dt >= 0
                invariant = false;
            end
        end
        
        if ~invariant
            break
        end
    end
    if ~invariant
        break
    end
end


if invariant
    fprintf('The safe set is invariant for u_lb = %0.2f and u_ub = %0.2f\n', u_lb, u_ub)
else
    fprintf('The safe set is not invariant for u_lb = %0.2f and u_ub = %0.2f\n', u_lb, u_ub)
end

cd(oldDir)

end
 


% Generates all states of the subsystem, with v taking only its extreme
% values
function [Nb_states, X] = states(sys_idx, var_idx, options)

% Voltage inputs
lb_v = options.lb_v_nei; % lower bound on the voltage of neighbors
ub_v = options.ub_v_nei; % upper bound on the voltage of neighbors
N_v = options.nb_pt_v; % number of points to test for v
% Frequency inputs
lb_w = options.lb_w; % lower bound on the frequency of neighbors
ub_w = options.ub_w; % upper bound on the frequency of neighbors
N_w = options.nb_pt_w; % number of points to test for w
% Angle inputs
lb_theta = options.lb_theta; % lower bound on the angle
ub_theta = options.ub_theta; % upper bound on the angle
N_theta = options.nb_pt_theta; % number of points to test for theta

% Range of values per state components
phase_angle = linspace(lb_theta, ub_theta, N_theta);

if mod(var_idx, 3) == 0 % voltage study
    voltage = [lb_v, ub_v];
    frequency = linspace(lb_w, ub_w, N_w);
    
    if sys_idx == 1
        Nb_states = 2*N_w;
    else
        Nb_states = 2*N_w*N_theta;
    end
else % frequency study
    voltage = linspace(lb_v, ub_v, N_v);
    frequency = [lb_w, ub_w];
    
    if sys_idx == 1
        Nb_states = 2*N_v;
    else
        Nb_states = 2*N_v*N_theta;
    end
end


if sys_idx == 1
    X = zeros(2, Nb_states);
    i = 1;
    for w = frequency
        for v = voltage
            X(:,i) = [w; v];
            i = i + 1;
        end
    end
else
    X = zeros(3, Nb_states);
    i = 1;
    for theta = phase_angle
        for w = frequency
            for v = voltage
                X(:,i) = [theta; w; v];
                i = i + 1;
            end
        end
    end
end


end




% Generates all neighbor states
function [Nb_nei_states, Y_nei] = neighbor_states(SysDesc, sys_idx, options)

neighbors = SysDesc(sys_idx).nei(2:end); % first of the list is current node
nb_nei = length(neighbors); % number of neighbors
if nb_nei == 0
    Nb_nei_states = 0;
    Y_nei = [];
    return
end

% Voltage inputs
lb_v = options.lb_v_nei; % lower bound on the voltage of neighbors
ub_v = options.ub_v_nei; % upper bound on the voltage of neighbors
N_v = options.nb_pt_v; % number of points to test for v
% Frequency inputs
lb_w = options.lb_w; % lower bound on the frequency of neighbors
ub_w = options.ub_w; % upper bound on the frequency of neighbors
N_w = options.nb_pt_w; % number of points to test for w
% Angle inputs
lb_theta = options.lb_theta; % lower bound on the angle
ub_theta = options.ub_theta; % upper bound on the angle
N_theta = options.nb_pt_theta; % number of points to test for theta

coupling_v = 1.01*options.coupling_v; % coupling of voltages with neighbors
coupling_w = 1.01*options.coupling_w; % coupling of frequency with neighbors

% Range of values per state components
voltage = linspace(lb_v, ub_v, N_v);
frequency = linspace(lb_w, ub_w, N_w);
phase_angle = linspace(lb_theta, ub_theta, N_theta);

% remove all voltages that are further than coupling distance from the lower or upper bound of v
if ~isempty(coupling_v)
    to_keep = (abs(voltage - options.lb_v) < coupling_v) + (abs(voltage - options.ub_v) < coupling_v);
    to_keep = (to_keep > 0.5); % if some voltage is within coupling distance of both the lower and upper bound of v, to_keep = 2, brought back to 1
    coupled_voltage = zeros(1, sum(to_keep));
    id_kept = 0;
    for id = 1:N_v
        if to_keep(id)
            id_kept = id_kept + 1;
            coupled_voltage(id_kept) = voltage(id); 
        end
    end
    voltage = coupled_voltage;
    N_v = id_kept;
end
% remove all frequencies that are further than coupling distance from the lower or upper bound of w
if ~isempty(coupling_w)
    to_keep = (abs(frequency - options.lb_w) < coupling_w) + (abs(frequency - options.ub_w) < coupling_w);
    to_keep = (to_keep > 0.5); % if some frequency is within coupling distance of both the lower and upper bound of w, to_keep = 2, brought back to 1
    coupled_frequency = zeros(1, sum(to_keep));
    id_kept = 0;
    for id = 1:N_w
        if to_keep(id)
            id_kept = id_kept + 1;
            coupled_frequency(id_kept) = frequency(id); 
        end
    end
    frequency = coupled_frequency;
    N_w = id_kept;
end



% Generates all possible states for all the neighbors
Nb_nei_states = 1;

for nei = neighbors
    if nei == 1 % Generates states for neighbor 1
        Nb_states_1 = N_v*N_w;
        Nb_nei_states = Nb_nei_states * Nb_states_1;
        Y_1 = zeros(3, Nb_states_1);
        i = 1;
        for w = frequency
            for v = voltage
                Y_1(:,i) = [0; w; v];
                i = i + 1;
            end
        end
        
    else % Generates states for a neighbor not 1
        eval(sprintf('Nb_states_%i = N_v*N_w*N_theta;', nei))
        eval(sprintf('Nb_nei_states = Nb_nei_states * Nb_states_%i;', nei))
        eval(sprintf('Y_%i = zeros(3, Nb_states_%i);', nei, nei))
        i = 1;
        for theta = phase_angle
            for w = frequency
                for v = voltage
                    eval(sprintf('Y_%i(:,i) = [theta; w; v];', nei))
                    i = i + 1;
                end
            end
        end
    end
end



if nb_nei == 1
    eval(sprintf('Y_nei = Y_%i;', neighbors(1) ))
    
elseif nb_nei == 2
    Y_nei = zeros(3*nb_nei, Nb_nei_states);
    
    for i = 1:eval(sprintf('Nb_states_%i', neighbors(1) ))
        for j = 1:eval(sprintf('Nb_states_%i', neighbors(2) ))
            % repeat column i of Y_neighbor(1)
            eval(sprintf('Y_nei(1:3, (i-1)*Nb_states_%i + j) = Y_%i(:,i);', neighbors(2), neighbors(1) ))
        end
        eval(sprintf('Y_nei(4:6, (i-1)*Nb_states_%i+1:i*Nb_states_%i) = Y_%i;', neighbors(2), neighbors(2), neighbors(2) ))
    end
else
    error('Too many neighbors.')
end

end



% Write the functions sysodej_k_percent
function fname = genSimModel(SysDesc, sys_idx, resul_dir, droop_percentage)
 
fname = sprintf('sysode%i_%i_percent', sys_idx, droop_percentage);
path = strcat(resul_dir,'/',fname,'.m');

if exist(path) == 2
    return
end

fid = fopen(path,'w');
 
fprintf(fid,sprintf('function f = %s(t,x,y,u)\n', fname));
 
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
    fexpr = string(sum(SysDesc(sys_idx).f(i_var,:)));
    if i_var == length(SysDesc(sys_idx).x)-1
        fprintf(fid,strcat('f(',string(i_var),') = ',fexpr,' + u(1);\n'));
    elseif i_var == length(SysDesc(sys_idx).x)
        fprintf(fid,strcat('f(',string(i_var),') = ',fexpr,' + u(2);\n'));
    else
        fprintf(fid,strcat('f(',string(i_var),') = ',fexpr,';\n'));
    end
end
 
end


% var_idx is the id of the variable of sys_idx to be kept in safe set,
% var_idx in {1,2,3}.
% Determine the upper bound and lower of the controls needed to keep state
% in safe set.
function [min_u_lb, max_u_ub] = sos_constraint(SysDesc, sys_idx, options, var_idx)

% x_id in [[1,12]] refers to the x variable of interest
if sys_idx == 1
    x_id = var_idx + 1;
else
    x_id = 3*(sys_idx-1) + var_idx;
end

% x_dot = sum(SysDesc(sys_idx).f(var_idx,:)); % differential equation
fprintf('Using tayDynamics and not sysode in sos_constraint\n')
f = tayDynamics('val',1,1);
x_dot = f(x_id-1);



%%%% neighbor nodes
neighbors = SysDesc(sys_idx).nei(2:end); % first of the list is current node
nb_nei = length(neighbors);

%%% Voltage inputs
lb_v = options.lb_v; % lower bound on the voltage
ub_v = options.ub_v; % upper bound on the voltage
lb_v_nei = options.lb_v_nei; % lower bound on the voltage of neighbors
ub_v_nei = options.ub_v_nei; % upper bound on the voltage of neighbors
coupling_v = options.coupling_v; % coupling of voltages with neighbors

%%% Frequency inputs
lb_w = options.lb_w; % lower bound on the frequency
ub_w = options.ub_w; % upper bound on the frequency
lb_w_nei = options.lb_w_nei; % lower bound on the frequency of neighbors
ub_w_nei = options.ub_w_nei; % upper bound on the frequency of neighbors
coupling_w = options.coupling_w; % coupling of frequency with neighbors

%%% Angle inputs
lb_theta = options.lb_theta; % lower bound on the angle
ub_theta = options.ub_theta; % upper bound on the angle


syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12
sos_options.params.fid = 0; % no print
% sos_options.solver = 'sdpa';
sos_options.solver = 'sedumi';

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
        ineq_ub = [ineq_ub, ineq_theta, eval(sprintf('x%i - (ub_v_nei-coupling_v)', v_nei)), eval(sprintf('ub_v_nei - x%i', v_nei)) ];
        ineq_lb = [ineq_lb, ineq_theta, eval(sprintf('x%i - lb_v_nei', v_nei)), eval(sprintf('(lb_v_nei+coupling_v) - x%i', v_nei)) ];
    end
    
    %%% upper bound of v
    v_dot_ub = subs(x_dot, eval(sprintf('x%i', x_id)), ub_v); % evaluate v at its upper bound
    [max_u_ub, ~, ~] = findbound(-v_dot_ub, ineq_ub, equality_constraints, 4, sos_options);
    
    
    %%% lower bound of v
    v_dot_lb = subs(x_dot, eval(sprintf('x%i', x_id)), lb_v); % evaluate v at its lower bound
    [min_u_lb, ~, ~] = findbound(v_dot_lb, ineq_lb, equality_constraints, 4, sos_options);
    min_u_lb = -min_u_lb;
    
else % frequency study
    % 4 cases here, because both v and w are coupled
    equality_constraints = [];
    
    %%% inequality constraints
    [ineq_theta, ineq_w_ub, ineq_w_lb, ineq_v_ub, ineq_v_lb] = deal([]);
    
    if sys_idx >= 2
        ineq_theta = [ eval(sprintf('x%i-lb_theta', x_id-1)), eval(sprintf('ub_theta-x%i', x_id-1)) ];
    end
    ineq_v_ub = [ineq_v_ub, eval(sprintf('x%i - (ub_v_nei-coupling_v)', x_id+1)), eval(sprintf('ub_v_nei - x%i', x_id+1)) ];
    ineq_v_lb = [ineq_v_lb, eval(sprintf('x%i - lb_v_nei', x_id+1)), eval(sprintf('(lb_v_nei+coupling_v) - x%i', x_id+1)) ];
    
    
    for nei = neighbors
        
        if nei >= 2
            theta_nei = 3*(nei-1) + 1; % id of the angle of the neighbor
            ineq_theta = [ineq_theta, eval(sprintf('x%i-lb_theta', theta_nei)), eval(sprintf('ub_theta-x%i', theta_nei)) ];
        end
        v_nei = nei*3; % id of the voltage of the neighbor
        ineq_v_ub = [ineq_v_ub, eval(sprintf('x%i - (ub_v_nei-coupling_v)', v_nei)), eval(sprintf('ub_v_nei - x%i', v_nei)) ];
        ineq_v_lb = [ineq_v_lb, eval(sprintf('x%i - lb_v_nei', v_nei)), eval(sprintf('(lb_v_nei+coupling_v) - x%i', v_nei)) ];
        
        w_nei = 3*(nei-1) + 2; % id of the frequency of the neighbor
        ineq_w_ub = [ineq_w_ub, eval(sprintf('x%i - (ub_w_nei-coupling_w)', w_nei)), eval(sprintf('ub_w_nei - x%i', w_nei)) ];
        ineq_w_lb = [ineq_w_lb, eval(sprintf('x%i - lb_w_nei', w_nei)), eval(sprintf('(lb_w_nei+coupling_w) - x%i', w_nei)) ];
        
    end
    
    %%% upper bound of w
    w_dot_ub = subs(x_dot, eval(sprintf('x%i', x_id)), ub_w); % evaluate w at its upper bound
    [max_u_ub_v_lb, ~, ~] = findbound(-w_dot_ub, [ineq_theta, ineq_w_ub, ineq_v_lb], equality_constraints, 4, sos_options);
    [max_u_ub_v_ub, ~, ~] = findbound(-w_dot_ub, [ineq_theta, ineq_w_ub, ineq_v_ub], equality_constraints, 4, sos_options);
    max_u_ub = min(max_u_ub_v_lb, max_u_ub_v_ub);
    
    %%% lower bound of w
    w_dot_lb = subs(x_dot, eval(sprintf('x%i', x_id)), lb_w); % evaluate w at its lower bound
    [min_u_lb_v_lb, ~, ~] = findbound(w_dot_lb, [ineq_theta, ineq_w_lb, ineq_v_lb], equality_constraints, 4, sos_options);
    [min_u_lb_v_ub, ~, ~] = findbound(w_dot_lb, [ineq_theta, ineq_w_lb, ineq_v_ub], equality_constraints, 4, sos_options);
    min_u_lb = max(-min_u_lb_v_lb, -min_u_lb_v_ub);    
end


end