clear variables; clc; close all

% directory where the model, and results are stored
resul_dir = sprintf('%s/results',pwd);


%% LOAD/READ the necessary model files, and control barrier functions 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% load model description (ODEs, variables, network interconnections, etc.)
%
% x: set of symbols denoting ALL system variables (angles, freq, voltages)
% f: set of symbolic ODEs representing the dynamics of each variable
%
% SysDesc: a structure array containing the interconnected network model,
% with each entry in SysDesc representing one subsystem
%
% SysDesc(i).x: set of symbolic variables associated with subsystem-i
% SysDesc(i).f: set of symbolic ODEs associated with subsystem-i
% SysDesc(i).nei: index set of the neighbors to subsystem-i (INCL. i)
% SysDesc(i).h: *** IGNORE ***
%
% NOTE that subsystem-1 is chosen as the reference bus for modeling
% purpose, which is why subsystem-1 has 2 state variables (frequency and
% votlage, while all other three subsystems have 3 states each (phase
% angle, frequency and voltage)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load(sprintf('%s/SysModel.mat',resul_dir),'x','f','SysDesc')

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% load subsystem barrier and the associated safety control functions
%
% barScaled: set of barrier functions for each subsystem, scaled to unity
% conFun:   a cell array of the active and reactive power controls for each
%           subsystem, with each conFun{i} having two symbolic entries -
%           the first entry corresponds to the active power control, while
%           the the second is for reactive power - both as functions of the
%           local state variables.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load(sprintf('%s/barrier.mat',resul_dir),'barScaled','conFun')


%% PERFORM time-domain simulations to test different control policies
 
% choose an inverter node/bus for simulation
sys_idx = 4;    % OPTIONS: [1], [2], [3] or [4]
% select control option
options.cflag = 1; % [0]: no control applied, [1]: safety control (pre-computed), [2]: arbitrary/custom control    

% simulation time-steps in seconds
options.tFinal = 2;     % total simulation duration in [s]
options.tSteps = 0.01;  % simulation time-step in [s]
options.tContr = 0.01;  % control changes at this time-step [s]
options.tDistu = 0.01;  % neighboring states change at this time-step [s]
 
% select whether to use pre-saved initial conditions and disturbance data
use_preSeededData = 0;  % [0]: new ICs and disturbances, [1]: pre-saved IC
 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NOW ready to run the time-domain simulations
%
% T: an array of time-stamps at the specified tSteps interval
% X: an array of the time-stamped values of the subsystem state-variables
% yval: an array of the time-stamped values of the neighbor states
% uval: an array of the time-stamped values of the applied control inputs
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if use_preSeededData == 0
    [T,X,yval,uval] = runSim(SysDesc,barScaled,sys_idx,conFun,resul_dir,options);
else
    % select prior results file index to re-use ICs and disturbances
    preSeed_fid = 1;
    [T,X,yval,uval] = runSim(SysDesc,barScaled,sys_idx,conFun,resul_dir,options,preSeed_fid);
end
 

%% SAVE AND PLOT results

if 0    % 1: save results, 0: don't save results
    result_idx = 1;
    while exist(strcat(resul_dir,'/sim_results',string(sys_idx),'_',string(result_idx),'.mat'),'file')
        result_idx = result_idx + 1;
    end
    save(strcat(resul_dir,'/sim_results',string(sys_idx),'_',string(result_idx),'.mat'),'T','X','yval','uval','options')
end
 
if 1    % plot time-domain state trajectories   
    % voltage plot
    figure
    plot(T,X(:,end),'LineWidth',1.5)
    hold on
    plot(T,0.2*ones(size(T)),'r--','LineWidth',1.5)
    plot(T,-0.4*ones(size(T)),'r--','LineWidth',1.5)
    xlabel('time [s]')
    ylabel('(shifted) voltage [p.u.]')
   
    % frequency plot
    figure
    plot(T,X(:,end-1),'LineWidth',1.5)
    xlabel('time [s]')
    ylabel('(shifted) frequency [Hz]')
    
    % reactive power control input
    figure
    stairs(T,uval(:,2),'k-','LineWidth',1.5)
    xlabel('time [s]')
    ylabel('control input, u^q [p.u.]')
   
    % active power control input
    figure
    stairs(T,uval(:,1),'k-','LineWidth',1.5)
    xlabel('time [s]')
    ylabel('control input, u^p [p.u.]')
end
 
 
%% FUNCTIONS
 
function [tPoint,X,yval,uval] = runSim(SysDesc,barScaled,sys_idx,conFun,resul_dir,options,varargin)
 
if nargin > 8 && exist(strcat(resul_dir,'/sim_results',string(sys_idx),'_',string(varargin{1}),'.mat'),'file') % re-run from a saved results directory
    new_options = options; % save the supplied options, before loading the saved results   
    use_saved_data = 1;
   
    fprintf('\n*** Loading initial conditions and disturbance from prior result file #%d ...\n',varargin{1})   
    load(strcat(resul_dir,'/sim_results',string(sys_idx),'_',string(varargin{1}),'.mat'),'X','yval','options')
    initX = X(1,:);
    ySave = yval;
   
    % update the saved timing parameters
    tFinal = options.tFinal;
    tSteps = options.tSteps;
    tContr = options.tContr;
    tDistu = options.tDistu;
    clear X yval options
   
    options = new_options;
else
    use_saved_data = 0;
    fprintf('\n*** Generating new conditions and disturbance for simulation.\n')   
    tFinal = options.tFinal;
    tSteps = options.tSteps;
    tContr = options.tContr;
    tDistu = options.tDistu;
end
 
if rem(tContr,tSteps)~=0
    error('control timestep has to be in multiple of the simulation timestep')
elseif rem(tDistu,tSteps)~=0
    error('disturbance timestep has to be in multiple of the simulation timestep')
end
 
tPoint = (0:tSteps:tFinal)';
 
% initialize the time-series vectors
X = NaN(length(tPoint),length(SysDesc(sys_idx).x));
 
uval = NaN(length(tPoint),2);
 
if ~use_saved_data
    [~,xSamples] = genSamples(SysDesc,barScaled,sys_idx);
end


u0 = conFun{sys_idx};
fname = genSimModel(SysDesc,sys_idx,resul_dir); 
yval = cell(length(SysDesc(sys_idx).nei)-1,1);
yTemp = cell(length(SysDesc(sys_idx).nei)-1,1);
 
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

if use_saved_data
    X(1,:)  = initX; % initatize X

    for i_nei = 2:length(SysDesc(sys_idx).nei)
        yval{i_nei-1} = NaN(length(tPoint),length(SysDesc(SysDesc(sys_idx).nei(i_nei)).x));
        yTemp{i_nei-1} = ySave{i_nei-1}(1,:);    % initialize Y
        yval{i_nei-1}(1,:) = yTemp{i_nei-1};
    end

    % compute control bounds
    uval(1,:) = updateControl(1,SysDesc,X,yval,varstr,valstr,u0,options);

    for iT = 2:length(tPoint)

        eval(sprintf('[tTemp,xTemp] = ode45(@(t,x) %s(t,x,yTemp,uval(iT-1,:)),tPoint(iT-1:iT),X(iT-1,:));',fname))

        X(iT,:) = xTemp(end,:);
        if rem(iT-1,tDistu/tSteps)==0    % update yTemp
            for i_nei = 2:length(SysDesc(sys_idx).nei)
                yTemp{i_nei-1} = ySave{i_nei-1}(iT,:);    % initialize Y
                yval{i_nei-1}(iT,:) = yTemp{i_nei-1};
            end
        else            % re-use yTemp from pevious time
            for i_nei = 2:length(SysDesc(sys_idx).nei)
                yval{i_nei-1}(iT,:) = yTemp{i_nei-1};
            end
        end

        % compute/update control
        if rem(iT-1,tContr/tSteps)==0   % update uval
            uval(iT,:) = updateControl(iT,SysDesc,X,yval,varstr,valstr,u0,options);
        else    % re-use from previous time-step
            uval(iT,:) = uval(iT-1,:);
        end
    end
else
    X(1,:)  = xSamples{1}(randperm(size(xSamples{1},1),1),:); % initatize X
    for i_nei = 2:length(SysDesc(sys_idx).nei)
        yval{i_nei-1} = NaN(length(tPoint),length(SysDesc(SysDesc(sys_idx).nei(i_nei)).x));
        yTemp{i_nei-1} = xSamples{i_nei}(randperm(size(xSamples{i_nei},1),1),:);    % initialize Y
        yval{i_nei-1}(1,:) = yTemp{i_nei-1};
    end

    % compute control bounds
    uval(1,:) = updateControl(1,SysDesc,X,yval,varstr,valstr,u0,options);
   
    for iT = 2:length(tPoint)
        eval(sprintf('[tTemp,xTemp] = ode45(@(t,x) %s(t,x,yTemp,uval(iT-1,:)),tPoint(iT-1:iT),X(iT-1,:));',fname))
       
        X(iT,:) = xTemp(end,:);
        if isnan(sum(X(iT,:)))
            error('reduce simulation time-step, or decrease the control effort')
        end
        if rem(iT-1,tDistu/tSteps)==0    % update yTemp
            for i_nei = 2:length(SysDesc(sys_idx).nei)
                yTemp{i_nei-1} = xSamples{i_nei}(randperm(size(xSamples{i_nei},1),1),:);    % initialize Y
                yval{i_nei-1}(iT,:) = yTemp{i_nei-1};
            end
        else            % re-use yTemp from pevious time
            for i_nei = 2:length(SysDesc(sys_idx).nei)
                yval{i_nei-1}(iT,:) = yTemp{i_nei-1};
            end
        end

        % compute/update control
        if rem(iT-1,tContr/tSteps)==0   % update uval
            uval(iT,:) = updateControl(iT,SysDesc,X,yval,varstr,valstr,u0,options);
        else    % re-use from previous time-step
            uval(iT,:) = uval(iT-1,:);
        end
    end
end
cd(oldDir)

end

 

function fname = genSimModel(SysDesc,sys_idx,resul_dir)

fname = strcat('sysode',string(sys_idx));
fid = fopen(strcat(resul_dir,'/',fname,'.m'),'w');
fprintf(fid,strcat('function f = sysode',string(sys_idx),'(t,x,y,u)\n'));
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

 

function uval = updateControl(iT,SysDesc,X,yval,varstr,valstr,u0,options)

eval(strcat('u0val = double(subs(u0,{',varstr,'},{',valstr,'}));'))

switch options.cflag
    case 0      % no control
        uval = zeros(size(u0val));
    case 1      % computed control from CBF
        uval = u0val;
    case 2      % any other control input
        uval = -1 + 2*rand;
end

end

 

function [valrange,xSamples] = genSamples(SysDesc,barScaled,sys_idx)

valrange = cell(length(SysDesc(sys_idx).nei),1);
N = 1e4;
xSamples = cell(length(SysDesc(sys_idx).nei),1);

for i = 1:length(SysDesc(sys_idx).nei)

    i_sys = SysDesc(sys_idx).nei(i);   
    sys_var = SysDesc(i_sys).x;
    valrange{i} = NaN(length(sys_var),2);

    % generate the seed sample points
    sPoints = lhsdesign(N,length(sys_var));
    xPoints = NaN(size(sPoints));
    repStr = 'barsamples = double(subs(barScaled(i_sys),SysDesc(i_sys).x,{';

    for i_var = 1:length(sys_var)
        var_coeff = -diff(diff(barScaled(i_sys),sys_var(i_var)),sys_var(i_var))/2; % divide by 2 since it is quadratic
        valrange{i}(i_var,:) = [-1 1]/sqrt(var_coeff);

        % now generate the state-space samples
        xPoints(:,i_var) = valrange{i}(i_var,1) + diff(valrange{i}(i_var,:)) * sPoints(:,i_var);
        repStr = strcat(repStr,sprintf('xPoints(:,%d),',i_var));
    end
    eval(strcat(repStr,'}));'))

   
    % select the sample points that are inside the 0 level-set
    xSamples{i} = xPoints(barsamples.*(barsamples-0.05)<=0,:);

end

end