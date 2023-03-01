clear variables
clc
close all

% droop coefficients
np_value = 10*0.2*2*pi/5; % = 2.513 rad/s
nq_value = 10*0.1/5; % = 0.2 p.u.

%%%%%% S_theta = +-pi/2  S_omega = +-0.5  delta_v = 0.05
% lambda_p_star = np_value*[3.6, 6.8, 2.9, 3]/100
% lambda_q_star = nq_value*[44.8, 85.7, 10, 11.5]/100

%%%%%% S_theta = +-pi/6  S_omega = +-3  delta_v = 0.02
% lambda_p_star = np_value*[0.4884, 0.9267, 0.3483, 0.54412]
% lambda_q_star = nq_value*[1.1407, 2.172, 0.8055, 1.205]

%%%%%%%%%%% INPUTS %%%%%%%%%

node = 2;
voltage = 0;
frequency = 1;
theta_ub = pi/6;
coupling = 0.02;


%%%%%%%%%%%%%%%%%%%%%% Lambda* graph


figure(1)
hold on
grid on


coupling_range = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6];

if node == 1 && voltage
    % S_theta = [-pi/6, pi/6] and S_v = [-0.4, 0.2]
    lambda_star_pi6 = [1.1407, 1.1407, 1.1407, 1.1407, 1.1407, 1.1407, 1.1407, 1.1407, 1.1407];
    plot(coupling_range, nq_value*lambda_star_pi6, 'LineWidth', 3)

elseif node == 2 && voltage   
%     % S_theta = [-pi/6, pi/6] and S_v = [-0.4, 0.2]
%     lambda_star_pi6 = [2.1723, 2.1723, 2.1723, 2.1723, 2.1723, 2.1723, 2.1723, 2.1723, 2.1723];
%     plot(coupling_range, nq_value*lambda_star_pi6, 'LineWidth', 3)
    
    % S_theta = [-pi/2, pi/2] and S_v = [-0.4, 0.2]
    lambda_star_pi2 = [0.90047, 0.89149, 0.85730, 0.79853, 0.63735, 0.53030, 0.45405, 0.39696, 0.35263];
    plot(coupling_range, nq_value*lambda_star_pi2, 'LineWidth', 3)
     
%     % S_theta = [-pi, pi] and S_v = [-0.4, 0.2]
%     lambda_star_pi = [0.12637, 0.12439, 0.11705, 0.10900, 0.09583, 0.08550, 0.07718, 0.07034, 0.06461];
%     plot(coupling_range, nq_value*lambda_star_pi, 'LineWidth', 3)
%     legend('S_\theta = [-\pi/6, \pi/6]', 'S_\theta = [-\pi/2, \pi/2]', 'S_\theta = [-\pi, \pi]')
    
    
    % S_theta = [-pi/2, pi/2] and S_v = [-0.3, 0.1]
    lambda_star_shrink = [0.33877, 0.33368, 0.31478, 0.29396, 0.25962, 0.23246, 0.21045, 0.21045, 0.21045];
    plot(coupling_range, nq_value*lambda_star_shrink, 'LineWidth', 3)
    % S_theta = [-pi/2, pi/2] and S_v = [-0.2, 0]
    lambda_star_shrinker = [0.10960, 0.10859, 0.10472, 0.10026, 0.09239, 0.09239, 0.09239, 0.09239, 0.09239];
    plot(coupling_range, nq_value*lambda_star_shrinker, 'LineWidth', 3)
    legend('S_v = [-0.4, 0.2] p.u.', 'S_v = [-0.3, 0.1] p.u.', 'S_v = [-0.2, 0] p.u.')
    ylim(nq_value*[0.08 1])
    
    
    xlim([coupling_range(1) coupling_range(end)])
    ylabel('\lambda_2^q* [p.u./p.u.]')
end

xlabel('allowable neighbor uncertainty \Delta_v [p.u.]')
set(gca,'FontSize',20)

















%%%%%%%%% Graph of invariance as a function of droop-coefficients %%%%%%%%%

%%%%%%%%%%%%%%%% JULIA %%%%%%%%%%%%%%%


figure(2)
hold on
grid on

if node == 1 && voltage && theta_ub == pi/6 && coupling == 0.02
    %%%% Node 1 voltage  lambda^q
    MAX_U_UB = [0.2678285144492031 0.21617830438498234 0.1645280907265324 0.11287787761188409 0.06122766328276439 0.009577449313333205 -0.0420727660053422 -0.09372298117331723 -0.14537320357896108 -0.19702341359354242];
    MIN_U_LB = [-0.6669467727989197 -0.6074095948165482 -0.5478724201497553 -0.4883352415562904 -0.4287980748829158 -0.369260892079799 -0.3097237138739564 -0.25018654003614715 -0.1906493620624223 -0.13111218345300465];
    droop = nq_value*(0.3:0.1:1.2);
    lambda_star = nq_value*1.1407;
    safe_u = -0.1664055;
    ylim(nq_value*[28 122]/100)
    xlim([-0.7 0.3])
    ylabel('\lambda_1^q [p.u./p.u.]')
    
elseif node == 1 && voltage && theta_ub == pi/2 && coupling == 0.05
    %%%% Node 1 voltage  lambda^q
    MAX_U_UB = [0.37112894449034467,   0.3194787444719229,  0.26782852855973976,  0.21617831536726714,  0.16452811175531712,  0.11287788947065387,  0.06122768407841993,  0.009577474289089899, -0.042072740695855544, -0.09372295089836052];
    MIN_U_LB = [-0.61411  -0.382662  -0.151213  0.0802357  0.311685  0.543133  0.774582  1.00603  1.23748  1.46893];
    droop = nq_value*(0.1:0.1:1);
    lambda_star = nq_value*0.44812;
    safe_u = 0.191356;
    ylim(nq_value*[8 52]/100)
    xlim([-0.65 0.4])
    ylabel('\lambda_1^q [p.u./p.u.]')
    
elseif node == 1 && frequency && theta_ub == pi/2 && coupling == 0.05
    %%%% Node 1 frequency  lambda^p
    MAX_U_UB = [-0.9947424029818656 -3.046433134534633 -5.098123624936316 -7.149814288275068 -9.201504782613608 -11.25319545749717 -13.3048855869775 -15.356576689343473 -17.408266947974102 -19.459957688591498];
    MIN_U_LB = [2.73721559455594 6.5313797093393235 10.325543191796616 14.11970728844696 17.913870797697975 21.708034390477547 25.50219852926462 29.296362165188622 33.09052587052473 36.88468972889245];
    droop = np_value*(0.1:0.1:1);
    lambda_star = np_value*0.03614;
    safe_u = 0.31506;
    ylabel('\lambda_1^p [rad/s/p.u.]')
    
elseif node == 2 && voltage && theta_ub == pi/6 && coupling == 0.02
    %%%% Node 2 voltage  lambda^q
    MAX_U_UB = [0.34292730154883005 0.2620082814802607 0.18108925778230203 0.10017022595182916 0.019251231617535905 -0.061667785893922966 -0.1425868124802646 -0.2235058344171359];
    MIN_U_LB = [-0.738076984749844 -0.6545725194732801 -0.5710680492526542 -0.48756357650806087 -0.40405910904092884 -0.32055463933317957 -0.23705016834012527 -0.15354571121106989];
    droop = nq_value*(0.2:0.3:2.3);
    lambda_star = nq_value*2.1723;
    safe_u = -0.18907;
    ylim(nq_value*[18 232]/100)
    xlim([-0.75 0.35])
    ylabel('\lambda_2^q [p.u./p.u.]')
        
elseif node == 2 && frequency && theta_ub == pi/6 && coupling == 0.02
    %%%% Node 2 frequency  lambda^p
    MAX_U_UB = [5.368677182491231 4.784251847052688 4.199829322929603 3.6154034794394265 3.0309792757338903 2.446555119536775 1.8621309074374695 1.2777066234005174 0.693282483858554 0.10885798626502062];
    MIN_U_LB = [-5.252758815919596 -4.5524165162331 -3.85207387629451 -3.1517315696255133 -2.451389321362152 -1.751045913343945 -1.0507034191843372 -0.3503610911022062 0.34998131738362476 1.0503237260778568];
    droop = np_value*(0.1:0.1:1);
    lambda_star = np_value*0.9267;
    safe_u = 0.5371;
    ylim(np_value*[8 102]/100)
    xlim([-5.5 5.5])
    ylabel('\lambda_2^p [rad/s/p.u.]')
    
elseif node == 2 && voltage && theta_ub == pi/2 && coupling == 0.05
    %%%% Node 2 voltage  lambda^q
    MIN_U_LB = [-0.6818402848500731 -0.5699338131473111 -0.4580272984540959 -0.34612074203124876 -0.23421423536688643 -0.12230771557037384 -0.010401188805245855 0.10150534561667246 0.21341187122073518 0.32531839129842277];
    MAX_U_UB = [0.3699003258295095 0.34292731894855827 0.31595431299351145 0.2889813058765005 0.26200829982297597 0.2350352938297018 0.20806228789493164 0.18108928197731997 0.15411627602113198 0.12714327024683678];
    droop = nq_value*(0.1:0.1:1);
    lambda_star = nq_value*0.85717;
    safe_u = 0.1654;
    xlim([-0.7 0.4])
    ylim(nq_value*[8 92]/100)
    ylabel('\lambda_2^q [p.u./p.u.]')
     
elseif node == 2 && frequency && theta_ub == pi/2 && coupling == 0.05
    %%%% Node 2 frequency  lambda^p
    MAX_U_UB = [0.8935084958842141 0.7948336401652255 0.6961588015991969 0.5974839595936282 0.498809115776717 0.4001342649102642 0.30145942050927377 0.0 0.0 0.0];
    MIN_U_LB = [-0.8007736629391603 -0.6093640593802034 -0.41795426718001216 -0.22654447695374266 -0.035134627798260466 0.1562751877605419 0.3476849555000219 0.0 0.0 0.0];
    droop = np_value*(0.01:0.01:0.1);
    lambda_star = np_value*0.0684;
    safe_u = 0.3172480714;
    xlim([-0.9 1])
    ylim(np_value*[0.5 7.5]/100)
    ylabel('\lambda_2^p [rad/s/p.u.]')
      
end



i = 1;
min_u_lb = MIN_U_LB(i);
max_u_ub = MAX_U_UB(i);
    
while i < length(MIN_U_LB) && min_u_lb < max_u_ub
    
    plot([min_u_lb, max_u_ub], [droop(i) droop(i)], 'blue', 'LineWidth', 3)
    
    i = i + 1;
    min_u_lb = MIN_U_LB(i);
    max_u_ub = MAX_U_UB(i);
end
plot( [MAX_U_UB(1) max_u_ub], [droop(1) droop(i)], '-.green', 'LineWidth', 3)
plot( [MIN_U_LB(1) min_u_lb], [droop(1) droop(i)], '-.green', 'LineWidth', 3)

plot([max_u_ub min_u_lb], [droop(i) droop(i)], ':r', 'LineWidth', 3)
scatter(safe_u, lambda_star, 20, 'blue', 'filled')


xlabel('safety admissible controls [p.u.]')
set(gca,'FontSize',20)
