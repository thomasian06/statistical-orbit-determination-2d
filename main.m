%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 5044 Final Project
%
% Purpose:
%
% Author(s): Ian Thomas, 
%
% Last Modified: 12/9/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
clc;

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

set(0,'DefaultAxesLooseInset',[0,0,0,0])

colors = [0    0.4471    0.7412
          0.8510    0.3255    0.0980
          0.9294    0.6941    0.1255
          0.4941    0.1843    0.5569
          0.4667    0.6745    0.1882
          0.3020    0.7451    0.9333
          0.6353    0.0784    0.1843
          0.9882    0.6431    0.6431
          0.0588    1.0000    1.0000
          1.0000    0.0745    0.6510
          0     0     1
          0.5020    0.5020    0.5020];
      
format long

pos = [100, 100, 800, 600];

%% Part 1

% Declaration of Parameters
mu = 398600; % km^3/s^2
RE = 6378; % km

n = 4; % number of states
p = 3; % number of measurements

X0 = 6678; % km
Y0 = 0; % km
r0 = sqrt(X0^2+Y0^2); % km
dX0 = 0; % km/s
dY0 = r0*sqrt(mu/r0^3); % km/s

T = 2*pi*sqrt(r0^3/mu); % orbital period [s]
omega = 2*pi/T; % angular velocity of sc [rad/s]
omegaE = 2*pi/86400; % angular velocity of earth [rad/s]

Dt = 10; % time discretization [s] 

nS = 12; % number of stations

% Nominal Solution
thetanom = @(t) omega * t;
Xnom = @(t) r0*cos(thetanom(t));
Ynom = @(t) r0*sin(thetanom(t));
dXnom = @(t) -r0*omega*sin(thetanom(t));
dYnom = @(t) r0*omega*cos(thetanom(t));
xnom = @(t) [Xnom(t), dXnom(t), Ynom(t), dYnom(t)];

Xi = @(t, i) RE*cos(omegaE*t+(i-1)*pi/6);
Yi = @(t, i) RE*sin(omegaE*t+(i-1)*pi/6);
dXi = @(t, i) -RE*omegaE*sin(omegaE*t+(i-1)*pi/6);
dYi = @(t, i) RE*omegaE*cos(omegaE*t+(i-1)*pi/6);
xi = @(t) [Xi(t), dXi(t), Yi(t), dYi(t)];

% Measurements
rho = @(t, x, i) sqrt((x(:,1)-Xi(t, i)).^2+(x(:,3)-Yi(t,i)).^2);
drho = @(t, x, i) ((x(:,1)-Xi(t,i)).*(x(:,2)-dXi(t,i))+...
    (x(:,3)-Yi(t,i)).*(x(:,4)-dYi(t,i))) ./ rho(t,x,i);
phi = @(t, x, i) atan2((x(:,3)-Yi(t,i)),(x(:,1)-Xi(t,i)));
mask = @(t, x, i) (absangdiff(atan2((x(:,3)-Yi(t,i)),(x(:,1)-Xi(t,i))),...
    atan2(Yi(t,i),Xi(t,i))) < (pi/2));
% mask = @(t, x, i) ~(abs(angdiff(atan2((x(:,3)-yi(t,i)),(x(:,1)-xi(t,i))),...
%     atan2(yi(t,i),xi(t,i)))) < (pi/2));

% Initial Perturbation
d_x0 = [0, 0.075, 0, -0.021]'; % perturbation

%% Nonlinear Simulation

% State
T = 1401;
tspan = (0:T-1)*10;
% tspan = [0, 14000];
x0_nl = [X0, dX0, Y0, dY0]' + d_x0;
tol = 1e-12;
opts = odeset('RelTol', tol);
[t_nl, x_nl] = ode45(@(t, S) nl_sim_nonoise(t, S, mu),...
    tspan, x0_nl, opts);


F_nl_states = figure('DefaultAxesPosition', [0.05, 0.05, 0.95, 0.95]);
set(gcf, 'Position',  pos)
x_nl_ylabels = {'$X$ [km]','$\dot X$ [km/s]',...
    '$Y$ [km]','$\dot Y$ [km/s]'};
for i = 1:n
    subplot(n, 1, i);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(x_nl_ylabels{i})
    plot(t_nl, x_nl(:, i), 'LineWidth', 1.5)
    set(gca,'FontSize',12)
end
sgtitle('Nonlinear Simulated States')
saveas(gcf, '1', 'epsc')

rho_nl = rho(t_nl, x_nl, 1:nS);
drho_nl = drho(t_nl, x_nl, 1:nS);
phi_nl = phi(t_nl, x_nl, 1:nS);
mask_nl = mask(t_nl, x_nl, 1:nS);
mask_nl = permute(repmat(mask_nl, 1, 1, 3), [1, 3, 2]);

% make time first dimension, measurement index second dimension, and
% station index third dimension
y_nl = permute(cat(3, rho_nl, drho_nl, phi_nl), [1, 3, 2]);
y_nl(mask_nl) = NaN;

F_nl_meas = figure('DefaultAxesPosition', [0.05, 0.05, 0.95, 0.95]);
set(gcf, 'Position',  pos)
y_nl_ylabels = {'$\rho$ [km]','$\dot \rho$ [km/s]','$\phi$ [rad]'};
for j = 1:p
    ax = subplot(p, 1, j);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(y_nl_ylabels{j})
    for i = 1:nS
        plot(t_nl, y_nl(:, j, i), '.', 'Color', colors(i, :)); 
    end
    set(gca,'FontSize',12)
end
sgtitle('Nonlinear Simulated Measurements')
saveas(gcf, '2', 'epsc')


%% Linearized Simulation

% A tilde
A = @(a, r) [0, 1, 0, 0
          mu*(2*cos(a)^2-sin(a)^2)/r^3, 0, mu*(3*sin(2*a))/(2*r^3), 0
          0, 0, 0, 1
          mu*(3*sin(2*a))/(2*r^3),0,mu*(2*sin(a)^2-cos(a)^2)/r^3,0];
B = [0, 0; 1, 0; 0, 0; 0, 1]; % B tilde
Y = B; % Gamma tilde

% C tilde partial derivatives
C11 = @(t, x, i) (x(1)-Xi(t,i))/sqrt((x(1)-Xi(t,i))^2+(x(3)-Yi(t,i))^2);
C13 = @(t, x, i) (x(3)-Yi(t,i))/sqrt((x(1)-Xi(t,i))^2+(x(3)-Yi(t,i))^2);
C21 = @(t, x, i) (x(2)-dXi(t,i))*...
    sqrt((x(1)-Xi(t,i))^2+(x(3)-Yi(t,i))^2)/...
    ((x(1)-Xi(t,i))^2+(x(3)-Yi(t,i))^2) - ...
    ((x(1)-Xi(t,i))*(x(2)-dXi(t,i))+...
    (x(3)-Yi(t,i))*(x(4)-dYi(t,i)))*...
    (1/sqrt((x(1)-Xi(t,i))^2+(x(3)-Yi(t,i))^2)*...
    (x(1)-Xi(t,i)))/((x(1)-Xi(t,i))^2+(x(3)-Yi(t,i))^2);
C22 = @(t, x, i) (x(1)-Xi(t,i))*...
    sqrt((x(1)-Xi(t,i))^2+(x(3)-Yi(t,i))^2)/...
    ((x(1)-Xi(t,i))^2+(x(3)-Yi(t,i))^2);
C23 = @(t, x, i) (x(4)-dYi(t,i))*...
    sqrt((x(1)-Xi(t,i))^2+(x(3)-Yi(t,i))^2)/...
    ((x(1)-Xi(t,i))^2+(x(3)-Yi(t,i))^2) - ...
    ((x(1)-Xi(t,i))*(x(2)-dXi(t,i))+...
    (x(3)-Yi(t,i))*(x(4)-dYi(t,i)))*...
    (1/sqrt((x(1)-Xi(t,i))^2+(x(3)-Yi(t,i))^2)*...
    (x(3)-Yi(t,i)))/((x(1)-Xi(t,i))^2+(x(3)-Yi(t,i))^2);
C24 = @(t, x, i) (x(3)-Yi(t,i))*...
    sqrt((x(1)-Xi(t,i))^2+(x(3)-Yi(t,i))^2)/...
    ((x(1)-Xi(t,i))^2+(x(3)-Yi(t,i))^2);
C31 = @(t, x, i) -1/(1+((x(3)-Yi(t,i))/(x(1)-Xi(t,i)))^2)*...
    (x(3)-Yi(t,i))/(x(1)-Xi(t,i))^2;
C33 = @(t, x, i) 1/(1+((x(3)-Yi(t,i))/(x(1)-Xi(t,i)))^2)*...
    1/(x(1)-Xi(t,i));

% C tilde
C = @(t, x, i) [C11(t, x, i), 0, C13(t, x, i), 0
     C21(t, x, i), C22(t, x, i), C23(t, x, i), C24(t, x, i)
     C31(t, x, i), 0, C33(t, x, i), 0];

Fk = @(k) eye(n)+Dt*A(thetanom(k*Dt), r0);
Gk = Dt*B;
Omegak = Dt*Y;
Hk = @(k, i) C(k*Dt, xnom(Dt*k), i);


T = 1400;
d_x_lin = zeros(n, T+1);
d_xprev = d_x0;
k = (0:T)';

for i = 0:T
    d_x_lin(:, i+1) = d_xprev;
    d_xprev = Fk(i)*d_xprev;
end

d_x_lin = d_x_lin';
t_lin = k*Dt;

F_lin_delta = figure('DefaultAxesPosition', [0.05, 0.05, 0.95, 0.95]);
set(gcf, 'Position',  pos)
d_x_lin_ylabels = {'$\delta X$ [km]','$\delta \dot X$ [km/s]',...
    '$\delta Y$ [km]','$\delta \dot Y$ [km/s]'};
for i = 1:n
    subplot(n, 1, i);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(d_x_lin_ylabels{i})
    plot(t_lin, d_x_lin(:, i), 'LineWidth', 1.5)
    set(gca,'FontSize',12)
end
sgtitle('Linearized Perturbations')
saveas(gcf, '3', 'epsc')

x_nom = [Xnom(k*Dt), dXnom(k*Dt), Ynom(k*Dt), dYnom(k*Dt)];
x_lin = d_x_lin+x_nom;

F_lin_states = figure('DefaultAxesPosition', [0.05, 0.05, 0.95, 0.95]);
set(gcf, 'Position',  pos)
x_lin_ylabels = {'$X$ [km]','$\dot X$ [km/s]',...
    '$Y$ [km]','$\dot Y$ [km/s]'};
for i = 1:n
    subplot(n, 1, i);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(x_lin_ylabels{i})
    plot(t_lin, x_lin(:, i), 'LineWidth', 1.5)
    set(gca,'FontSize',12)
end
sgtitle('Linearized State Simulation')
saveas(gcf, '4', 'epsc')

d_y_lin = zeros(T+1, nS, p);

for j = 0:T
    for i = 1:nS
        d_y_lin(j+1, i, :) = Hk(j, i)*d_x_lin(j+1, :)';
    end
end

rho_nom = rho(t_lin, x_nom, 1:nS);
drho_nom = drho(t_lin, x_nom, 1:nS);
phi_nom = phi(t_lin, x_nom, 1:nS);
ynom = permute(cat(3, rho_nom, drho_nom, phi_nom), [1, 3, 2]);

mask_lin = mask(t_lin, x_lin, 1:nS);
mask_lin = permute(repmat(mask_lin, 1, 1, 3), [1, 3, 2]);
d_y_lin = permute(d_y_lin, [1, 3, 2]);
y_lin = ynom + d_y_lin;
y_lin(mask_lin) = NaN;

% Linearized Measurements
F_lin_meas = figure('DefaultAxesPosition', [0.05, 0.05, 0.95, 0.95]);
set(gcf, 'Position',  pos)
y_lin_ylabels = {'$\rho$ [km]','$\dot \rho$ [km/s]','$\phi$ [rad]'};
for j = 1:p
    ax = subplot(3, 1, j);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(y_lin_ylabels{j})
    if j == 3
        for i = 1:nS
            plot(t_lin, wrapToPi(y_lin(:, j, i)), '.', 'Color', colors(i, :)); 
        end 
    else
        for i = 1:nS
            plot(t_lin, y_lin(:, j, i), '.', 'Color', colors(i, :)); 
        end 
    end
    set(gca,'FontSize',12)
end
sgtitle('Linearized Measurements')
saveas(gcf, '5', 'epsc')

%% Part 2: Filtering

rng(100)

load orbitdeterm_finalproj_KFdata

% Define MC Simulation Parameters
n_trials = 100; % number of trials to generate
T = 1401; % a little over two orbital periods
k = (0:T-1)';
t = k*Dt;

% Generate Nominal Solution 

x_nom = xnom(t)';
rho_nom = rho(t, x_nom', 1:nS);
drho_nom = drho(t, x_nom', 1:nS);
phi_nom = phi(t, x_nom', 1:nS);
mask_nom = mask(t, x_nom', 1:nS);
mask_nom = permute(repmat(mask_nom, 1, 1, 3), [3, 2, 1]);
y_nom = permute(cat(3, rho_nom, drho_nom, phi_nom), [3, 2, 1]);
mu0_nom = [X0, dX0, Y0, dY0]';

% Generate Monte Carlo TMT Data

% P0 = diag([100, 1, 100, 1]); % good trial for EKF
P0 = diag([10, 0.01, 10, 0.000000001]);
% P0 = diag([0, 0, 0, 0]);
d_x0 = [10, 0.00, 10, 0.0]';
mu0 = [X0, dX0, Y0, dY0]'+d_x0;
u = [0, 0]';

x_TMT = zeros(n, T, n_trials);
y_TMT = zeros(p, nS, T, n_trials);
d_y_TMT = y_TMT;

fprintf('Running Monte Carlo TMT Simulations: ')
tic
for i = 1:n_trials
   
    S0 = mvnrnd(mu0, P0)';
%     S0 = mu0;
    w = mvnrnd([0; 0], Qtrue, T)';
%     w = mvnrnd([0; 0], zeros(2), T)';
    [~, x] = ode45(@(t, S) nl_sim(t, S, mu, u, w, Dt), t, S0, opts);
    
    rho_mc = rho(t, x, 1:nS);
    drho_mc = drho(t, x, 1:nS);
    phi_mc = phi(t, x, 1:nS);
    mask_mc = mask(t, x, 1:nS);
    mask_mc = permute(repmat(mask_mc, 1, 1, 3), [3, 2, 1]);
    
    % make time first dimension, measurement index second dimension, and
    % station index third dimension
    y_mc = permute(cat(3, rho_mc, drho_mc, phi_mc), [3, 2, 1]);
    y_mc(mask_mc) = NaN;
    v = permute(reshape(...
        mvnrnd(zeros(p, 1), Rtrue, T*nS), [T, nS, p]), [3, 2, 1]);
    d_y_mc = y_mc - y_nom;
    y_mc = y_mc + v;
    d_y_mc = d_y_mc + v;
    
    x_TMT(:, :, i) = x';
    y_TMT(:, :, :, i) = y_mc;
    d_y_TMT(:, :, :, i) = d_y_mc;
    
end
toc

% Generate Perturbation Data
x_nom_n = repmat(x_nom, [1, 1, n_trials]);
y_nom_n = repmat(y_nom, [1, 1, 1, n_trials]);

d_x_TMT = x_TMT - x_nom_n;

d_phi = d_y_TMT(3, :, :, :);
d_phi = wrapToPi(d_phi);
d_y_TMT(3, :, :, :) = d_phi;

% Generate NEES/NIS bounds
alpha = 0.05;
r1_NEES = chi2inv(alpha/2, n_trials*n)/n_trials;
r2_NEES = chi2inv(1-alpha/2, n_trials*n)/n_trials;
r1_NIS = chi2inv(alpha/2, n_trials*p)/n_trials;
r2_NIS = chi2inv(1-alpha/2, n_trials*p)/n_trials;


% Truth State
F_truth_states = figure('DefaultAxesPosition', [0.05, 0.05, 0.95, 0.95]);
set(gcf, 'Position',  pos)
x_lin_ylabels = {'$X$ [km]','$\dot X$ [km/s]',...
    '$Y$ [km]','$\dot Y$ [km/s]'};
for i = 1:n
    subplot(n, 1, i);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(x_nl_ylabels{i})
    plot(t, d_x_TMT(i, :, 1), 'LineWidth', 1.5)
    set(gca,'FontSize',12)
end
sgtitle('Truth Perturbations')
saveas(gcf, 'd_ts_LKF', 'epsc')

% Truth Measurements
F_truth_meas = figure('DefaultAxesPosition', [0.05, 0.05, 0.95, 0.95]);
set(gcf, 'Position',  pos)
y_lin_ylabels = {'$\rho$ [km]','$\dot \rho$ [km/s]','$\phi$ [rad]'};
for j = 1:p
    ax = subplot(3, 1, j);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(y_lin_ylabels{j})
    for i = 1:nS
        plot(t, reshape(d_y_TMT(j, i, :, 1), T, []), '.',...
            'Color', colors(i, :)); 
    end
    set(gca,'FontSize',12)
end
sgtitle('Truth Measurement Perturbations')
saveas(gcf, 'd_tm_LKF', 'epsc')

trial_plot = 2;

%% Linearized Kalman Filter

Qmag = permute(logspace(-3, 1, T), [1, 3, 2]);
f = 1e3; 
Qframe = repmat(f*diag([1, 1e-3, 1, 1e-3]), [1, 1, T]);
Qguess = Qmag.*Qframe;

% Qguess = 1e-10*eye(2);

% P0 = 1e-4*eye(4);

x_LKF = zeros(n, T, n_trials);
P_LKF = zeros(n, n, T, n_trials);
y_LKF = zeros(p, nS, T, n_trials);
S_LKF = zeros(p*nS, p*nS, T, n_trials);

fprintf('Running LKF over MC Simulations: ')
tic
for i = 1:n_trials
    d_y_i = d_y_TMT(:, :, :, i);
%     d_y_i(:, :, 400:1200) = NaN;
    [d_x_t, P_t, d_y_m, S_t] = LKF(...
        n, d_y_i, Qguess, Rtrue,...
        mu0-mu0_nom, P0, Fk, Omegak, Hk);
    x_LKF(:, :, i) = d_x_t;
    P_LKF(:, :, :, i) = P_t;
    y_LKF(:, :, :, i) = d_y_m;
    S_LKF(:, :, :, i) = S_t;
end
toc

x_LKF = x_LKF + x_nom_n;
y_LKF = y_LKF + y_nom_n;

% phi_LKF = y_LKF(3, :, :, :);
% phi_LKF = wrapToPi(phi_LKF);
% y_LKF(3, :, :, :) = phi_LKF;

e_x_LKF = x_TMT - x_LKF;
e_y_LKF = y_TMT - y_LKF;

% Run NEES/NIS Tests

NEES_LKF = zeros(T, n_trials);
NIS_LKF = NaN * zeros(T, n_trials);

fprintf('Running NEES/NIS on LKF Data: ')
tic
for i = 1:n_trials
    e_x_dat = e_x_LKF(:, :, i);
    P_dat = P_LKF(:, :, :, i);
    e_y_dat = e_y_LKF(:, :, :, i);
    S_dat = S_LKF(:, :, :, i);
    for j = 1:T
        
        % Run NEES Calculation
        e_x = e_x_dat(:, j);
        eps_x = e_x' * (P_dat(:, :, j) \ e_x);
        NEES_LKF(j, i) = eps_x;
        
        % Run NIS Calculation
        S_i = S_dat(:, :, j);
        S_inds = find(~isnan(S_i(1, :)), 1, 'last');
        S_i = S_i(1:S_inds, 1:S_inds);
        e_y = e_y_dat(:, :, j);
        e_y = e_y(:);
        e_y = e_y(~isnan(e_y));
        [r, ~] = size(S_i);
        eps_y = NaN;
        if r == length(e_y)
            e_y_norm_one_meas = S_i \ e_y;
            if r > 3
                eps_y = e_y(1:3)'*e_y_norm_one_meas(1:3);
            else
                eps_y = e_y' * e_y_norm_one_meas; 
            end 
        end
        NIS_LKF(j, i) = eps_y;
        
    end
end
toc

figure(1)
hold on;
grid on;
xlabel('$t$ [s]')
ylabel('$\epsilon_x$')
title('NEES Test Results, LKF')
plot(t, sum(NEES_LKF, 2)/n_trials, '.', 'Color', colors(1, :));
plot([t(1), t(end)], [r1_NEES, r1_NEES], '--',...
    'LineWidth', 3, 'Color', colors(2, :));
plot([t(1), t(end)], [r2_NEES, r2_NEES], '--',...
    'LineWidth', 3, 'Color', colors(2, :), 'HandleVisibility', 'off');
axis([-Inf, Inf, 0, 10*r2_NEES]);
set(gca,'FontSize',14)
saveas(gcf, 'NEES_LKF', 'epsc')


figure(2)
hold on;
grid on;
xlabel('$t$ [s]')
ylabel('$\epsilon_y$')
title('NIS Test Results, LKF')
plot(t, sum(NIS_LKF, 2)/n_trials, '.', 'Color', colors(1, :));
plot([t(1), t(end)], [r1_NIS, r1_NIS], '--',...
    'LineWidth', 3, 'Color', colors(2, :));
plot([t(1), t(end)], [r2_NIS, r2_NIS], '--',...
    'LineWidth', 3, 'Color', colors(2, :), 'HandleVisibility', 'off');
axis([-Inf, Inf, 0, 1.5*r2_NIS]);
set(gca,'FontSize',14)
saveas(gcf, 'NIS_LKF_full', 'epsc')

% Sample Run State Estimate Errors
F3 = figure('DefaultAxesPosition', [0.05, 0.05, 0.95, 0.95]);
set(gcf, 'Position',  pos)
x_nl_ylabels = {'$X$ [km]','$\dot X$ [km/s]','$Y$ [km]','$\dot Y$ [km/s]'};
for i = 1:n
    subplot(n, 1, i);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(x_nl_ylabels{i})
    plot(t, e_x_LKF(i, :, trial_plot), 'LineWidth', 1.5)
    P = permute(P_LKF(:, :, :, trial_plot), [1, 3, 2]);
    plot(t, 2*sqrt(P(i, :, i)), '--', 'LineWidth', 1.5,...
        'Color', colors(2, :));
    plot(t, -2*sqrt(P(i, :, i)), '--', 'LineWidth', 1.5,...
        'Color', colors(2, :), 'HandleVisibility', 'off');
    set(gca,'FontSize',12)
end
sgtitle('LKF State Estimate Errors')

saveas(gcf, 'se_LKF_full', 'epsc')

axes_bounds = [-Inf, Inf, -Inf, Inf; -Inf, Inf,...
    -Inf, Inf; -Inf, Inf, -0.4, 0.4];

% Sample Run Measurements
F4 = figure('DefaultAxesPosition', [0.05, 0.05, 0.95, 0.95]);
set(gcf, 'Position',  pos)
y_lin_ylabels = {'$\rho$ [km]','$\dot \rho$ [km/s]','$\phi$ [rad]'};
for j = 1:p
    ax = subplot(3, 1, j);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(y_lin_ylabels{j})
    if j == 4
        for i = 1:nS
            plot(t, wrapToPi(reshape(y_TMT(j, i, :, trial_plot), T, [])-...
                reshape(y_LKF(j, i, :, trial_plot), T, [])), '.',...
                'Color', colors(i, :)); 
        end
    else
        for i = 1:nS
            plot(t, reshape(y_TMT(j, i, :, trial_plot), T, [])-...
                reshape(y_LKF(j, i, :, trial_plot), T, []), '.',...
                'Color', colors(i, :)); 
        end 
    end
    set(gca,'FontSize',12)
    axis(axes_bounds(j, :))
end
sgtitle('Linearized Measurement Innovations')
saveas(gcf, 'mi_LKF_full', 'epsc')

F5 = figure('DefaultAxesPosition', [0.05, 0.05, 0.95, 0.95]);
set(gcf, 'Position',  [100, 100, 1000, 800])
y_lin_ylabels = {'$\rho$ [km]','$\dot \rho$ [km/s]','$\phi$ [rad]'};
for j = 1:p
    ax = subplot(3, 1, j);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(y_lin_ylabels{j})
    for i = 1:nS
%         plot(t, reshape(y_TMT(j, i, :, 1), T, []), '.',...
%             'Color', colors(1, :)); 
        plot(t, reshape(y_LKF(j, i, :, trial_plot), T, []), '.',...
            'Color', colors(i, :)); 
    end
end
sgtitle('Linearized Measurements')

%% Extended Kalman Filter

diags = 1.2e-10;
offdiags = -0.1*1e-10;
Qguess = [diags, offdiags; offdiags, diags];

x_EKF = zeros(n, T, n_trials);
P_EKF = zeros(n, n, T, n_trials);
y_EKF = zeros(p, nS, T, n_trials);
S_EKF = zeros(p*nS, p*nS, T, n_trials);

Fx = @(x) eye(n)+Dt*A(atan2(x(3), x(1)), sqrt(x(1)^2+x(3)^2));
Hx = @(k, x, i) C(k*Dt, x, i);

fprintf('Running EKF over MC Simulations: ')
tic
for i = 1:n_trials
    [x_t, P_t, y_m, S_t] = EKF(...
        n, y_TMT(:, :, :, i), Qguess, 0.9*Rtrue, mu0, P0, Fx, Omegak,...
        Hx, Dt, rho, drho, phi, mu);
    x_EKF(:, :, i) = x_t;
    P_EKF(:, :, :, i) = P_t;
    y_EKF(:, :, :, i) = y_m;
    S_EKF(:, :, :, i) = S_t;
end
toc

e_x_EKF = x_TMT - x_EKF;
e_y_EKF = y_TMT - y_EKF;

NEES_EKF = zeros(T, n_trials);
NIS_EKF = NaN * zeros(T, n_trials);

fprintf('Running NEES/NIS over EKF Data: ')
tic
for i = 1:n_trials
    e_x_dat = e_x_EKF(:, :, i);
    P_dat = P_EKF(:, :, :, i);
    e_y_dat = e_y_EKF(:, :, :, i);
    S_dat = S_EKF(:, :, :, i);
    for j = 1:T
        
        % Run NEES Calculation
        e_x = e_x_dat(:, j);
        eps_x = e_x' * (P_dat(:, :, j) \ e_x);
        NEES_EKF(j, i) = eps_x;
        
        % Run NIS Calculation
        S_i = S_dat(:, :, j);
        S_inds = find(~isnan(S_i(1, :)), 1, 'last');
        S_i = S_i(1:S_inds, 1:S_inds);
        e_y = e_y_dat(:, :, j);
        e_y = e_y(:);
        e_y = e_y(~isnan(e_y));
        [r, ~] = size(S_i);
        eps_y = NaN;
        if r == length(e_y)
            e_y_norm_one_meas = S_i \ e_y;
            if r > 3
                eps_y = e_y(1:3)'*e_y_norm_one_meas(1:3);
            else
                eps_y = e_y' * e_y_norm_one_meas; 
            end 
        end
        NIS_EKF(j, i) = eps_y;
        
    end
end
toc

%%
figure(5)
hold on;
grid on;
xlabel('$t$ [s]')
ylabel('$\epsilon_x$')
title('NEES Test Results, EKF')
plot(t, sum(NEES_EKF, 2)/n_trials, '.', 'Color', colors(1, :));
plot([t(1), t(end)], [r1_NEES, r1_NEES], '--',...
    'LineWidth', 3, 'Color', colors(2, :));
plot([t(1), t(end)], [r2_NEES, r2_NEES], '--',...
    'LineWidth', 3, 'Color', colors(2, :), 'HandleVisibility', 'off');
axis([-Inf, Inf, 0, 1.5*r2_NEES]);
saveas(gcf, 'NEES_EKF', 'epsc')
set(gca,'FontSize',12)

figure(6)
hold on;
grid on;
xlabel('$t$ [s]')
ylabel('$\epsilon_y$')
title('NIS Test Results, EKF')
plot(t, sum(NIS_EKF, 2)/n_trials, '.', 'Color', colors(1, :));
plot([t(1), t(end)], [r1_NIS, r1_NIS], '--',...
    'LineWidth', 3, 'Color', colors(2, :));
plot([t(1), t(end)], [r2_NIS, r2_NIS], '--',...
    'LineWidth', 3, 'Color', colors(2, :), 'HandleVisibility', 'off');
axis([-Inf, Inf, 0, 1.5*r2_NIS]);
saveas(gcf, 'NIS_EKF', 'epsc')
set(gca,'FontSize',12)


axes_bounds = [-Inf, Inf, -0.15, 0.15
    -Inf, Inf, -1.5e-3, 1.5e-3
    -Inf, Inf, -0.15, 0.15
    -Inf, Inf, -1.5e-3, 1.5e-3];

% Sample Run State Estimate Errors
F7 = figure('DefaultAxesPosition', [0.05, 0.05, 0.95, 0.95]);
set(gcf, 'Position',  pos)
x_nl_ylabels = {'$X$ [km]','$\dot X$ [km/s]','$Y$ [km]','$\dot Y$ [km/s]'};
for i = 1:n
    subplot(n, 1, i);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(x_nl_ylabels{i})
    plot(t, e_x_EKF(i, :, trial_plot), 'LineWidth', 1.5)
    P = permute(P_EKF(:, :, :, trial_plot), [1, 3, 2]);
    plot(t, 2*sqrt(P(i, :, i)), '--', 'LineWidth', 1.5,...
        'Color', colors(1, :));
    plot(t, -2*sqrt(P(i, :, i)), '--', 'LineWidth', 1.5,...
        'Color', colors(1, :), 'HandleVisibility', 'off');
    axis(axes_bounds(i, :))
    set(gca,'FontSize',12)
end
sgtitle('EKF State Estimate Errors')
saveas(gcf, 'se_EKF', 'epsc')

% Sample Run Measurements

axes_bounds = [-Inf, Inf, -0.5, 0.5
    -Inf, Inf, -Inf, Inf
    -Inf, Inf, -Inf, Inf];

F8 = figure('DefaultAxesPosition', [0.05, 0.05, 0.95, 0.95]);
set(gcf, 'Position',  pos)
y_lin_ylabels = {'$\rho$ [km]','$\dot \rho$ [km/s]','$\phi$ [rad]'};
for j = 1:p
    ax = subplot(3, 1, j);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(y_lin_ylabels{j})
    if j == 3
        for i = 1:nS
            plot(t, wrapToPi(reshape(y_TMT(j, i, :, trial_plot), T, [])-...
                reshape(y_EKF(j, i, :, trial_plot), T, [])), '.',...
                'Color', colors(i, :)); 
        end
    else
        for i = 1:nS
            plot(t, reshape(y_TMT(j, i, :, trial_plot), T, [])-...
                reshape(y_EKF(j, i, :, trial_plot), T, []), '.',...
                'Color', colors(i, :)); 
        end 
    end
    axis(axes_bounds(j, :))
    set(gca,'FontSize',12)
end
sgtitle('EKF Measurement Innovations')
saveas(gcf, 'mi_EKF', 'epsc')

%% Test On 'Actual' Data

T = length(ydata);
t = (0:T-1)'*Dt;

% format ydata

y_test = NaN*zeros(p, nS, T);
for i = 1:T
    Yi = ydata{i};
    if ~isempty(Yi)
        if ~isnan(Yi(4))
            y_test(:, Yi(4), i) = Yi(1:3);
        end
    end
end

% LKF
Qmag = permute(logspace(-5, -5, T), [1, 3, 2]);
f = 1;
Qframe = repmat(f*[1 0; 0 1], [1, 1, T]);
Qguess = Qmag.*Qframe;
% Qguess = 1e-10 * [1 0; 0 1];
P0 = diag([10, 0.1, 10, 0.1]);
[d_x_LKF, P_LKF, d_y_LKF, S_LKF] = LKF(n, y_test-y_nom, Qguess, Rtrue,...
        mu0-mu0_nom, P0, Fk, Omegak, Hk);

x_LKF = d_x_LKF + x_nom;
y_LKF = d_y_LKF + y_nom;
    

% EKF
Qguess = 1e-10 * [1 0; 0 1];
[x_EKF, P_EKF, y_EKF, S_EKF] = EKF(n, y_test, Qguess, Rtrue,...
    mu0, P0, Fx, Omegak, Hx, Dt, rho, drho, phi, mu);


% Sample Run State Estimate Errors
F_s_LKF = figure('DefaultAxesPosition', [0.05, 0.05, 0.95, 0.95]);
set(gcf, 'Position',  pos)
x_nl_ylabels = {'$X$ [km]','$\dot X$ [km/s]','$Y$ [km]','$\dot Y$ [km/s]'};
for i = 1:n
    subplot(n, 1, i);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(x_nl_ylabels{i})
    plot(t, x_LKF(i, :), 'LineWidth', 1.5)
    set(gca,'FontSize',14)
end
sgtitle('LKF State Estimate')
saveas(gcf, 'F_LKF_s', 'epsc')

F_c_LKF = figure('DefaultAxesPosition', [0.05, 0.05, 0.95, 0.95]);
set(gcf, 'Position',  pos)
x_nl_ylabels = {'$X$ [km]','$\dot X$ [km/s]','$Y$ [km]','$\dot Y$ [km/s]'};
for i = 1:n
    subplot(n, 1, i);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(x_nl_ylabels{i})
    P = permute(P_LKF, [1, 3, 2]);
    plot(t, 2*sqrt(P(i, :, i)), 'LineWidth', 1.5,...
        'Color', colors(1, :));
    set(gca,'FontSize',14)
end
sgtitle('LKF $2\sigma$ Bounds')
saveas(gcf, 'F_LKF_c', 'epsc')




F_s_EKF = figure('DefaultAxesPosition', [0.05, 0.05, 0.95, 0.95]);
set(gcf, 'Position',  pos)
x_nl_ylabels = {'$X$ [km]','$\dot X$ [km/s]','$Y$ [km]','$\dot Y$ [km/s]'};
for i = 1:n
    subplot(n, 1, i);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(x_nl_ylabels{i})
    plot(t, x_EKF(i, :), 'LineWidth', 1.5)
    set(gca,'FontSize',14)
end
sgtitle('EKF State Estimate')
saveas(gcf, 'F_EKF_s', 'epsc')

axes_bounds = [-Inf, Inf, 0, 1
    -Inf, Inf, 0, 5e-3
    -Inf, Inf, 0, 1
    -Inf, Inf, 0, 5e-3];

F_c_EKF = figure('DefaultAxesPosition', [0.05, 0.05, 0.95, 0.95]);
set(gcf, 'Position',  pos)
x_nl_ylabels = {'$X$ [km]','$\dot X$ [km/s]','$Y$ [km]','$\dot Y$ [km/s]'};
for i = 1:n
    subplot(n, 1, i);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(x_nl_ylabels{i})
    P = permute(P_EKF, [1, 3, 2]);
    plot(t, 2*sqrt(P(i, :, i)), 'LineWidth', 1.5,...
        'Color', colors(1, :));
    axis(axes_bounds(i, :))
    set(gca,'FontSize',14)
end
sgtitle('EKF $2\sigma$ Bounds')
saveas(gcf, 'F_EKF_c', 'epsc')



% Sample Run Measurements
figure(15)
y_lin_ylabels = {'$\rho$ [km]','$\dot \rho$ [km/s]','$\phi$ [rad]'};
for j = 1:p
    ax = subplot(3, 1, j);
    hold on;
    grid on;
    xlabel('$t$ [s]')
    ylabel(y_lin_ylabels{j})
    for i = 1:nS
        plot(t, reshape(y_EKF(j, i, :), T, [])-reshape(y_LKF(j, i, :), T, []), '.',...
            'Color', colors(i, :)); 
%         plot(t, reshape(y_LKF(j, i, :), T, []), 'o',...
%             'Color', colors(i, :)); 
    end
end
sgtitle('Filtered Measurements')


%% Functions

% Nonlinear Simulation Function
function dS = nl_sim(t, S, mu, u, w, Dt)

    X = S(1);
    Y = S(3);
    dX = S(2);
    dY = S(4);
    r = sqrt(X^2+Y^2);
    i = round(t/Dt) + 1;
    ddX = -mu*X/r^3 + u(1) + w(1, i);
    ddY = -mu*Y/r^3 + u(2) + w(2, i);
    
    dS = [dX; ddX; dY; ddY];

end

function dS = nl_sim_nonoise(~, S, mu)

    X = S(1);
    Y = S(3);
    dX = S(2);
    dY = S(4);
    r = sqrt(X^2+Y^2);
    ddX = -mu*X/r^3;
    ddY = -mu*Y/r^3;
    
    dS = [dX; ddX; dY; ddY];

end

% absolute angular difference
function d = absangdiff(a, b)
    d = abs(b - a);
    d = abs(mod(d, 2*pi) - pi);
end

% Linearized Kalman Filter 
function [dx, P, dy, S] = LKF(...
    n, dy_meas, Q, Rtrue, mu0, P0, Fk, Omega, Hk)
    
    [p, nS, T] = size(dy_meas);
    
    dx = zeros(n, T);
    dy = NaN*zeros(p, nS, T);
    P = zeros(n, n, T);
    S = NaN*zeros(p*nS, p*nS, T);
    
    k = (0:T-1)';
    
    mu_prev = mu0;
    P_prev = P0;
    
    dx(:, 1) = mu_prev;
    P(:, :, 1) = P_prev;
    
    dyframe = dy_meas(:, :, 1);
    yind = find(~isnan(dyframe(1, :)));
    n_meas = length(yind);
    if n_meas > 0
        H = zeros(p*n_meas, n);
        for j = 1:n_meas
            H((p*j-(p-1)):(p*j), :) = Hk(k(1), yind(j)); 
        end
        R = kron(eye(n_meas), Rtrue);
        S(1:p*n_meas, 1:p*n_meas, 1) = H*P0*H'+R;
    end
    
    [rQ,~,nQ] = size(Q);
    
    for i = 2:T
        F = Fk(i-2);
        
        dxi = F*mu_prev;
        Pi = F*P_prev*F';
        if nQ == T
            if rQ == 2
                Pi = Pi + Omega*Q(:, :, i)*Omega';
            else
                Pi = Pi + Q(:, :, i);
            end
        else
            if rQ == 2
                Pi = Pi + Omega*Q(:, :, 1)*Omega';
            else
                Pi = Pi + Q(:, :, 1)';
            end
        end
%         Pi = F*P_prev*F'+Omega*Q*Omega';
        
%         dxi = [0; 0; 0; 0];
        
        dyframe = dy_meas(:, :, i);
        yind = find(~isnan(dyframe(1, :)));
        n_meas = length(yind);
        if n_meas > 0
            H = zeros(p*n_meas, n);
            for j = 1:n_meas
                H((p*j-(p-1)):(p*j), :) = Hk(k(i), yind(j)); 
            end
            
            dyi = dyframe(:, yind);
            dyi = dyi(:);
            
            dym = H*dxi;
            ddy = dyi-dym;
%             ddy = zeros(size(ddy));
%             dym = reshape(dym, [p, n_meas]);
%             ddy = reshape(ddy, [p, n_meas]);
%             dphi = ddy(3, :);
%             inds = find(abs(dphi) < pi);
%             dym = dym(:, inds);
%             ddy = ddy(:, inds);
%             n_meas = length(inds);
%             yind = yind(inds);
%             ddy = ddy(:);
%             dym = dym(:);
%             
%             H = zeros(p*n_meas, n);
%             for j = 1:n_meas
%                 H((p*j-(p-1)):(p*j), :) = Hk(k(i), yind(j)); 
%             end
            
            
            R = kron(eye(n_meas), Rtrue);
            Si = H*Pi*H'+R;
            S(1:p*n_meas, 1:p*n_meas, i) = Si;
            K = Pi*H'*(Si\eye(p*n_meas));
            Pi = (eye(n)-K*H)*Pi;
            
            dxi = dxi+K*ddy;
           
            dy(:, yind, i) = reshape(dym, [p, n_meas]);
            
        end
        
        mu_prev = dxi;
        P_prev = Pi;
        dx(:, i) = mu_prev;
        P(:, :, i) = P_prev;
    
    end
    
end

% Extended Kalman Filter
function [x_t, P_t, y_m, S_t] = EKF(...
    n, y, Q, R_meas, mu0, P0, Fx, Omega, Hx, Dt, rho, drho, phi, mu)

    [p, nS, T] = size(y);
    
    x_t = zeros(n, T);
    P_t = zeros(n, n, T);
    
    y_m = NaN * zeros(p, nS, T);
    
    mu_prev = mu0;
    P_prev = P0;
    
    x_t(:, 1) = mu_prev;
    P_t(:, :, 1) = P_prev;
    
    k = (0:(T-1))';
    
    tspan = [0, 1, Dt];
    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    
    S_t = NaN * zeros(p*nS, p*nS, T);
    
    % Find first S
    yinds = find(~isnan(y(1, :, 1)));
    n_meas = length(yinds);
    if n_meas > 0 
        H = zeros(n_meas*p, n);
        for j = 1:n_meas
            H((p*j-2):(p*j),:) = Hx(k(1), mu_prev, yinds(j));
        end
        R = kron(eye(n_meas), R_meas);
        S_t(1:n_meas*p, 1:n_meas*p, 1) = H*P0*H'+R;
    end

    for i = 2:T

        [~, x] = ode45(@(t, S) nl_sim_nonoise(t, S, mu),...
            tspan, mu_prev, opts);
        
        x = x(end, :);
        
        F = Fx(mu_prev);
        P = F*P_prev*F' + Omega*Q*Omega';
        
        yframe = y(:, :, i);
        yinds = find(~isnan(yframe(1, :)));
        n_meas = length(yinds);
        mu_prev = x';
        P_prev = P;
        if n_meas > 0
            R = kron(eye(n_meas), R_meas);
            H = zeros(p*n_meas, n);
            for j = 1:n_meas
                H((p*j-2):(p*j), :) = Hx(k(i), x, yinds(j)); 
            end
            
            yi = yframe(:, yinds);
            yi = yi(:);
            
            rho_mi = rho(k(i)*Dt, x, yinds);
            drho_mi = drho(k(i)*Dt, x, yinds);
            phi_mi = phi(k(i)*Dt, x, yinds);
            ymi2 = [rho_mi; drho_mi; phi_mi];
            ymi = ymi2(:);
            
            e_yi = yi - ymi;
            
            S = H*P*H'+R;
            K = P*H'*(S\eye(p*n_meas));
            mu_prev = x' + K*e_yi;
            P_prev = (eye(n)-K*H)*P;
            
            S_t(1:p*n_meas, 1:p*n_meas, i) = S;
            
            y_m(:, yinds, i) = ymi2;
        end
        
        x_t(:, i) = mu_prev;
        P_t(:, :, i) = P_prev;
        
    end

end

