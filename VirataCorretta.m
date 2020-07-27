clc; clear all; close all;
%% Declaring global variables
global ...
g ... % gravity acceleration
rho_0 V_0 q_0 gamma_0... % air density at z_0
delta_T_0 delta_s_0 ...
vTime_bp vCL_bp ...
myAC sound_speed_0 ...% the aircraft object, populated outside this func.
M0All M01All M02All ...
M03All M04All M05All ...
alpha_0

M0All=datcomimport([pwd,'\DATCOM\FokkerF50(Mach0).out'],true,0);
M01All=datcomimport([pwd,'\DATCOM\FokkerF50(Mach01).out'],true,0);
M02All=datcomimport([pwd,'\DATCOM\FokkerF50(Mach02).out'],true,0);
M03All=datcomimport([pwd,'\DATCOM\FokkerF50(Mach03).out'],true,0);
M04All=datcomimport([pwd,'\DATCOM\FokkerF50(Mach04).out'],true,0);
M05All=datcomimport([pwd,'\DATCOM\FokkerF50(Mach05).out'],true,0);

aircraftDataFileName = 'DSV_Aircraft_data.txt';
myAC = DSVAircraft(aircraftDataFileName);

%% Initial conditions
z_0     = -7620.;                  % a.s.l. altitude (m)
V_0     = 145;                    % flight speed
q_0     = 0.;                      % pitching angular speed (rad/s(
gamma_0 = convang(0,'deg','rad'); % path climb angle (rad)

%% Using Matlab built-in ISA model for density
[air_Temp_0, sound_speed_0, air_pressure_0, rho_0] = atmosisa(-z_0);

%% Gravity acceleration
g = 9.81; % (m/s^2)

%% Design vector for trim
% alpha   = xi(1);
% delta_e = xi(2);
% delta_s = xi(3);
% delta_T = xi(3);

%% initial guess for the design vector
xi0 = [  ...
    0;  ... % alpha_0
    0;  ... % delta_e_0
    0;  ... % delta_s_0
    0.5 ... % delta_T_0
    ];

%% Cost function minimization

% Aeq, in Aeq*x=beq linear constraint
Aeq       = zeros(4,4);
beq       = zeros(4,1);
%Aeq(2,2)  = 1;
Aeq(3,3)  = 1; % pick delta_s
delta_s_0 = convang(-1,'deg','rad');
beq(3,1)  = delta_s_0; % keep delta_s fix
%beq(2,1)  = convang(-0.55,'deg','rad');

% bounds
lb =[convang(-15,'deg','rad'), ... % minimum alpha
    convang(-20,'deg','rad'), ... % minimum elevator deflection
    convang(-5,'deg','rad'),  ... % minimum stabilizer incidence
    0.2                       ... % minimum thrust fraction
    ];
ub =[convang( 15,'deg','rad'), ... % maximum alpha
    convang( 13,'deg','rad'), ... % maximum elevator deflection
    convang(  2,'deg','rad'), ... % maximum stabilizer incidence
    1.0                       ... % maximum thrust fraction
    ];

options = optimset( ...
    'tolfun',1e-9, ...               % threshold
    'Algorithm','interior-point' ... % algor. type
    );
[xi,fval] = ...
    fmincon(@costLongEquilibriumStaticStickFixed, ...
    xi0,   ...
    [],    ... % A, A*x<=b
    [],    ... % b
    Aeq,   ... % Aeq, Aeq*x=beq
    beq,   ... % beq
    lb,ub, ...
    @myNonLinearConstraint, ...
    options);

alpha_0       = xi(1);
alpha_0_deg   = convang(alpha_0,'rad','deg');
theta_0       = gamma_0 + alpha_0 + myAC.mu_x;
theta_0_deg   = convang(theta_0,'rad','deg');
delta_e_0     = xi(2);
delta_e_0_deg = convang(delta_e_0,'rad','deg');
delta_s_0     = xi(3);
delta_s_0_deg = convang(delta_s_0,'rad','deg');
delta_T_0     = xi(4);

%% Solve the trim problem first, wings-level and zero R/C
% load trim results
% Data read from file:
% V_0, alpha_0_deg, alpha_0, delta_e_0_deg, delta_e_0
% delta_s_0_deg, delta_s_0, delta_T_0
% initial heading
psiGT_0 = convang(0.0,'deg','rad');
% initial load factor
fza_0 = 1.0;
% initial CL
[CL_eq,~,~,~,~] = AeroCoeff(alpha_0_deg, delta_e_0_deg, V_0/sound_speed_0, 0, 0,V_0);
%% Setting up the manoeuver
t_fin = 30.0; % Simulation final time
% Breakpoints in time law of CL
vTime_bp = [0, t_fin*(1/10), t_fin*(1/4), t_fin*(2/4), t_fin*(2/3), t_fin];
vCL_bp = [CL_eq, CL_eq+0.01, CL_eq+0.10, CL_eq+0.38, CL_eq+0.395, CL_eq+0.4];
%% Solve the problem (17.61)
% Mass matrix
M = diag([1,1,0,0]');
%-- initial state vector
x0 = [ ...
V_0; ... % x1 <-- V
psiGT_0; ... % x2 <-- psiGT
fza_0; ... % x3 <-- fza
delta_e_0]; ... % x4 <-- delta_e
options_ODE = odeset('AbsTol',1e-7,'Mass',M);
[vTime, vX] = ode15s(@correctedTurnCLAssigned, [0 t_fin], x0, options_ODE);
%% auxiliary variables
vY(:,1) = ... % bank angle
acos(1./vX(:,3));
vY(:,2) = ... % turn radius
(vX(:,1).^2)./(g.*sqrt(vX(:,3).^2 - 1));
options= optimset('fzero');
q=g./vX(:,1).*(vX(:,3)-1./vX(:,3));
CL=interp1(vTime_bp,vCL_bp,vTime,'pchip');
for i=1:length(vTime)
    vY(i,3) = fzero(@ZeroAlpha, alpha_0, options,CL(i),vX(i,4),vX(i,1),q(i));
end
for i=1:length(vTime)
    [~,CD(i),~,~,~] = AeroCoeff(vY(i,3)*180/pi,vX(i,4)*180/pi,vX(i,1)/sound_speed_0,0,q(i)*180/pi,vX(i,1));
end

vY(:,4) = CD;
%% trajectory - Solution of navigation equations
% RHS of navigation equations
dPosEdt = @(t,Pos) ...
[ ...
interp1(vTime,vX(:,1).*cos(vX(:,2)),t); ... 
interp1(vTime,vX(:,1).*sin(vX(:,2)),t); ... 
0];
options = odeset( ...
'RelTol', 1e-3, ...
'AbsTol', 1e-3*ones(3,1) ...
);
PosE0 = [0;0;z_0];
[vTime2, vPosE] = ode45(dPosEdt, vTime, PosE0, options);

%% Euler angles
% Rotation of angle a about 3rd axis
Rot3 = @(a) ...
[ cos(a),sin(a),0; ...
-sin(a),cos(a),0; ...
0, 0,1];
% Rotation of angle a about 2nd axis
Rot2 = @(a) ...
[ cos(a), 0, -sin(a); ...
0, 1, 0; ...
sin(a), 0, cos(a)];
% Rotation of angle a about 1st axis
Rot1 = @(a) ...
[ 1, 0, 0; ...
0, cos(a), sin(a); ...
0, -sin(a), cos(a)];
% Pre-assign Euler angles arrays
vPsi = zeros(length(vTime),1);
vTheta = zeros(length(vTime),1);
vPhi = zeros(length(vTime),1);
% Apply transformations (17.2) from known . gt;
for k=1:length(vTime)
    psigt = vX(k,2); nu = vY(k,1); alpha = vY(k,3);
    % Rotation sequence: 3-1-2
    Tbe = Rot2(alpha)*Rot1(nu)*Rot3(psigt);
    [psi,theta,phi] = dcm2angle(Tbe);
    vPsi(k) = psi; vTheta(k) = theta; vPhi(k) = phi;
end
% Time sequence of aircraft orientation quaternion
vQuat = angle2quat(vPsi,vTheta,vPhi);


%% Plot section
h_fig=figure(1);
vXe=vPosE(:,1);
vYe=vPosE(:,2);
vZe=vPosE(:,3)+5500;
scale_factor=0.002;
step=[1,30,40,48,53];
theView=[1,1,0.5];
plotTrajectoryAndBody2(h_fig,vXe,vYe,vZe,vQuat,scale_factor,step,theView)

T_fitto=linspace(0,t_fin,1000);

figure (2);
CL_fitto=interp1(vTime,CL,T_fitto,'pchip');
plot(T_fitto,CL_fitto,'-k');
xlabel('$t (s)$','Interpreter','Latex','FontSize',12);
ylabel('$C_{L}$','Interpreter','Latex','FontSize',12);
set(get(gca,'YLabel'),'Rotation',0,'Position',[-1.3 0.67]);
grid on;
grid minor;

figure(3)
subplot(2,2,1)
V_fitto=interp1(vTime,vX(:,1),T_fitto,'pchip');
plot(T_fitto,V_fitto,'-k');
xlabel('$t (s)$','Interpreter','Latex','FontSize',12);
ylabel('$V (m/s)$','Interpreter','Latex','FontSize',12);
set(get(gca,'YLabel'),'Rotation',0,'Position',[-2.5 131]);
grid on;
grid minor;

subplot(2,2,2)
Psi_fitto=interp1(vTime,vX(:,2)*180/pi,T_fitto,'pchip');
plot(T_fitto,Psi_fitto,'-k');
xlabel('$t (s)$','Interpreter','Latex','FontSize',12);
ylabel('$\psi (deg)$','Interpreter','Latex','FontSize',12);
set(get(gca,'YLabel'),'Rotation',0,'Position',[-2.5 65]);
grid on;
grid minor;

subplot(2,2,3)
fza_fitto=interp1(vTime,vX(:,3),T_fitto,'pchip');
plot(T_fitto,fza_fitto,'-k');
xlabel('$t (s)$','Interpreter','Latex','FontSize',12);
ylabel('$f_{zA}$','Interpreter','Latex','FontSize',12);
set(get(gca,'YLabel'),'Rotation',0,'Position',[-2.5 1.42]);
grid on;
grid minor;

subplot(2,2,4)
deltae_fitto=interp1(vTime,vX(:,4),T_fitto,'pchip');
plot(T_fitto,deltae_fitto*180/pi,'-k');
xlabel('$t (s)$','Interpreter','Latex','FontSize',12);
ylabel('$\delta_{e} (deg)$','Interpreter','Latex','FontSize',12);
set(get(gca,'YLabel'),'Rotation',0,'Position',[-2.5 -3.8]);
grid on;
grid minor;

figure(4)
subplot(2,2,1)
bank_fitto=interp1(vTime,vY(:,1),T_fitto,'pchip');
plot(T_fitto,bank_fitto*180/pi,'-k');
xlabel('$t (s)$','Interpreter','Latex','FontSize',12);
ylabel('$\nu (deg)$','Interpreter','Latex','FontSize',12);
set(get(gca,'YLabel'),'Rotation',0,'Position',[-2.5 36]);
grid on;
grid minor;

subplot(2,2,2)
radius_fitto=interp1(vTime,vY(:,2),T_fitto,'pchip');
plot(T_fitto,radius_fitto,'-k');
xlabel('$t (s)$','Interpreter','Latex','FontSize',12);
ylabel('$R (m)$','Interpreter','Latex','FontSize',12);
set(get(gca,'YLabel'),'Rotation',0,'Position',[-2.5 3300]);
grid on;
grid minor;

subplot(2,2,3)
alpha_fitto=interp1(vTime,vY(:,3),T_fitto,'pchip');
plot(T_fitto,alpha_fitto*180/pi,'-k');
xlabel('$t (s)$','Interpreter','Latex','FontSize',12);
ylabel('$\alpha (deg)$','Interpreter','Latex','FontSize',12);
set(get(gca,'YLabel'),'Rotation',0,'Position',[-2.5 6.3]);
grid on;
grid minor;


subplot(2,2,4)
CD_fitto=interp1(vTime,vY(:,4),T_fitto,'pchip');
plot(T_fitto,CD_fitto,'-k');
xlabel('$t (s)$','Interpreter','Latex','FontSize',12);
ylabel('$C_{D}$','Interpreter','Latex','FontSize',12);
set(get(gca,'YLabel'),'Rotation',0,'Position',[-2.5 0.088]);
grid on;
grid minor;

