function dx = correctedTurnCLAssigned(t,x)
% pass this to ode15s
%% Declaring global variables
global ...
g ... % gravity acceleration
rho_0 ... % air density at z_0
delta_T_0 delta_s_0 ...
vTime_bp vCL_bp ...
myAC alpha_0 sound_speed_0 % the aircraft object, populated outside this func.
%% Give the state vector components proper names
V = x(1);
psiGT = x(2);
fza = x(3);
delta_e = x(4);
% give dimensions to the return vector
dx = zeros(4,1);
%% Inertia settings 
Ixz= 3000;
%% Right-hand-sides of equations of motion
rho = rho_0;
delta_T = delta_T_0;
CL = interp1(vTime_bp,vCL_bp,t,'pchip');
ToW = delta_T*myAC.T/myAC.W;
WoS = myAC.W/myAC.S;
options= optimset('fzero');
q=g/V*(fza-1/fza);
a_ = fzero(@ZeroAlpha, alpha_0, options,CL,delta_e,V,q);
b_ = (rho*V^2)/(2*WoS);
[~,CD,~,~,~] = AeroCoeff(a_*180/pi,delta_e*180/pi,V/sound_speed_0,0,q*180/pi,V);
dx(1,1) = g*( ... 
ToW*( cos(myAC.mu_T) - a_*sin(myAC.mu_T) ) ...
- b_*(CD) ...
);
dx(2,1) = (g/V)*sqrt(fza^2 - 1); 
dx(3,1) = -fza ... 
+ ToW*( a_*cos(myAC.mu_T) + sin(myAC.mu_T) ) ...
+ b_*CL; 
M_ob=-Ixz*((dx(2,1))^2)/fza^2;
CM_ob=M_ob/(0.5*rho_0*V^2*myAC.S*myAC.mac);
[~,~,CM,~,~] = AeroCoeff(a_*180/pi,delta_e*180/pi,V/sound_speed_0,0,q*180/pi,V);
dx(4,1)=CM_ob-(CM+myAC.Cm_delta_s*delta_s_0);
end