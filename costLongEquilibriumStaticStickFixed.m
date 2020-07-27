function f = costLongEquilibriumStaticStickFixed(x)

% Cost function for trimming the 3-DoF aircraft motion.
% The right-hand-sides of equations of motion are calculated
% then squared and summed up to form a non-negative function.
% A minimum condition for the cost function means that a trim 
% condition has been reached.
% Here x is the 'design variable' in the design space

%% Declaring global variables
global ...
    g ...                   % gravity acceleration
    z_0 V_0 q_0 gamma_0... % initial conditions
    rho_0 ...               % air density at z_0
    myAC sound_speed_0      % the aircraft object, populated outside this func.

%% Aircraft relative density
mu_rel = (myAC.W/g)/(rho_0*myAC.S*myAC.b);

%% Give the design vector components proper names
alpha   = x(1);
delta_e = x(2);
delta_s = x(3);
delta_T = x(4);
[CL,CD,CM,CL_q,~] = AeroCoeff(alpha*180/pi,delta_e*180/pi,V_0/sound_speed_0,0,q_0*180/pi,V_0);
%% Equation of motion right-hand-sides, for steady flight

F1 = (delta_T*myAC.T/myAC.W)*cos(alpha - myAC.mu_x + myAC.mu_T) ...
       -sin(gamma_0) ...
       -((rho_0*V_0^2)/(2*(myAC.W/myAC.S))) ...
         *(CD);

F2 = ( 1. - CL_q*180/pi*(myAC.mac/myAC.b)/(4*mu_rel) )*q_0 ...
       -(g/V_0)*(delta_T*myAC.T/myAC.W)*sin(alpha - myAC.mu_x + myAC.mu_T) ...
       +(g/V_0)*cos(gamma_0) ...
       -((rho_0*V_0*g)/(2*(myAC.W/myAC.S))) ...
         *(CL+myAC.CL_delta_s*delta_s);
 
F3 = CM+myAC.Cm_T_0+myAC.Cm_delta_s*delta_s;

%% Build the cost function
f = F1*F1 + F2*F2 + F3*F3;
end
