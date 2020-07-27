function CL_Diff = ZeroAlpha(alpha,CL_ASS,delta_e,V,q)
global sound_speed_0 delta_s_0 myAC
[CL_AC,~,~,~,~] = AeroCoeff(alpha*180/pi,delta_e*180/pi,V/sound_speed_0,0,q*180/pi,V);
CL_Diff=CL_ASS-(CL_AC+myAC.CL_delta_s*delta_s_0);
end