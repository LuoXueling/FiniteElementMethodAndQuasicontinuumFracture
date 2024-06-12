function [disp_n1,velo_n1,acce_n1,totalDisp,loadTime,P]=initialize(elemCfg)
% Initialize field variables.
disp_n1=sparse(elemCfg.meshDispDOF,1);  % Displacement
velo_n1=sparse(elemCfg.meshDispDOF,1);  % Velocity
acce_n1=sparse(elemCfg.meshDispDOF,1);  % Acceleration
totalDisp=0;
loadTime=0;
P=0;