function [t,e,x49]=FFESMVA(t,is)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOTION INPUT - FLEXIBLE ANALYSIS OF BT SOURCE IN RING APPLICATOR 
% SUMMARY - This file determines the input position, velocity and acceleration 
% INPUT   - t:   time
%         - 'Simulation.bin': final time and insertion velocity
% OUTPUT  - x145:  position, velocity and acceleration on input node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = t;
e = [];

% Open simulation parameters
fid   = fopen('Simulation.bin', 'r');
tmat  = fread(fid,'double'); fclose(fid);

% Insert parameters
tlim     = tmat(1,:);     % final time
vel      = tmat(2,:);     % insertion velocity

if t<=tlim
    x49 = -0.12 + vel*t;
    v49 =  vel;
    a49 =  0;
else
    x49 =  0;
    v49 =  0;
    a49 =  0;
end

x49   = [49 1 x49 v49 a49];

end