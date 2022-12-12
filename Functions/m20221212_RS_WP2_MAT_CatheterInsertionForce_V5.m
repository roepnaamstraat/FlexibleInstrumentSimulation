%% DYNAMIC MULTIBODY MODELS OF CATHETER INSERTIONS IN S-SHAPED CHANNELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Robin Straathof, 2022
%  r.straathof-1@tudelft.nl

%  SUMMARY - This program runs a simplified rigid/flexible multibody analysis of a BT 
%  catheter and obturator in S-shaped channels. It first constructs the
%  centreline of the channel based on a set of design parameters. All data
%  is saved in a struct: h.  
%
%  The software then exports data to bin files such that these can be read
%  in the user-defined routines and runs SPACAR. The output of the dynamic
%  analysis is saved in a struct: sim.
%
%  As an example the insertion force versus insertion depth is plotted. 

%  REQUIREMENTS - MATLAB R2021b
%  INPUTS - In same folder
%   - RIG_CAT_S.dat:  containing SPACAR inputs for rigid analysis    
%   - FLEX_CAT_S.dat: containing SPACAR inputs for flexible analysis
%   - FFESFVA.m:      containing user-defined force routine for flexible analysis
%   - FFESMVA.m:      containing user-defined motion routine for flexible analysis
%   - FRESFVA.m:      containing user-defined force routine for rigid analysis
%   - FRESMVA.m:      containing user-defined motion routine for rigid analysis
%   - distance2curve.m:   computing distance to curve by D'Errico(2012)
%   - LiePose.m       computing centreline curve using Lie group / algebra theory

%   ... In other folder:
%   - SPACAR 2017 needs to be installed and linked to MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clc 
clear all

% Add path to SPACAR:
addpath('...\Toolbox\spacar')   % Define path to SPACAR here

%% Specification of model parameters and generation of centreline curve
%  General simulation parameters
stiffn   = 1000000;  % [N/m]:   wall stiffness
c_w      = 10;       % [Ns/m]:  wall damping coefficient
mu_s_SLS = 0.094;    % []:      static friction coefficient SLS
mu_k_SLS = 0.074;    % []:      kinetic friction coefficient SLS
v_brk    = 0.001;    % [m/s]:   breakaway velocity

ro       = 0.00099;  % [m]:     outer radius catheter 

vel_ins  = 0.005;    % [m/s]:   insertion velocity
vel_sim  = 0.1;      % [m/s]:   simulation velocity

%  Simulation specific parameters
N        =  20;      % []:      discretisation of points in channel
x0       = -0.19;    % [m]:     x-coordinate of start of channel
d_1      =  0.01;    % [m]:     length first straight section
d_3      =  0.03;    % [m]:     length third straight section
ra       =  0.0012;  % [m]:     inner transition region
rb       =  0.0013;  % [m]:     outer transition region
tau      =  0;       % []:      geometric torsion
r        =  0.05;    % [m]:     radius of curvature
w        =  0.03;    % [m]:     planar width of channel
ht       =  0.12;    % [m]:     planar height of channel
depth    =  0.12;    % [m]:     planned insertion depth

% Generate S-shaped channel using Lie algebra
syms d_2_sym t_1_sym
grey     = [0.5,0.5,0.5];

[d_2_sol,t_1_sol] = solve(d_2_sym*sin(t_1_sym) + 2*r*(1-cos(t_1_sym)) == w, d_1 + d_3 + d_2_sym*cos(t_1_sym) + 2*r*sin(t_1_sym) == ht);
[d_2,in] = max(d_2_sol);
d_2      = double(d_2);
t_1      = double(t_1_sol(in));

% Generate S-shaped channel using Lie algebra
X_t   = [1,0, 0,x0;       % original pose of channel
         0,1, 0, 0;
         0,0, 1, 0;
         0,0, 0, 1];

u_t   = [(-x0+d_1)/N   , 1*10^(-20) , 0               ; % control inputs
    t_1*r/N       , 1/r        , 0               ;
    d_2/N         , 1*10^(-20) , (pi-tau)/d_2    ;
    t_1*r/N       , 1/r        , 0               ;
    d_3/N         , 1*10^(-20) , 0              ];

clear p_t
p_t     = [x0;0;0];
for j = 1:(5*N)                                      % execute controls
    X_t = LiePose(X_t,u_t(floor((j-1)/N)+1,:));
    p_t(:,j+1) = X_t(1:3,4);
end
x_c     = p_t';

% Save to structure
h.chadata.x_c   = x_c;                   % centreline data
h.chadata.ra    = ra;                    % transition region
h.chadata.rb    = rb;                    % lumen radius
h.chadata.tau   = tau;                   % geometric torsion
h.chadata.r     = r;                     % channel curvature

%% Run simulation
% Run simulation for all channels in template
sim_mod = 'flex';        % 'rig': Rigid analysis, 'flex': Flexible analysis

% Input for SPACAR analysis
x_c = h.chadata.x_c;                     % centreline data
ra  = h.chadata.ra;                      % transition region
rb  = h.chadata.rb;                      % lumen radius
tau = h.chadata.tau;                     % geometric torsion
r   = h.chadata.r;                       % channel curvature
tlim  = depth/vel_sim;                   % time till end of insertion

% Write data to binary files
fileID1 = fopen('CentreLine.bin','w');
fwrite(fileID1, x_c, 'double');fclose(fileID1);

fileID2 = fopen('Channel.bin','w');
fwrite(fileID2, [rb,ra,0;stiffn,0,0;c_w,0,0;mu_s_SLS,mu_k_SLS,v_brk;...
    vel_sim/vel_ins,0,0], 'double');fclose(fileID2);

fileID3 = fopen('Instrument.bin','w');
fwrite(fileID3,ro, 'double');fclose(fileID3);

fileID4 = fopen('Simulation.bin','w');
fwrite(fileID4, [tlim;vel_sim], 'double');fclose(fileID4);

% Limit end time through .dat file
% Run rigid analysis
if matches('rig',sim_mod)
    A = regexp( fileread('RIG_CAT_S_3D_V2.dat'), '\n', 'split');
    A{229} = sprintf(append('TIMESTEP ',num2str(tlim),' 8000          %% time and number of steps'));
    fid = fopen('RIG_CAT_S_3D_V3.dat', 'w');
    fprintf(fid, '%s\n', A{:});
    fclose(fid);

    tstart = tic;
    spacar(1,'RIG_CAT_S_3D_V3');                       % forward dynamic analysis
    tend   = toc(tstart);
    spavisual('RIG_CAT_S_3D_V3');                      % visualise outcome

    % SPACAR output
    irig         = 1:4:49;                               % nodes
    noderig      = lnp(irig,1:3);                        % nodal matrix
    tval         = min(abs(time-tlim));
    t_inx        = find(abs(time-tlim)==tval);

    sim.rig.x_c   = x_c;
    sim.rig.ra    = ra;
    sim.rig.rb    = rb;
    sim.rig.tau   = tau;
    sim.rig.r     = r;

    sim.rig.t_inx = t_inx;                    % end index of simulation
    sim.rig.t     = time;                     % time matrix
    sim.rig.x_1   = x(:,lnp(1,1:3));          % position tip
    sim.rig.x_in  = x(:,lnp(49,1:3));         % input position end of source cable
    sim.rig.x_end = [x(t_inx,noderig(:,1))',x(t_inx,noderig(:,2))',x(t_inx,noderig(:,3))']; % configuration at distal position
    sim.rig.x     = reshape(x(:,noderig),length(x),length(noderig),3);  % full configuration matrix
    sim.rig.v_1   = xd(:,lnp(1,1:3));         % velocity of tip
    sim.rig.v_in  = xd(:,lnp(49,1:3));        % input velocity
    sim.rig.f_1   = fx(:,lnp(1,1:3));         % force tip
    sim.rig.f_in  = fxtot(:,lnp(49,1:3));     % input force
    sim.rig.f_fil = movmean(sim.rig.f_in(:,1),[20 20]);% filtered force characteristic
    sim.rig.t_sim = tend;                     % time for simulation
else
    % Run flexible analysis
    A = regexp( fileread('FLEX_CAT_S_3D_V2.dat'), '\n', 'split');
    A{201} = sprintf(append('TIMESTEP ',num2str(tlim),' 8000          %% time and number of steps'));
    fid = fopen('FLEX_CAT_S_3D_V3.dat', 'w');
    fprintf(fid, '%s\n', A{:});
    fclose(fid);

    tlim  = depth/(2*vel_sim);                           % time till end of insertion
    fileID2 = fopen('Channel.bin','w');
    fwrite(fileID2, [rb,ra,0;stiffn,0,0;c_w,0,0;mu_s_SLS,mu_k_SLS,v_brk;...
        vel_sim*2/vel_ins,0,0], 'double');fclose(fileID2);

    fileID4 = fopen('Simulation.bin','w');
    fwrite(fileID4, [tlim;vel_sim*2], 'double');fclose(fileID4);

    tstart = tic;
    spacar(1,'FLEX_CAT_S_3D_V3');                        % forward dynamic analysis
    tend   = toc(tstart);
    spavisual('FLEX_CAT_S_3D_V3');                       % visualise outcome

    % SPACAR output
    iflex     = 1:2:49;                                  % nodes
    nodeflex  = lnp(iflex,1:3);
    tval      = min(abs(time-tlim));
    t_inx     = find(abs(time-tlim)==tval);

    sim.flex.x_c   = x_c;
    sim.flex.ra    = ra;
    sim.flex.rb    = rb;
    sim.flex.tau   = tau;
    sim.flex.r     = r;

    sim.flex.t_inx = t_inx;                    % end index of simulation
    sim.flex.t     = time;                     % time matrix
    sim.flex.x_1   = x(:,lnp(1,1:3));          % position tip
    sim.flex.x_in  = x(:,lnp(49,1:3));         % input position end of source cable
    sim.flex.x_end = [x(t_inx,nodeflex(:,1))',x(t_inx,nodeflex(:,2))',x(t_inx,nodeflex(:,3))']; % configuration at distal position
    sim.flex.x     = reshape(x(:,nodeflex),length(x),length(nodeflex),3);  % full configuration matrix
    sim.flex.v_1   = xd(:,lnp(1,1:3));         % velocity of tip
    sim.flex.v_in  = xd(:,lnp(49,1:3));        % input velocity
    sim.flex.f_1   = fx(:,lnp(1,1:3));         % force tip
    sim.flex.f_in  = fxtot(:,lnp(49,1:3));     % input force
    sim.flex.f_fil = movmean(sim.flex.f_in(:,1),[20 20]);% filtered force characteristic
    sim.flex.t_sim = tend;                     % time for simulation
end

%% Example of plotting results
% Plot results of SPACAR simulation
figure(5)
plot(sim.(sim_mod).x_in(:,1)+0.12,sim.(sim_mod).f_fil(:,1),'LineWidth',2,'Color','r');
hold on
xlim([0 0.12]);
ylim([0 30]);
xlabel('Insertion depth (m)');
ylabel('Insertion force (N)');
set(gca, 'FontName', 'Bahnschrift light','FontSize',12)