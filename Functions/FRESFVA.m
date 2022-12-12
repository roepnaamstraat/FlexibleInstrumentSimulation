function [time,sig,f] = FRESFVA(t,ne,le,e,ep,nx,lnp,x,xp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORCE INPUT - RIGID ANALYSIS OF CATHETER IN S-SHAPED CHANNEL 
% SUMMARY - This file takes a configuration and determines nodal loads. 
% INPUT   - x:   coordinate matrix of nodes
%         - xp:  velocity matrix of nodes
%         - lnp: location matrix of nodes
%         - 'CentreLine.bin': centreline of channel
%         - 'Channel.bin':    channel and contact properties
%         - 'Instrument.bin': instrument properties
% OUTPUT  - f:   matrix of prescribed loads on the nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Channel properties
fid = fopen('CentreLine.bin', 'r');
x_c  = fread(fid,'double'); fclose(fid);
x_c  = reshape(x_c,[length(x_c)/3,3]); % centreline

fid2 = fopen('Channel.bin', 'r');
Amat  = fread(fid2,'double'); fclose(fid2);
Amat  = reshape(Amat,[length(Amat)/3,3]);

rb      = Amat(1,1);            % radius channel
ra      = Amat(1,2);            % radius transition zone
stiffn  = Amat(2,1);            % wall stiffness
c_w     = Amat(3,1);            % wall damping coefficient
mu_s    = Amat(4,1);            % static friction coefficient
mu_k    = Amat(4,2);            % kinetic friction coefficient
v_brk   = Amat(4,3);            % breakaway velocity
vel_rat = Amat(5,1);            % velocity ratio

% Instrument properties
fid3 = fopen('Instrument.bin', 'r');
ro  = fread(fid3,'double'); fclose(fid3);

% Contact regions
a       = ra-ro;   
b       = rb-ro;   

% Compute forces and moments
% Preallocation
fx  = [(1:4:45)',1*ones(12,1),zeros(12,1)];   
fy  = [(1:4:45)',2*ones(12,1),zeros(12,1)]; 
fz  = [(1:4:45)',3*ones(12,1),zeros(12,1)]; 
   
M1  = [(2:4:46)',1*ones(12,1),zeros(12,1)]; 
M2  = [(2:4:46)',2*ones(12,1),zeros(12,1)];    
M3  = [(2:4:46)',3*ones(12,1),zeros(12,1)]; 
M4  = [(2:4:46)',4*ones(12,1),zeros(12,1)];  

% Determine nodal loads for all nodes 
j   = 1;
for i = 1:12
    node     = i*4-3;
    
    % Positions and velocities
    pos       = x(lnp(node,1:3));         % contains x,y,z coordinates
    pos2      = x(lnp((i+1)*4-3,1:3));
    v_o       = xp(lnp(node,1:3))/vel_rat;% contains linear velocity in global frame
    q         = x(lnp(node+1,:));         % contains angular coordinates in Euler parameters
    qp        = xp(lnp(node+1,:));        % contains angular velocity in Euler parameters 
    qtilde    = [ 0   -q(4)  q(3);
                  q(4) 0    -q(2);
                 -q(3) q(2)  0   ];
    Q_bar     = [q(1)    -q(2:4);
                 q(2:4)'  q(1)*eye(3)-qtilde];
    omega     = 2*Q_bar'*qp';             % angular velocity in global fixed ref. frame
    omega     = omega(2:4)/vel_rat;
    
    % Determine contact triad  
    [p_p,d_c,~] = distance2curve(x_c,pos);
    nvec      = (pos-p_p)/norm(pos-p_p);  % normal vector
    v_n       = (v_o*nvec')*nvec;         % normal velocity at CP
    tvec      = (v_o-v_n)/norm(v_o-v_n);  % tangent vector
    
    % Correct for angle at tip
    t_insvec  = (pos2-pos)/norm(pos2-pos);% tip vector
    bvec      = cross(t_insvec,nvec);
    t_pervec  = cross(bvec,t_insvec);
    t_pervec  = t_pervec/norm(t_pervec');
    ro_tip    = t_pervec*nvec'*ro;

    % Definitions of contact zones
    if i == 1        % tip
        a  = ra-ro_tip;
        b  = rb-ro_tip;
        ro = ro_tip;
    end
    
    % Tangential velocity
    v_omega   = cross(omega',ro*nvec');
    v_t       = v_o*tvec'+v_omega*tvec';
    
    % Define a custom force based on contact zone
    ksi       = (d_c-a)/(b-a);       % normalisation
    if d_c < a                       % zone I
        fn_val = 0;
        ff_val = 0;
    elseif a < d_c && d_c < b        % zone II
        fn_val   = (stiffn/2)*(b-a)*ksi^2+c_w*(3-2*ksi)*ksi^2*v_n*nvec'; % normal force
        ff_val   = sqrt(2*exp(1))*(-fn_val*mu_k + fn_val*mu_s)*exp(-(v_t/(v_brk*sqrt(2)))^2)*v_t/(v_brk*sqrt(2)) +...
                   fn_val*mu_k*tanh(v_t/(v_brk/10));                     % friction force
    else                             % zone III
        fn_val   = stiffn*(b-a)*(ksi-1/2)+c_w*v_n*nvec';                 % normal force
        ff_val   = sqrt(2*exp(1))*(-fn_val*mu_k + fn_val*mu_s)*exp(-(v_t/(v_brk*sqrt(2)))^2)*v_t/(v_brk*sqrt(2)) +...
                   fn_val*mu_k*tanh(v_t/(v_brk/10));                     % friction force
    end
    
    F         = -fn_val*nvec-ff_val*tvec;  % force at node
    M_o       = cross(ro*nvec,F);          % moment at node
    T_bar     = 2*Q_bar*[0;M_o'];          % moment at node in Euler parameters
        
    fx(j,:)   = [node   1 F(1)];           % apply force in x-direction
    fy(j,:)   = [node   2 F(2)];           % apply force in y-direction
    fz(j,:)   = [node   3 F(3)];           % apply force in z-direction
    M1(j,:)   = [node+1 1 T_bar(1)];       % apply moment in 0-direction
    M2(j,:)   = [node+1 2 T_bar(2)];       % apply moment in 1-direction
    M3(j,:)   = [node+1 3 T_bar(3)];       % apply moment in 2-direction
    M4(j,:)   = [node+1 4 T_bar(4)];       % apply moment in 3-direction
    
    j = j+1;
end

f    = [fx;fy;fz;M1;M2;M3;M4];
f(isnan(f))=0;

sig  = []; % can be left empty (is used for specifying generalized stresses)
time = t;  % this line needs to be there; don't modify.
if rem(t,0.001) <= 1*10^(-5)
    t
end


end