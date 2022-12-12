function [X_td] = LiePose(X_t,u_t)
%Function to calculate new pose from current pose and control inputs
delta = 1;

vee_t = u_t(1);
k_t   = u_t(2);
t_t   = u_t(3);
    
w_t   = [vee_t*t_t;0;vee_t*k_t];
v_t   = [vee_t;0;0];
w_hat = [ 0     ,-w_t(3), w_t(2);
    w_t(3), 0     ,-w_t(1);
    -w_t(2), w_t(1), 0     ];

phi_t = eye(3)+w_hat/norm(w_t)*sin(norm(w_t)*delta)+...
    w_hat^2/(norm(w_t)^2)*(1-cos(norm(w_t)*delta));

exp_ksi = [phi_t     ,1/(norm(w_t)^2)*((eye(3)-phi_t)*cross(w_t,v_t)+w_t*w_t'*v_t*delta);
           zeros(1,3),1];
X_td  = X_t*exp_ksi;
end

