function [dT, dtheta] = feedback3(X, X_ref, K)

% Etats de reference
phi_ref = X_ref(3);
u_ref   = X_ref(4);
r_ref   = X_ref(6);

% Ecart de position
e_pos = [X(1) - X_ref(1);
         X(2) - X_ref(2)];

R = [cos(phi_ref),  sin(phi_ref);
         -sin(phi_ref), cos(phi_ref)];

e_traj = R * e_pos;
e_l = e_traj(1);
e_n = e_traj(2);

% Ecart etats
dphi = X(3) - phi_ref;
du   = X(4) - u_ref;
dv   = X(5);
dr   = X(6) - r_ref;

[K_lon, K_lat] = eval_K(K,X(4));

dT = -K_lon*[e_l; du];
dtheta = -K_lat * [e_n; dphi; dv; dr];

end