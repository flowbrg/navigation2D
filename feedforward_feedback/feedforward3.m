function [T_ff, theta_ff] = feedforward3(t, sigma_r, params)
% Inversion analytique du modèles sous les hypothèses
% H1 : alpha ~ 0
% H2 : theta << 1

m   = params(1);
Yg  = params(3);
I   = params(4);
Lg  = params(5);
rho = params(6);
S   = params(7);
Fx  = params(8);

V_min   = 1e-2;     % Saturation vitesse pour atan2
theta_max = pi/3;   % Saturation angle de gouverne, limite physique pi/2

% Derivees de sigma_r
[dx, dy, ddx, ddy, dddx, dddy] = eval_poly_derivatives(t, sigma_r);

% Vitesse
V2  = dx^2 + dy^2;
V   = sqrt(U2);
V   = max(V, V_min);

% Etats de reference
% phi_ref = atan2(dy,dx) % non requis pour T_ff et theta_ff
u_ref   = V;

du_ref  = (dx*ddx + dy*ddy)/ V;

r_ref   = (dx*ddy - dy*ddx)/V^2;

num_r   = dx*dddy - dy*dddx;
dnum_r_dt_approx = 0; % terme croisé négligé si trajectoire lisse
dot_V2  = 2*(dx*ddx + dy*ddy);
dr_ref  = (num_r*V2 - (dx*ddy - dy*ddx)*dot_V2) / V2^2;

% Feedforward T
T_ff = -m*du_ref + Fx*u_ref*abs(u_ref);

% Feedforward theta
Fg_ref   = (I*dr_ref + Yg*r_ref) / Lg;
theta_ff = -Fg_ref / (2*rho*S*V2);

theta_ff = max(-theta_max, min(theta_max, theta_ff));

end