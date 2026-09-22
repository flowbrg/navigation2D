function [U_ff, X_ref] = feedforward3(sigma_r, params)
% Inversion analytique du modèles sous les hypothèses
% H1 : alpha ~ 0
% H2 : theta << 1

m   = params.m;
Yg  = params.Yg;
I   = params.I;
Lg  = params.Lg;
rho = params.rho;
S   = params.S;
Fx  = params.Fx;

V_min   = 1e-2;     % Saturation vitesse pour atan2

% Derivees de x_r et y_r
x_ref = sigma_r(1);
y_ref = sigma_r(2);
dx    = sigma_r(3);
dy    = sigma_r(4);
ddx   = sigma_r(5);
ddy   = sigma_r(6);
dddx  = sigma_r(7);
dddy  = sigma_r(8);

% Vitesse
V2  =  max(dx^2 + dy^2, V_min^2);
V   = sqrt(V2);

% Etats de reference
phi_ref = atan2(dy,dx); % non requis pour T_ff et theta_ff
u_ref   = V;
v_ref   = 0; % H1

du_ref  = (dx*ddx + dy*ddy)/ V;

r_ref   = (dx*ddy - dy*ddx)/V^2;

num_r   = dx*dddy - dy*dddx;
%num_r_dt_approx = 0; % terme croisé négligé si trajectoire lisse
dot_V2  = 2*(dx*ddx + dy*ddy);
dr_ref  = (num_r*V2 - (dx*ddy - dy*ddx)*dot_V2) / V2^2;

% Feedforward T
T_ff = +m*du_ref + Fx*u_ref*abs(u_ref);

% Feedforward theta
Fg_ref   = (I*dr_ref + Yg*r_ref) / Lg;
theta_ff = -Fg_ref / (2*rho*S*V2);

% Return
U_ff = [T_ff, theta_ff];
X_ref = [x_ref, y_ref, phi_ref, u_ref, v_ref, r_ref];

end