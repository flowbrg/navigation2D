function dxdt = dyn3(X, U, params)
% Modèle en repère corps : X = [x, y, phi, u, v, r]
% U = [T, theta]
% Forces exprimées dans Rb

x   = X(1);
y   = X(2);
phi = X(3);
u   = X(4);   % vitesse longitudinale dans Rb
v   = X(5);   % vitesse latérale dans Rb
r   = X(6);   % vitesse de lacet

T     = U(1);
theta = U(2);

m   = params(1);
%f   = params(2);
Yg  = params(3);
I   = params(4);
Lg  = params(5);
rho = params(6);
S   = params(7);
Fx  = params(8);
Fy  = params(9);

% Angles
alpha = atan2(v, u);   % angle d'incidence sur le bateau
beta  = alpha - theta; % angle d'incidence sur la gouverne

cphi = cos(phi);
sphi = sin(phi);

% Matrice de rotation de Rb dans R0
R = [cphi -sphi; sphi cphi];

% Vitesses
V = R*[u; v]; % [vx; vy]
V_norm = sqrt(u^2 + v^2);

% Forces
Fg = rho*S*sin(2*beta)*V_norm^2; % Force de gouverne (Norme)
FT_vec = T*[1; 0];  % Force de gouverne dans Rb
Fg_vec = Fg*[sin(theta); -cos(theta)];
Ff_vec = [-Fx*u*abs(u); -Fy*v*abs(v)];

% Acceleration
% m*(du/dt - v*r) = F_u  =>  du/dt = F_u/m + v*r
% m*(dv/dt + u*r) = F_v  =>  dv/dt = F_v/m - u*r
a = (FT_vec + Fg_vec + Ff_vec)/m + [v*r; -u*r];

% Equations d'etat
dx   = V(1);
dy   = V(2);
dphi = r;
du   = a(1);
dv   = a(2);
dr   = (Fg * Lg * cos(theta) - Yg * r) / I;

dxdt = [dx; dy; dphi; du; dv; dr];
end