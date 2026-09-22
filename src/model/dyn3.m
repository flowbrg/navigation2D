function dxdt = dyn3(X, U, params)
% Modele en repere corps : X = [x, y, phi, u, v, r]
% U = [T, theta]
% Forces exprimees dans Rb

x   = X(1);
y   = X(2);
phi = X(3);
u   = X(4);   % vitesse longitudinale dans Rb
v   = X(5);   % vitesse laterale dans Rb
r   = X(6);   % vitesse de lacet

T     = U(1);
theta = U(2);

m   = params.m;
Yg  = params.Yg;
I   = params.I;
Lg  = params.Lg;
rho = params.rho;
S   = params.S;
Fx  = params.Fx;
Fy  = params.Fy;

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