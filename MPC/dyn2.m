function dxdt = dyn2(X, U, params)
% Modèle en repère cartésien : X = [x, y, phi, vx, vy, r]
% U = [T, theta]

x   = X(1);
y   = X(2);
phi = X(3);
vx  = X(4);
vy  = X(5);
r   = X(6);

T     = U(1);
theta = U(2);

m   = params(1);
f   = params(2);
Yg  = params(3);
I   = params(4);
Lg  = params(5);
rho = params(6);
S   = params(7);
Fx  = params(8);
Fy  = params(9);

% Norme et direction de la vitesse
u   = sqrt(vx^2 + vy^2);
psi = atan2(vy, vx);    % alpha + phi

% Angles
alpha = psi - phi;
beta  = alpha - theta;

cphi = cos(phi);
sphi = sin(phi);
calpha = cos(alpha);
salpha = sin(alpha);

% Matrice de rotation de Rb dans R0
R = [cphi -sphi; sphi cphi];

% Base du repere de Frenet
% eT = [cos(psi); sin(psi)];
% eN = [-sin(psi); cos(psi)];

% Vitesses
v_lon = u*calpha;
v_lat = u*salpha;

% Norme de la force de la gouverne
Fg = rho*S*sin(2*beta)*u^2;

% Forces
FT = T*[1; 0]; % Poussee dans b
%Fg_vec = Fg*[sin(phi-theta); -cos(phi-theta)];
Fg_vec = Fg*[sin(theta); -cos(theta)]; % Force de la gouverne dans Rb
%Ff = -f*u^2*eT;
Ff = [-Fx*v_lon*abs(v_lon); -Fy*v_lat*abs(v_lat)];

% Acceleration
a = R*(FT + Fg_vec + Ff)/m;

% Equations d'etat
dx  = vx;
dy  = vy;
dphi = r;
dvx = a(1);
dvy = a(2);
dr  = (Fg*Lg*cos(theta) - Yg*r)/I;

dxdt = [dx; dy; dphi; dvx; dvy; dr];
end