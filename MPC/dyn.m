function dxdt = dyn(X, U, params)
x       = X(1);
y       = X(2);
phi     = X(3);
alpha   = X(4);
u       = X(5);
r       = X(6);

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

beta = alpha - theta;   % direction de la vitesse
psi = alpha + phi;

calpha = cos(alpha);
salpha = sin(alpha);

% Matrice de rotation de Frenet dans b (eT, eN)' = R*(xb, yb)'
R = [calpha salpha; -salpha calpha];

v_lon = u*cos(alpha);
v_lat = u*sin(alpha);

% Forces (normes)
Fg = rho*S*sin(2*beta)*u^2; % Norme de la force de la gouverne
% Forces (vecteurs)
FT = T*[1; 0]; % Poussee dans b
Fg_vec = Fg*[sin(theta); -cos(theta)]; % Force de la gouverne dans b
%Ff = -f*u^2*[1; 0]; % Trainee dans le repere de Frenet
Ff = [-Fx*v_lon*abs(v_lon); -Fy*v_lat*abs(v_lat)];

% Acceleration
a = (R*(FT + Fg_vec + Ff))/m;

dx      = u*cos(psi);
dy      = u*sin(psi);
dphi    = r;
dalpha  = a(2)/u-r;
du      = a(1);
dr      = (Fg*Lg*cos(theta)-Yg*r)/I;

dxdt = [dx; dy; dphi; dalpha; du; dr];
end