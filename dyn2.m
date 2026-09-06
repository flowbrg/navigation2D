function dxdt = dyn2(t, x, T_cmd, theta_cmd, m, I, f, Lg, g, rho, S)
phi = x(3);
vx  = x(4);
vy  = x(5);
r   = x(6);

T     = T_cmd(t);
theta = theta_cmd(t);

% Norme et direction de la vitesse
u   = sqrt(vx^2 + vy^2);
psi = atan2(vy, vx);    % alpha + phi

% Angles
alpha = psi - phi;
beta  = alpha - theta;

% Force de la gouverne
Fg = rho*S*sin(2*beta)*u^2;

% Base du repere de Frenet
eT = [cos(psi); sin(psi)];
eN = [-sin(psi); cos(psi)];

% Forces
FT = T*[cos(phi); sin(phi)];
Fg_vec = Fg*[sin(phi-theta); -cos(phi-theta)];
Ff = -f*u^2*eT;

% Acceleration (calcul avec des tableaux)
a = (FT + Fg_vec + Ff)/m;

% Equations d'etat
dx  = vx;
dy  = vy;
dphi = r;
dvx = a(1);
dvy = a(2);
dr  = (Fg*Lg*cos(theta) - g*r)/I;

dxdt = [dx; dy; dphi; dvx; dvy; dr];
end