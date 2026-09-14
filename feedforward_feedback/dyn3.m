function dxdt = dyn3(t, x, T_cmd, theta_cmd, m, I, f, Lg, g, rho, S)
phi = x(3);
u  = x(4);
v  = x(5);
r   = x(6);

T     = T_cmd(t);
theta = theta_cmd(t);

% Norme et direction de la vitesse
V   = sqrt(v^2 + u^2);
psi = atan2(v, u);    % alpha + phi

% Angles
alpha = psi - phi;
beta  = alpha - theta;

% Force de la gouverne
Fg = rho*S*sin(2*beta)*u^2;

% Base du repere de Frenet
eT = [cos(alpha); sin(alpha)];
eN = [-sin(alpha); cos(alpha)];

% Forces
FT = T*[1; 0];
Fg_vec = Fg*[sin(theta); -cos(theta)];
Ff = -f*V^2*eT;

% Acceleration (calcul avec des tableaux)
a = (FT + Fg_vec + Ff)/m;

% Equations d'etat
dx  = u*cos(phi)+v*sin(phi);
dy  = u*sin(phi)+v*cos(phi);
dphi = r;
du = a(1)+r*v;
dv = a(2)-r*u;
dr  = (Fg*Lg*cos(theta) - g*r)/I;

dxdt = [dx; dy; dphi; du; dv; dr];
end