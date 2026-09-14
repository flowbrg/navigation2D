% Paramètres physiques du bateau - valeurs arbitraires
Cx  = 0.09;     % Demi corps profilé
Sx  = 1*0.5;    % Maitre-couple, largeur * tirant d'eau, largeur ~= longueur/3
rho = 1025;     % Eau de mer [kg/m3]
m   = 150;      % Masse [kg]
f   = (Cx/2)*rho*Sx;  % Frottement longitudinal [kg/m]
Yg  = 1e3;      % Amortissement lacet [kg.m2/rad/s]
I   = 80;       % Moment d'inertie vertical [kg.m2]
Lg  = 1.2;      % Bras de levier CG-gouverne [m]
S   = 0.03;     % Surface gouverne [m2]
Fx  = f;
Fy  = f*1e2;
params = [m, f, Yg, I, Lg, rho, S, Fx, Fy];