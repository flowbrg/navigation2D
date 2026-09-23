% Paramètres physiques du bateau - valeurs arbitraires
params = struct();

Cx  = 0.09;     % Demi corps profilé
Sx  = 1*0.5;    % Maitre-couple, largeur * tirant d'eau, largeur ~= longueur/3
params.rho = 1025;     % Eau de mer [kg/m3]
params.m   = 150;      % Masse [kg]
f   = (Cx/2)*params.rho*Sx;  % Frottement longitudinal [kg/m]
params.Yg  = 1e3;      % Amortissement lacet [kg.m2/rad/s]
params.I   = 80;       % Moment d'inertie vertical [kg.m2]
params.Lg  = 1.2;      % Bras de levier CG-gouverne [m]
params.S   = 0.03;     % Surface gouverne [m2]
params.Fx  = f;        % Coefficient de frottement longitudinal [kg/m]
params.Fy  = f*1e3;    % Coefficient de frottement lateral [kg/m]

% Paramètres descommandes
params.T_max     = 700;    % Saturation de la poussée [N]
params.theta_max = pi/4;   % Saturation de l'angle de gouverne [rad]

save("./params.mat", "params");