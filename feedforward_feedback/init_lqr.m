%% Calcul du gain feedback
% Linéarisation du modèe autour de
% Xep = [x0, y0, phi0, u0, 0, 0]'
% Ueq = [T0, 0]
% Avec T0 = Fx u0^2

clear; clc; close all;

init_params

u0 = 5; % [m/s]

% Parametres
m   = params(1);
Fx  = params(2);
Yg  = params(3);
I   = params(4);
Lg  = params(5);
rho = params(6);
S   = params(7);

% Coefficients
af  = 2*Fx*u0/m;
av  = 2*rho*S*u0/m;
ar  = Yg/I;
bv  = 2*rho*S*u0^2/m;
br  = 2*rho*S*Lg*u0^2/I;

% Sous-systeme longitudinal
A_lon = [0,  1;
         0, -af];
B_lon = [0; 1/m];

Q_lon = diag([1/4, 1/4]);
R_lon = (1/200)^2;

K_lon = lqr(A_lon, B_lon, Q_lon, R_lon);

% Sous-systeme lateral ---
A_lat = [0,    1,    0;
         0,   -av,  -u0;
         0,  br/Lg, -ar];
B_lat = [0; bv; -br];

Q_lat = diag([1, 1, 1/0.09]);
R_lat = (6/pi)^2;   % = 1/(pi/6)^2

K_lat = lqr(A_lat, B_lat, Q_lat, R_lat);

% Affichage des pôles en boucle fermée
fprintf('Pôles BF longitudinal : '); 
disp(eig(A_lon - B_lon*K_lon)');
fprintf('Pôles BF latéral      : '); 
disp(eig(A_lat - B_lat*K_lat)');