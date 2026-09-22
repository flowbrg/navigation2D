%% Calcul du gain feedback
% Linéarisation du modèe autour de
% Xeq = [x0, y0, phi0, u0, 0, 0]'
% Ueq = [T0, 0]
% Avec T0 = Fx u0^2
% Calcul d'un gain LQR sur modele augmente de ksi = int(en)


clear; clc; close all;

load("./data/params.mat")

u0 = linspace(0.1,8,20); % [m/s]
K_du = ones(1, length(u0));


Q_lon = diag([1/4, 1/4]);
R_lon = (1/400)^2;

for i = 1:length(u0)
    
    u = u0(i);
    % Coefficients
    af  = 2*params.Fx*u/params.m;
    
    % Sous-systeme longitudinal
    A_lon = [0,  1;
             0, -af];
    B_lon = [0; 1/params.m];
    
    K_lon = lqr(A_lon, B_lon, Q_lon, R_lon);
    
    K_du(i) = K_lon(2);
end

K = struct();

K.K_el = 200;
K.K_du = K_du;
K.K_lat = [-0.01, -2.32, -0.26, -0.22];
K.u0 = u0;

save("./lqr.mat", "K");

