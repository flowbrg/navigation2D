%% Calcul du gain feedback
% Linéarisation du modèe autour de
% Xeq = [x0, y0, phi0, u0, 0, 0]'
% Ueq = [T0, 0]
% Avec T0 = Fx u0^2
% Calcul d'un gain LQR sur modele augmente de ksi = int(en)


%clear; clc; close all;

load("./data/params.mat")

Fx  = params.Fx;
m   = params.m;
S   = params.S;
rho = params.rho;
Lg  = params.Lg;
I   = params.I;
Yg  = params.Yg;


u0 = linspace(0.1,8,200); % [m/s]
K = ones(6, length(u0));


Q_lon = diag([1/4, 1/4]);
R_lon = (1/400)^2;
Q_lat = diag([1e-3, (2/pi)^2, 1/4, 1/0.09]);
R_lat   = (6/pi)^2;   % = 1/(pi/6)^2

for i = 1:length(u0)
    
    u = u0(i);
    % Coefficients
    af  = 2*Fx*u/m;
    av  = 2*rho*S*u/m;
    ar  = Yg/I;
    bv  = 2*rho*S*u^2/m;
    br  = 2*rho*S*Lg*u^2/I;
    
    % Sous-systeme longitudinal
    A_lon = [0,  1;
             0, -af];
    B_lon = [0; 1/m];
    
    K_lon = lqr(A_lon, B_lon, Q_lon, R_lon);
    
    % Sous-systeme lateral
    A_lat   = [0,  u,   1,    0;
               0,   0,   0,    1;
               0,   0,  -av,  -u;
               0,   0, br/Lg, -ar];
    
    B_lat   = [0;   0;   bv;  -br];
   
    K_lat = lqr(A_lat, B_lat, Q_lat, R_lat);
    
    K(:,i) = [K_lon'; K_lat'];
end

figure;
subplot(3,2,1)
plot(u0,K(1,:));
subplot(3,2,2)
plot(u0,K(2,:));
subplot(3,2,3)
plot(u0,K(3,:));
subplot(3,2,4)
plot(u0,K(4,:));
subplot(3,2,5)
plot(u0,K(5,:));
subplot(3,2,6)
plot(u0,K(6,:));

% Affichage des poles en boucle fermee
%fprintf('Pôles BF longitudinal : '); 
%disp(eig(A_lon - B_lon*K_lon)');
%fprintf('Pôles BF latéral      : '); 
%disp(eig(A_lat_a- B_lat_a*K_lat_a)');