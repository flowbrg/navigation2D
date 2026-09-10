%% Paramètres de la scène
R_obs = 0.3;    % rayon des poteaux
ecart = 1;    % distance de sécurité aux obstacles

% Départ et cible
start_pos  = [0, 0];
target_pos = [10, 10];

% Centres des obstacles
obs = [4.0, 8.0; 6.0, 10; 5.0, 0.0; 5.0, 0.8; 5.0, 1.6; 5.0, 2.4;
    5.0, 3.2; 5.0, 4.0; 8.0, 8.0; 6.0, 5.0; 7.0, 9.0; 6.0, 6.0];
n_obs = size(obs, 1);

%% Signaux de commande
% Paramètres des commandes
% Poussée, échelon unitaire
T0  = 500; 

% Couple, trapèze
% Montée linéaire de 0 à theta_max entre t1 et t2
% Palier à theta_max entre t2 et t3
% Descente linéaire de theta_max à 0 entre t3 et t4
theta_max   = pi/12;
t1 = 2; t2 = 4; t3 = 6; t4 = 8;

% Définition des commandes
T_cmd       = @(t) T0;  % Echelon unitaire
theta_cmd = @(t) ...
    (t >= t1 & t < t2)  .* (theta_max * (t - t1)/(t2 - t1)) + ...
    (t >= t2 & t < t3)  .* theta_max + ...
    (t >= t3 & t <= t4) .* (theta_max * (t4 - t)/(t4 - t3));
%theta_cmd = @(t) ...
%    (t >= t1 & t < t4)  .* theta_max;
