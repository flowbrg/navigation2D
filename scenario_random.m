%% Paramètres de la scène
R_obs = 0.3;    % rayon des poteaux
ecart = 1;    % distance de sécurité aux obstacles

% Départ et cible
start_pos  = [0, 0];
target_pos = [10, 10];

N=12;

% Centres des obstacles
obs = 10 * rand(N, 2);
n_obs = N;

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
