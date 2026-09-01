%% Visualisation de la zone — Bateau 2D
clear; clc; close all;

% --- Paramètres de la scène ---
R_obs = 0.3;    % rayon des poteaux
ecart = 1;    % distance de sécurité aux obstacles

% Départ et cible
start_pos  = [0, 0];
target_pos = [10, 10];

% Centres des obstacles
obs = [4.0, 8.0; 6.0, 10; 5.0, 0.0; 5.0, 0.8; 5.0, 1.6; 5.0, 2.4;
    5.0, 3.2; 5.0, 4.0; 8.0, 8.0; 6.0, 5.0; 7.0, 9.0; 6.0, 6.0];
n_obs = size(obs, 1);

% Paramètres physiques du bateau - valeurs arbitraires
m   = 150;      % Masse [kg]
f   = 15;       % Frottement longitudinal [kg/m]
g   = 10;       % Amortissement lacet [kg·m²/rad/s]
I   = 80;       % Moment d'inertie vertical [kg·m²]
Lg  = 1.2;      % Bras de levier CG-gouverne [m]
rho = 1025;     % Eau de mer [kg/m³]
S   = 0.03;     % Surface gouverne [m²]