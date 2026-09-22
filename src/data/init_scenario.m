%% Paramètres de la scène


%clear; clc; close all;


% Départ et cible
start_pos  = [0, 0];
target_pos = [50, 50];

N=10;

% Obstacles

r_obs = 1;    % radius of obstacles
ecart = 1;    % safe distance

obs = 5*ones(N,2) + 40 * rand(N, 2); % [5;45] by [5;45]
n_obs = N;

scen = struct();
scen.start_pos = start_pos;
scen.target_pos = target_pos;
scen.r_obs = r_obs;
scen.ecart = ecart;
scen.obs = obs;
scen.n_obs = n_obs;

save('./scenario.mat', 'scen')
