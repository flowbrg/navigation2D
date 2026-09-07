%% Paramètres de la scène
R_obs = 1;    % rayon des poteaux
ecart = 1;    % distance de sécurité aux obstacles

% Départ et cible
start_pos  = [0, 0];
target_pos = [50, 50];

N=10;

% Centres des obstacles
obs = 5*ones(N,2) + 40 * rand(N, 2);
n_obs = N;