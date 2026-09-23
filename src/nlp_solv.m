%% Optimisation de trajectoire — Bateau 2D 
%  CasADi + IPOPT
% Modele d'etat vitesses cartesiennes
% Ne fonctionne pas correctement

clear; clc; close all;
import casadi.*

load("data/params.mat")
p = params;
load("data/scenario.mat")
s = scen;
addpath("./model")

%% --- Paramètres du NLP ---
Te_max   = 0.4;   % Resolution temporelle de la simulation [s]
Tf_init  = 30;  % initialisation Tf (s)
Tf_min   = 5;
Tf_max   = 30;
N        = 50; % nombre d'intervalles

% Bornes commandes
T_max     = p.T_max;
theta_max = p.theta_max;
u_min     = 1e-2;
u_max     = sqrt(T_max/p.Fx);    % vitesse max theorique = T_max/f 

% Matrice de pondération des commandes
W = eye(2); % T et theta normalises
Q= diag([0.01 0.01 1/4]);

% Cible
xt = scen.target_pos(1); yt = scen.target_pos(2);

%% --- Déclaration des variables symboliques CasADi ---

% Tf libre
Tf = MX.sym('Tf');
h  = Tf / N;

% États : [x, y, phi, vx, vy, r]  taille 6 x (N+1)
X  = MX.sym('X',  6, N+1);

% Commandes : [T_norm, theta_norm]  taille 2 x N
U  = MX.sym('U',  2, N);

%% --- Fonction de dynamique ---
% Entrée : état x (6x1), commande uc (2x1)
% Sortie : xdot (6x1)

x_s   = MX.sym('x_s',  6);
uc_s  = MX.sym('uc_s', 2);

phi_s   = x_s(3);
u_s    = x_s(4);
v_s    = x_s(5);
r_s     = x_s(6);
T_s     = uc_s(1)*T_max;
theta_s = uc_s(2)*theta_max;

% Norme et direction de la vitesse
alpha_s   = atan2(v_s, u_s);
beta_s    = alpha_s - theta_s;

cphi = cos(phi_s);
sphi = sin(phi_s);
R_s  = [cphi -sphi; sphi cphi];

% Vitesse
V        = R_s*[u_s; v_s];
V_norm_s = sqrt(u_s^2+v_s^2);

% Forces
Fg = p.rho*p.S*sin(2*beta_s)*V_norm_s^2; % Force de gouverne (Norme)
FT_vec = T_s*[1; 0];  % Force de gouverne dans Rb
Fg_vec = Fg*[sin(theta_s); -cos(theta_s)];
Ff_vec = [-p.Fx*u_s*abs(u_s); -p.Fy*v_s*abs(v_s)];
% Acceleration (calcul avec des tableaux)
a_s = (FT_vec + Fg_vec + Ff_vec)/p.m + [v_s*r_s; -u_s*r_s];

f_dyn = Function('f_dyn', {x_s, uc_s}, { ...
    vertcat( ...
        cphi*u_s - sphi*v_s, ...
        sphi*u_s + cphi*v_s, ...
        r_s, ...
        a_s(1), ...
        a_s(2), ...
        (Fg * p.Lg * cos(theta_s) - p.Yg * r_s) / p.I  ...
    )});

%% --- Construction du NLP ---

J    = 0;           % critère
g    = {};          % contraintes d'égalité (collocation)
g_lb = {};          % bornes inférieures
g_ub = {};          % bornes supérieures

% Pour un indice i, la valeur de g(i) doit etre comprise entre
% g_lb(i) et g_ub(i)

for k = 1:N
    xk  = X(:, k);
    xk1 = X(:, k+1);
    uk  = U(:, k);

    zk  = [xk(1)-xt; xk(2)-yt; atan2(xk(5),xk(4))-xk(3)];

    % Critère : distance à la cible + régularisation
    J = J + (zk'*Q*zk + ...
             uk'*W*uk) * h;

    % Collocation trapézoïdale
    fk  = f_dyn(xk,  uk);
    fk1 = f_dyn(xk1, uk);
    col = xk1 - xk - (h/2) * (fk + fk1);

    g{end+1}    = col;
    g_lb{end+1} = zeros(6,1);
    g_ub{end+1} = zeros(6,1);
end

% Terme terminal sur la cible
J = J + ((X(1,N+1)-xt)^2 + (X(2,N+1)-yt)^2);

% Contraintes d'obstacles (tous les noeuds)
for k = 1:N+1
    xk = X(:,k);
    for i = 1:scen.n_obs
        dist2 = (xk(1)-scen.obs(i,1))^2 + (xk(2)-scen.obs(i,2))^2;
        g{end+1}    = dist2;
        g_lb{end+1} = (scen.r_obs+scen.ecart)^2;
        g_ub{end+1} = inf;
    end
end

%% Conditions aux limites (égalités)
% État initial
g{end+1}    = X(:,1) - [0; 0; pi/4; 0; 0; 0];
g_lb{end+1} = zeros(6,1);
g_ub{end+1} = zeros(6,1);

% État final : position uniquement
g{end+1}    = X(1:2, N+1) - [xt; yt];
g_lb{end+1} = zeros(2,1);
g_ub{end+1} = zeros(2,1);

% Assemblage vecteur de décision
%w     = {Tf, reshape(X, [], 1), reshape(U, [], 1)};
w_vec = vertcat(Tf, vec(X), vec(U)); %vertcat(w{:});
g_vec = vertcat(g{:});
g_lb_vec = vertcat(g_lb{:});
g_ub_vec = vertcat(g_ub{:});

%% Bornes sur les variables de décision

% Tf
w_lb = Tf_min;
w_ub = Tf_max;

% États X : [x, y, phi, vx, vy, r] x (N+1)
x_lb = [-inf; -inf; -inf; -u_max; -u_max; -inf];
x_ub = [ inf;  inf;  inf;  u_max; u_max;  inf];
w_lb = [w_lb; repmat(x_lb, N+1, 1)];
w_ub = [w_ub; repmat(x_ub, N+1, 1)];

% Commandes U : [T, theta] x N
u_lb = [0; -1];
u_ub = [1;  1];
w_lb = [w_lb; repmat(u_lb, N, 1)];
w_ub = [w_ub; repmat(u_ub, N, 1)];

%% Point initial
% Commande d'init : poussée constante, gouverne nulle
uc_init = [0.5; 0];   % [T_norm, theta_norm]

% Fonction d'état avec commande figée (dénormalisée dans f_dyn)
ode_init = @(t, x) full(f_dyn(x, uc_init));

% Cap initial vers la cible
phi0_init = atan2(yt, xt);
x0_init   = [0; 0; phi0_init; 1e-2; 0; 0];

t_grid = linspace(0, Tf_init, N+1);   % grille NLP

[~, X_fwd] = ode45(ode_init, t_grid, x0_init, ...
                   odeset('RelTol',1e-4,'AbsTol',1e-6));
X_fwd = X_fwd';   % (6 x N+1)


s_vec = linspace(0, 1, N+1);
X_fwd(1,:) = s_vec * xt;
X_fwd(2,:) = s_vec * yt;
X_fwd(3,:) = atan2(yt - X_fwd(2,:), xt - X_fwd(1,:));
X_fwd(3,end) = X_fwd(3,end-1);

u_init        = zeros(2, N);
u_init(1,:)   = 0.5;   % T_norm

w0 = [Tf_init; reshape(X_fwd, [], 1); reshape(u_init, [], 1)];

%% --- Solveur IPOPT ---
nlp  = struct('x', w_vec, 'f', J, 'g', g_vec);
opts = struct();
opts.ipopt.max_iter             = 2000;
opts.ipopt.tol                  = 1e-5;       % 1e-6 inutile ici, gain ~20% itérations
opts.ipopt.constr_viol_tol      = 1e-5;

opts.ipopt.nlp_scaling_method   = 'gradient-based';

% Stratégie de barrière : adaptive >> monotone sur NLP non-convexes
opts.ipopt.mu_strategy          = 'adaptive';
opts.ipopt.mu_init              = 1e-1;       % barrière initiale plus agressive

% Solveur linéaire
opts.ipopt.linear_solver        = 'mumps';    % remplacer par 'ma57' si HSL dispo
opts.ipopt.mumps_pivtol         = 1e-4;       % pivotage plus robuste (défaut 1e-6)

% Reconstruction de Hessien : exact par défaut via CasADi AD
% N'activez limited-memory QUE si >2000 variables N * (size(X) + size(U))
% opts.ipopt.hessian_approximation = 'limited-memory';

opts.ipopt.print_level          = 5;

solver = nlpsol('solver', 'ipopt', nlp, opts);

sol = solver('x0',  w0, ...
             'lbx', w_lb, 'ubx', w_ub, ...
             'lbg', g_lb_vec, 'ubg', g_ub_vec);

%% --- Extraction de la solution ---
w_sol  = full(sol.x);

Tf_sol = w_sol(1);
idx_X  = 2 : 6*(N+1)+1;
idx_U  = 6*(N+1)+2 : length(w_sol);

X_sol  = reshape(w_sol(idx_X), 6, N+1);
U_sol  = reshape(w_sol(idx_U), 2, N);
t_sol  = linspace(0, Tf_sol, N+1);

fprintf('Tf optimal : %.2f s\n', Tf_sol);