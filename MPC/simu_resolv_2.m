%% Optimisation de trajectoire — Bateau 2D 
%  CasADi + IPOPT
% Modele d'etat vitesses cartesiennes
% Ne fonctionne pas correctement

clear; clc; close all;
import casadi.*

init_params

%scenario
scenario_random2

%% --- Paramètres du NLP ---
Te_max   = 0.5;   % Resolution temporelle de la simulation [s]
Tf_init  = 30;  % initialisation Tf (s)
Tf_min   = 5;
Tf_max   = 50;
N        = Tf_max/Te_max; % nombre d'intervalles

% Bornes commandes
T_max     = 700;
theta_max = pi/3;
%u_min    = 0;
u_max     = sqrt(T_max/f);    % vitesse max theorique = T_max/f 

% Matrice de pondération des commandes
W = diag([(1/T_max)^2 (1/theta_max)^2]);
%W = diag([0.001 1]);      % poids régularisation commandes
Q= 0.01*eye(2);

% Cible
xt = target_pos(1); yt = target_pos(2);

%% --- Déclaration des variables symboliques CasADi ---

% Tf libre
Tf = MX.sym('Tf');
h  = Tf / N;

% États : [x, y, phi, vx, vy, r]  taille 6 x (N+1)
X  = MX.sym('X',  6, N+1);

% Commandes : [T, theta]  taille 2 x N
U  = MX.sym('U',  2, N);

%% --- Fonction de dynamique ---
% Entrée : état x (6x1), commande uc (2x1)
% Sortie : xdot (6x1)

x_s   = MX.sym('x_s',  6);
uc_s  = MX.sym('uc_s', 2);

phi_s   = x_s(3);
vx_s    = x_s(4);
vy_s    = x_s(5);
r_s     = x_s(6);
T_s     = uc_s(1);
theta_s = uc_s(2);

% Norme et direction de la vitesse
us_s    = sqrt(vx_s^2+vy_s^2);
psi_s   = atan2(vy_s, vx_s);

% Angles
alpha_s = psi_s - phi_s;
beta_s  = alpha_s - theta_s;

cphi_s = cos(phi_s);
sphi_s = sin(phi_s);
calpha_s = cos(alpha_s);
salpha_s = sin(alpha_s);

% Matrice de rotation de Rb dans R0
R_s = [cphi_s -sphi_s; sphi_s cphi_s];

% Vitesses
v_lon_s = us_s*calpha_s;
v_lat_s = us_s*salpha_s;

% Force de la gouverne
Fg_s = rho*S*sin(2*beta_s)*us_s^2;

% Forces
FT_s = T_s*[1; 0];
Fg_vec_s = Fg_s*[sin(theta_s); -cos(theta_s)];
Ff_s = [-Fx*v_lon_s*abs(v_lon_s); -Fy*v_lat_s*abs(v_lat_s)];

% Acceleration (calcul avec des tableaux)
a_s = R_s*(FT_s + Fg_vec_s + Ff_s)/m;

f_dyn = Function('f_dyn', {x_s, uc_s}, { ...
    vertcat( ...
        vx_s, ...
        vy_s, ...
        r_s, ...
        a_s(1), ...
        a_s(2), ...
        (Fg_s * Lg * cos(theta_s) - Yg * r_s) / I  ...
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

    zk  = [xk(1)-xt; xk(2)-yt];

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
    for i = 1:n_obs
        dist2 = (xk(1)-obs(i,1))^2 + (xk(2)-obs(i,2))^2;
        g{end+1}    = dist2;
        g_lb{end+1} = (R_obs+ecart)^2;
        g_ub{end+1} = inf;
    end
end

%% --- Conditions aux limites (égalités) ---
% État initial
g{end+1}    = X(:,1) - [0; 0; pi/4; 0; 0; 0];
g_lb{end+1} = zeros(6,1);
g_ub{end+1} = zeros(6,1);

% État final : position uniquement
g{end+1}    = X(1:2, N+1) - [xt; yt];
g_lb{end+1} = zeros(2,1);
g_ub{end+1} = zeros(2,1);

%% --- Assemblage vecteur de décision ---
%w     = {Tf, reshape(X, [], 1), reshape(U, [], 1)};
w_vec = vertcat(Tf, vec(X), vec(U)); %vertcat(w{:});
g_vec = vertcat(g{:});
g_lb_vec = vertcat(g_lb{:});
g_ub_vec = vertcat(g_ub{:});

%% --- Bornes sur les variables de décision ---

% Tf
w_lb = Tf_min;
w_ub = Tf_max;

% États X : [x, y, phi, vx, vy, r] x (N+1)
x_lb = [-inf; -inf; -inf; -u_max; -u_max; -inf];
x_ub = [ inf;  inf;  inf;  u_max; u_max;  inf];
w_lb = [w_lb; repmat(x_lb, N+1, 1)];
w_ub = [w_ub; repmat(x_ub, N+1, 1)];

% Commandes U : [T, theta] x N
u_lb = [0;        -theta_max];
u_ub = [T_max;     theta_max];
w_lb = [w_lb; repmat(u_lb, N, 1)];
w_ub = [w_ub; repmat(u_ub, N, 1)];

%% --- Point initial ---
% Interpolation linéaire en position, reste nul

% Vitesse cible estimée (heuristique)
dist_cible = sqrt(xt^2 + yt^2);
v_target_init = 1.2 * dist_cible / Tf_init;

x_init_traj = zeros(6, N+1);
for k = 1:N+1
    s = (k-1)/N;
    x_init_traj(1,k) = s * xt;
    x_init_traj(2,k) = s * yt;
    x_init_traj(3,k) = pi/4;
    x_init_traj(4,k) = v_target_init * cos(pi/4);
    x_init_traj(5,k) = v_target_init * sin(pi/4);

end
u_init = zeros(2, N);
u_init(1,:) = 0.5*T_max;             % poussée initiale modérée

w0 = [Tf_init; reshape(x_init_traj, [], 1); reshape(u_init, [], 1)];

%% --- Solveur IPOPT ---
nlp  = struct('x', w_vec, 'f', J, 'g', g_vec);
opts = struct();
opts.ipopt.max_iter        = 2000;
opts.ipopt.tol             = 1e-6;
opts.ipopt.print_level     = 5;

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

%% --- Figure 1 : Trajectoires idéale et réelle ---

theta_c = linspace(0, 2*pi, 100);

figure(1); clf;
hold on; axis equal; grid on;
xlim([-1,52]); ylim([-1,52]);

xlabel('x (m)');
ylabel('y (m)');
title(sprintf('Trajectoire optimale et trajectoire réelle  (Tf = %.1f s)', Tf_sol));

% Obstacles
for i = 1:n_obs
    fill(obs(i,1) + R_obs*cos(theta_c), ...
         obs(i,2) + R_obs*sin(theta_c), ...
         [0.8 0.2 0.2], ...
         'EdgeColor','k');
end

% Départ et cible
plot(0, 0, 'gs', ...
     'MarkerSize', 12, ...
     'MarkerFaceColor', 'g', ...
     'DisplayName', 'Départ');

plot(xt, yt, 'p', ...
     'MarkerSize', 16, ...
     'MarkerFaceColor', 'y', ...
     'MarkerEdgeColor', 'k', ...
     'DisplayName', 'Cible');

% Trajectoire idéale issue de la collocation
plot(X_sol(1,:), X_sol(2,:), ...
     'b-', 'LineWidth', 2, ...
     'DisplayName', 'Trajectoire idéale (collocation)');

% Commandes idéales interpolées
t_u = t_sol(1:end-1);

T_cmd_opt = @(t) interp1(t_u, U_sol(1,:), t, 'previous', 'extrap');
theta_cmd_opt = @(t) interp1(t_u, U_sol(2,:), t, 'previous', 'extrap');

% Simulation de la trajectoire réelle avec la dynamique
x0_real = [0; 0; pi/4; 0; 0; 0];

ode_fun_real = @(t,x) dyn2(x, ...
    [T_cmd_opt(t), theta_cmd_opt(t)], params);

options_real = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

[t_real, x_real] = ode45(ode_fun_real, [0 Tf_sol], x0_real, options_real);

% Trajectoire réelle
plot(x_real(:,1), x_real(:,2), ...
     'r-', 'LineWidth', 2, ...
     'DisplayName', 'Trajectoire réelle');

% Direction de la vitesse - trajectoire idéale
n_arr = 15;
idx = round(linspace(1, N+1, n_arr));

V_sol = sqrt(X_sol(4,idx).^2 + X_sol(5,idx).^2);

quiver(X_sol(1,idx), X_sol(2,idx), ...
       0.4 * X_sol(4,idx)./V_sol, ...
       0.4 * X_sol(5,idx)./V_sol, ...
       0, 'k', ...
       'HandleVisibility', 'off');

% Direction de la vitesse - trajectoire réelle
idx_real = round(linspace(1, length(t_real), n_arr));

V_real = sqrt(x_real(idx_real,4).^2 + x_real(idx_real,5).^2);

quiver(x_real(idx_real,1), x_real(idx_real,2), ...
       0.4 * x_real(idx_real,4)./V_real, ...
       0.4 * x_real(idx_real,5)./V_real, ...
       0, ...
       'Color', [0.8500 0.3250 0.0980], ...
       'HandleVisibility', 'off');

%legend('Location','northwest');

%% --- Figure 2 : Etats idéale vs réels ---

figure('Name', 'États');

% Etats idéaux interpolés sur le temps réel
phi_ideal = interp1(t_sol, X_sol(3,:), t_real, 'linear', 'extrap');
vx_ideal  = interp1(t_sol, X_sol(4,:), t_real, 'linear', 'extrap');
vy_ideal  = interp1(t_sol, X_sol(5,:), t_real, 'linear', 'extrap');
r_ideal   = interp1(t_sol, X_sol(6,:), t_real, 'linear', 'extrap');

% Vitesse
v_ideal = sqrt(vx_ideal.^2 + vy_ideal.^2);
v_real  = sqrt(x_real(:,4).^2 + x_real(:,5).^2);

% Angle d'attaque
alpha_ideal = atan2(vy_ideal, vx_ideal) - phi_ideal;
alpha_real  = atan2(x_real(:,5), x_real(:,4)) - x_real(:,3);

% Cap
subplot(4,1,1);
plot(t_real, rad2deg(phi_ideal), 'b', 'LineWidth', 1.5, ...
     'DisplayName', 'Idéal');
hold on;
plot(t_real, rad2deg(x_real(:,3)), 'r', 'LineWidth', 1.5, ...
     'DisplayName', 'Réel');
ylabel('\phi (°)');
grid on;
title('Cap');
legend('Location','best');

% Angle d'incidence
subplot(4,1,2);
plot(t_real, rad2deg(alpha_ideal), 'b', 'LineWidth', 1.5);
hold on;
plot(t_real, rad2deg(alpha_real), 'r', 'LineWidth', 1.5);
ylabel('\alpha (°)');
grid on;
title('Angle d''attaque');

% Vitesse
subplot(4,1,3);
plot(t_real, v_ideal, 'b', 'LineWidth', 1.5);
hold on;
plot(t_real, v_real, 'r', 'LineWidth', 1.5);
ylabel('v (m/s)');
grid on;
title('Vitesse longitudinale');

% Vitesse de lacet
subplot(4,1,4);
plot(t_real, r_ideal, 'b', 'LineWidth', 1.5);
hold on;
plot(t_real, x_real(:,6), 'r', 'LineWidth', 1.5);
ylabel('r (rad/s)');
grid on;
title('Vitesse de lacet');
xlabel('t (s)');

%% Figure 3 : Commandes
figure('Name', 'Commandes');

t_u = t_sol(1:end-1);
figure(3);
subplot(2,1,1); stairs(t_u, U_sol(1,:), 'b', 'LineWidth',1.5);
ylabel('T (N)'); grid on; title('Poussée');

subplot(2,1,2); stairs(t_u, rad2deg(U_sol(2,:)), 'r', 'LineWidth',1.5);
ylabel('\theta (°)'); grid on; title('Angle de gouverne');
xlabel('t (s)');