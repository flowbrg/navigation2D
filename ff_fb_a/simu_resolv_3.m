%% Simulation navigation 2D basique
% Modele d'etat avec vitesses cartesiennes dans le repere du bateau
% X = [x, y, phi, u, v, r]

clear; clc; close all;

init_params
init_lqr
scenario_random
%calc_traj
poly_traj_opt_t

%% Trajectoire
%Tf = traj_ref.Tf;
u0 = 1;
%p_x = traj_ref.cx;
%p_y = traj_ref.cy;

traj_fun = make_trajectory(p_x, p_y, Tf*1.1); % [x_ref, y_ref, dx_ref,... ]


% Initialisation du bateau à la bonne
X0 = [0; 0; pi/4; u0; 0; 0; 0];

ode_fun  = @(t,X) closed_loop3(X, traj_fun(t), K, params);
options  = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t_sol, X_sol] = ode45(ode_fun, [0, Tf], X0, options);

%% Affichage
theta_c = linspace(0, 2*pi, 100);

% Figure 1 : Trajectoire
figure('Name', 'T1 - Trajectoire'); hold on; axis equal; grid on;
%xlim([-15, 15]); ylim([-15, 6]);

for i = 1:n_obs
    fill(obs(i,1) + (R_obs+ecart)*cos(theta_c), ...
        obs(i,2) + (R_obs+ecart)*sin(theta_c), ...
        [1.0 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    fill(obs(i,1) + R_obs*cos(theta_c), ...
        obs(i,2) + R_obs*sin(theta_c), ...
        [0.8 0.2 0.2], 'EdgeColor', 'k');
end

xlabel('x (m)'); ylabel('y (m)');
plot(X_sol(:,1), X_sol(:,2), 'b-', 'LineWidth', 2, 'DisplayName', 'Trajectoire');
plot(polyval(p_x, t_sol), polyval(p_y, t_sol), 'r-', 'LineWidth', 2, 'DisplayName', 'Référence');

plot(xs, ys, 'gs', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', 'Départ');
plot(xt, yt, 'p',  'MarkerSize', 16, 'MarkerFaceColor', 'y', ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Cible');

% Orientation de la vitesse
n_arrows = 15;  % Nombre de flèches
idx = round(linspace(1, length(t_sol), n_arrows));
alpha_sol = atan2(X_sol(:,5), X_sol(:,4));
V_sol   = sqrt(X_sol(:,5).^2 + X_sol(:,4).^2);
quiver(X_sol(idx,1), X_sol(idx, 2), V_sol(idx).*cos(alpha_sol(idx)+X_sol(idx, 3)),  V_sol(idx).*sin(alpha_sol(idx)+X_sol(idx, 3)), 0, 'k', 'LineWidth', 1.2, 'DisplayName', 'Vitesse \alpha+\phi');
% Cap
quiver(X_sol(idx,1), X_sol(idx, 2), cos(X_sol(idx,3)),  sin(X_sol(idx,3)), 0, 'r', 'LineWidth', 1.2, 'DisplayName', 'Cap \phi');


%legend('Location', 'northwest');

% Figure 2 : Etats
figure('Name', 'T1 - États');

% cap
subplot(4,1,1);
plot(t_sol, rad2deg(X_sol(:,3)), 'r', 'LineWidth', 1.5);
ylabel('\phi (°)'); grid on; title('Cap');

% angle d'incidence
subplot(4,1,2);
plot(t_sol, rad2deg(alpha_sol), 'g', 'LineWidth', 1.5);
ylabel('\alpha (°)'); grid on; title("Angle d'incidence");

% vitesse longi
subplot(4,1,3);
plot(t_sol, V_sol, 'b', 'LineWidth', 1.5);
ylabel('u (m/s)'); grid on; title('Vitesse');

% vitesse de lacet
subplot(4,1,4);
plot(t_sol, X_sol(:,6), 'm', 'LineWidth', 1.5);
ylabel('r (rad/s)'); grid on; title('Vitesse de lacet');
xlabel('t (s)');