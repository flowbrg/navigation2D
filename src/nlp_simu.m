%% Calcul de trajectoire ideale a priori
%  CasADi + IPOPT

clear; clc; close all;

nlp_solv

%% Simulation de la solution sur le modèle "réel" (sans collocation)
% Commandes idéales interpolées
t_u = t_sol(1:end-1);

T_cmd_opt = @(t) interp1(t_u, U_sol(1,:), t, 'previous', 'extrap');
theta_cmd_opt = @(t) interp1(t_u, U_sol(2,:), t, 'previous', 'extrap');

% Simulation de la trajectoire réelle avec la dynamique
x0_real = [0; 0; pi/4; 0; 0; 0];

ode_fun_real = @(t,x) dyn3(x, ...
    [T_cmd_opt(t)*T_max, theta_cmd_opt(t)*theta_max], params);

options_real = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

[t_real, x_real] = ode45(ode_fun_real, [0 Tf_sol], x0_real, options_real);

%% Figure 1 : Trajectoires idéale et réelle

theta_c = linspace(0, 2*pi, 100);

figure(1); clf;
hold on; axis equal; grid on;
xlim([-1,52]); ylim([-1,52]);

xlabel('x (m)');
ylabel('y (m)');
title(sprintf('Trajectoire optimale et trajectoire réelle  (Tf = %.1f s)', Tf_sol));

% Obstacles
for i = 1:scen.n_obs
    fill(scen.obs(i,1) + scen.r_obs*cos(theta_c), ...
         scen.obs(i,2) + scen.r_obs*sin(theta_c), ...
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