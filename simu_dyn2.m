%% Simulation navigation 2D basique
% Modele d'etat avec vitesses cartesiennes
% X = [x, y, phi, vx, vy, r]

clear; clc; close all;

init_params

scenario

%% Conditions initiales
% [x, y, phi, vx, vy, r)
x0 = [0; 0; pi/4; 0; 0; 0];    % Cap initial de 45°, vers la cible

% Intégration ode45, explicit runge-kutta
t_span          = [0,20];

ode_fun         = @(t,x) dyn2(t, x, T_cmd, theta_cmd, m, I, f, Lg, g, rho, S);
options         = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t_sol, x_sol]  = ode45(ode_fun, t_span, x0, options);

T_sol           = arrayfun(T_cmd, t_sol);
theta_sol       = arrayfun(theta_cmd, t_sol);

%% Figure 1 : Trajectoire
figure(1); hold on; axis equal; grid on;
%xlim([-15, 15]); ylim([-15, 6]);

xlabel('x (m)'); ylabel('y (m)');
plot(x_sol(:,1), x_sol(:,2), 'b-', 'LineWidth', 2, 'DisplayName', 'Trajectoire');

% Orientation de la vitesse
n_arrows = 15;  % Nombre de flèches
idx = round(linspace(1, length(t_sol), n_arrows));
quiver(x_sol(idx,1), x_sol(idx, 2), x_sol(idx,4),  x_sol(idx,5), 0, 'k', 'LineWidth', 1.2, 'DisplayName', 'Vitesse \alpha+\phi');
% Cap
quiver(x_sol(idx,1), x_sol(idx, 2), cos(x_sol(idx,3)),  sin(x_sol(idx,3)), 0, 'r', 'LineWidth', 1.2, 'DisplayName', 'Cap \phi');


legend('Location', 'northwest');

%% Figure 2 : Etats
figure('Name', 'États');

% cap
subplot(4,1,1);
plot(t_sol, rad2deg(x_sol(:,3)), 'r', 'LineWidth', 1.5);
ylabel('\phi (°)'); grid on; title('Cap');

% angle d'incidence
subplot(4,1,2);
plot(t_sol, rad2deg(atan2(x_sol(:,5), x_sol(:,4))-x_sol(:,3)), 'g', 'LineWidth', 1.5);
ylabel('\alpha (°)'); grid on; title('Angle d attaque');

% vitesse longi
subplot(4,1,3);
plot(t_sol, sqrt(x_sol(:,4).^2+x_sol(:,5).^2), 'b', 'LineWidth', 1.5);
ylabel('v (m/s)'); grid on; title('Vitesse longitudinale');

% vitesse de lacet
subplot(4,1,4);
plot(t_sol, x_sol(:,6), 'm', 'LineWidth', 1.5);
ylabel('r (rad/s)'); grid on; title('Vitesse de lacet');
xlabel('t (s)');

%% Figure 3 : Commandes
figure('Name', 'Commandes');

subplot(2,1,1);
plot(t_sol, T_sol, 'b', 'LineWidth', 1.5);
ylabel('T (N)'); grid on; title('Poussée — échelon');

subplot(2,1,2);
plot(t_sol, rad2deg(theta_sol), 'r', 'LineWidth', 1.5);
ylabel('\theta (rad)'); grid on; title('Angle de gouverne — trapèze');
xlabel('t (s)');
