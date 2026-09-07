%% Simulation navigation 2D basique
% Modele d'etat avec vitesses cartesiennes dans le repere du bateau
% X = [x, y, phi, u, v, r]

clear; clc; close all;

init_params
init_lqr


%% Niveau 1 - Trajectoire en ligne droite
Tf = 20; u0 = 3;
p_x = [u0, 0]; % x(t) = u0*t
p_y = [0, 0];  % y(t) = 0
X0 = [0; 0; 0; u0; 0; 0];

traj_fun = make_trajectory(p_x, p_y, Tf); % [x_ref, y_ref, dx_ref,... ]
ode_fun  = @(t,X) closed_loop3(X, traj_fun(t), K, params);
options  = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t_sol, X_sol] = ode45(ode_fun, [0, Tf], X0, options);

% Figure 1 : Trajectoire
figure('Name', 'T1 - Trajectoire'); hold on; axis equal; grid on;
%xlim([-15, 15]); ylim([-15, 6]);

xlabel('x (m)'); ylabel('y (m)');
plot(X_sol(:,1), X_sol(:,2), 'b-', 'LineWidth', 2, 'DisplayName', 'Trajectoire');
plot(polyval(p_x, t_sol), polyval(p_y, t_sol), 'r-', 'LineWidth', 2, 'DisplayName', 'Référence');

% Orientation de la vitesse
n_arrows = 15;  % Nombre de flèches
idx = round(linspace(1, length(t_sol), n_arrows));
alpha_sol = atan2(X_sol(:,5), X_sol(:,4));
V_sol   = sqrt(X_sol(:,5).^2 + X_sol(:,4).^2);
quiver(X_sol(idx,1), X_sol(idx, 2), V_sol(idx).*cos(alpha_sol(idx)+X_sol(idx, 3)),  V_sol(idx).*sin(alpha_sol(idx)+X_sol(idx, 3)), 0, 'k', 'LineWidth', 1.2, 'DisplayName', 'Vitesse \alpha+\phi');
% Cap
quiver(X_sol(idx,1), X_sol(idx, 2), cos(X_sol(idx,3)),  sin(X_sol(idx,3)), 0, 'r', 'LineWidth', 1.2, 'DisplayName', 'Cap \phi');


legend('Location', 'northwest');

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

%% Niveau 2 - Virage à courbe constante
R = 30; Tf = 100; u0 = 5;
t_fit = linspace(0, Tf, 50);
x_fit = R*sin(t_fit/R);
y_fit = R*(1-cos(t_fit/R));
p_x = polyfit(t_fit, x_fit, 5);
p_y = polyfit(t_fit, y_fit, 5);
X0 = [0; 0; 0; u0; 0; 0];

traj_fun = make_trajectory(p_x, p_y, Tf);
ode_fun  = @(t,X) closed_loop3(X, traj_fun(t), K, params);
options  = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t_sol, X_sol] = ode45(ode_fun, [0, Tf], X0, options);

% Figure T2 : Trajectoire
figure('Name', 'T2 - Trajectoire'); hold on; axis equal; grid on;
%xlim([-15, 15]); ylim([-15, 6]);

xlabel('x (m)'); ylabel('y (m)');
plot(X_sol(:,1), X_sol(:,2), 'b-', 'LineWidth', 2, 'DisplayName', 'Trajectoire');
plot(polyval(p_x, t_sol), polyval(p_y, t_sol), 'r-', 'LineWidth', 2, 'DisplayName', 'Référence');

% Orientation de la vitesse
n_arrows = 15;  % Nombre de flèches
idx = round(linspace(1, length(t_sol), n_arrows));
alpha_sol = atan2(X_sol(:,5), X_sol(:,4));
V_sol   = sqrt(X_sol(:,5).^2 + X_sol(:,4).^2);
quiver(X_sol(idx,1), X_sol(idx, 2), V_sol(idx).*cos(alpha_sol(idx)+X_sol(idx, 3)),  V_sol(idx).*sin(alpha_sol(idx)+X_sol(idx, 3)), 0, 'k', 'LineWidth', 1.2, 'DisplayName', 'Vitesse \alpha+\phi');
% Cap
quiver(X_sol(idx,1), X_sol(idx, 2), cos(X_sol(idx,3)),  sin(X_sol(idx,3)), 0, 'r', 'LineWidth', 1.2, 'DisplayName', 'Cap \phi');


legend('Location', 'southwest');

% Figure T2 : Etats
figure('Name', 'T2 - États');

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

%% Niveau 3 - Virage à courbe constante
X0 = [0; 0.5; deg2rad(10); u0*0.8; 0; 0];

traj_fun = make_trajectory(p_x, p_y, Tf);
ode_fun  = @(t,X) closed_loop3(X, traj_fun(t), K, params);
options  = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t_sol, X_sol] = ode45(ode_fun, [0, Tf], X0, options);

% Figure T2 : Trajectoire
figure('Name', 'T3 - Trajectoire'); hold on; axis equal; grid on;
%xlim([-15, 15]); ylim([-15, 6]);

xlabel('x (m)'); ylabel('y (m)');
plot(X_sol(:,1), X_sol(:,2), 'b-', 'LineWidth', 2, 'DisplayName', 'Trajectoire');
plot(polyval(p_x, t_sol), polyval(p_y, t_sol), 'r-', 'LineWidth', 2, 'DisplayName', 'Référence');

% Orientation de la vitesse
n_arrows = 15;  % Nombre de flèches
idx = round(linspace(1, length(t_sol), n_arrows));
alpha_sol = atan2(X_sol(:,5), X_sol(:,4));
V_sol   = sqrt(X_sol(:,5).^2 + X_sol(:,4).^2);
quiver(X_sol(idx,1), X_sol(idx, 2), V_sol(idx).*cos(alpha_sol(idx)+X_sol(idx, 3)),  V_sol(idx).*sin(alpha_sol(idx)+X_sol(idx, 3)), 0, 'k', 'LineWidth', 1.2, 'DisplayName', 'Vitesse \alpha+\phi');
% Cap
quiver(X_sol(idx,1), X_sol(idx, 2), cos(X_sol(idx,3)),  sin(X_sol(idx,3)), 0, 'r', 'LineWidth', 1.2, 'DisplayName', 'Cap \phi');


legend('Location', 'southwest');

% Figure T2 : Etats
figure('Name', 'T3 - États');

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