%% Simulation navigation 2D basique
% La direction est gérée par la directino de l'hélice, il n'y a pas de
% gouverne

clear; clc; close all;
init

%% Signaux de commande
% Paramètres des commandes
% Poussée, échelon unitaire
T0  = 500; 

% Couple, trapèze
% Montée linéaire de 0 à theta_max entre t1 et t2
% Palier à theta_max entre t2 et t3
% Descente linéaire de theta_max à 0 entre t3 et t4
theta_max   = pi/6;
t1 = 2; t2 = 4; t3 = 6; t4 = 8;

% Définition des commandes
T_cmd       = @(t) T0;  % Echelon unitaire
theta_cmd = @(t) ...
    (t >= t1 & t < t2)  .* (theta_max * (t - t1)/(t2 - t1)) + ...
    (t >= t2 & t < t3)  .* theta_max + ...
    (t >= t3 & t <= t4) .* (theta_max * (t4 - t)/(t4 - t3));

%% Conditions initiales
% [x, y, phi, alpha, u, r)
x0 = [0; 0; pi/4; 0; 1e-2; 0];    % Cap initial de 45°, vers la cible

% Intégration ode45, explicit runge-kutta
t_span          = [0,20];

ode_fun         = @(t,x) dynamics(t, x, T_cmd, theta_cmd, m, I, f, g, Lg);
options         = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t_sol, x_sol]  = ode45(ode_fun, t_span, x0, options);

T_sol           = arrayfun(T_cmd, t_sol);
theta_sol       = arrayfun(theta_cmd, t_sol);

%% Figure 1 : Trajectoire
figure(1); hold on; axis equal; grid on;
%xlim([-15, 15]); ylim([-15, 6]);
xlabel('x (m)'); ylabel('y (m)');
plot(x_sol(:,1), x_sol(:,2), 'b-', 'LineWidth', 2, 'DisplayName', 'Trajectoire');

% Orientation du bateau
n_arrows = 15;  % Nombre de flèches
idx = round(linspace(1, length(t_sol), n_arrows));
quiver(x_sol(idx,1), x_sol(idx, 2), 2*cos(x_sol(idx,3)+x_sol(idx,4)),  2*sin(x_sol(idx,3)+x_sol(idx,4)), 0, 'k', 'LineWidth', 1.2, 'DisplayName', 'Cap \phi');

legend('Location', 'northwest');

%% Figure 2 : Etats
figure('Name', 'États');

% cap
subplot(4,1,1);
plot(t_sol, rad2deg(x_sol(:,3)), 'r', 'LineWidth', 1.5);
ylabel('\phi (°)'); grid on; title('Cap');

% angle d'incidence
subplot(4,1,2);
plot(t_sol, rad2deg(x_sol(:,4)), 'g', 'LineWidth', 1.5);
ylabel('\alpha (°)'); grid on; title('Angle d attaque');

% vitesse longi
subplot(4,1,3);
plot(t_sol, x_sol(:,5), 'b', 'LineWidth', 1.5);
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
plot(t_sol, theta_sol, 'r', 'LineWidth', 1.5);
ylabel('\theta (rad)'); grid on; title('Couple de gouverne — trapèze');
xlabel('t (s)');

%% Fonction dynamique
function dxdt = dynamics(t, x, T_cmd, theta_cmd, m, I, f, g, Lg)
    phi = x(3);
    alpha = x(4);
    u   = x(5);
    r   = x(6);

    T       = T_cmd(t);
    theta   = theta_cmd(t);
    
    psi = alpha + phi;
    beta = alpha-theta;

    dxdt    = zeros(5,1);
    dxdt(1) = u*cos(psi);
    dxdt(2) = u*sin(psi);
    dxdt(3) = r;
    dxdt(4) = -T*sin(beta)/(m*u)-r;
    dxdt(5) = (T*cos(beta)-f*u^2)/m;
    dxdt(6) = (-Lg*T*sin(theta)-g*r)/I;   % Hyp: theta << 1, sin(theta)=theta
end