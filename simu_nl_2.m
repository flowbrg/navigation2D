%% Simulation navigation 2D basique
% Ca derape

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
theta_max   = pi/12;
t1 = 2; t2 = 4; t3 = 6; t4 = 8;

% Définition des commandes
T_cmd       = @(t) T0;  % Echelon unitaire
theta_cmd = @(t) ...
    (t >= t1 & t < t2)  .* (theta_max * (t - t1)/(t2 - t1)) + ...
    (t >= t2 & t < t3)  .* theta_max + ...
    (t >= t3 & t <= t4) .* (theta_max * (t4 - t)/(t4 - t3));
%theta_cmd = @(t) ...
%    (t >= t1 & t < t4)  .* theta_max;

%% Conditions initiales
% [x, y, phi, vx, vy, r)
x0 = [0; 0; pi/4; 0; 0; 0];    % Cap initial de 45°, vers la cible

% Intégration ode45, explicit runge-kutta
t_span          = [0,20];

ode_fun         = @(t,x) dynamics(t, x, T_cmd, theta_cmd, m, I, f, Lg, g, rho, S);
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

%% Fonction
function dxdt = dynamics(t, x, T_cmd, theta_cmd, m, I, f, Lg, g, rho, S)
    phi = x(3);
    vx  = x(4);
    vy  = x(5);
    r   = x(6);
    
    T     = T_cmd(t);
    theta = theta_cmd(t);
    
    % Norme et direction de la vitesse
    u   = sqrt(vx^2 + vy^2);
    psi = atan2(vy, vx);    % alpha + phi
    
    % Angles
    alpha = psi - phi;
    beta  = alpha - theta;
    
    % Force de la gouverne
    Fg = rho*S*sin(2*beta)*u^2;
    
    % Base du repere de Frenet
    eT = [cos(psi); sin(psi)];
    eN = [-sin(psi); cos(psi)];
    
    % Forces
    FT = T*[cos(phi); sin(phi)];
    Fg_vec = Fg*[sin(phi-theta); -cos(phi-theta)];
    Ff = -f*u^2*eT;
    
    % Acceleration (calcul avec des tableaux)
    a = (FT + Fg_vec + Ff)/m;
    
    % Equations d'etat
    dx  = vx;
    dy  = vy;
    dphi = r;
    dvx = a(1);
    dvy = a(2);
    dr  = (Fg*Lg*cos(theta) - g*r)/I;
    
    dxdt = [dx; dy; dphi; dvx; dvy; dr];
end