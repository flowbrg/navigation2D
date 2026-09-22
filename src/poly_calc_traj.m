%% Planificateur polynomial — Bateau 2D
% Chemin géométrique p(t) = [x(t), y(t)], polynôme de degré n
% Critère : k1*int(norm(ddp/dt2)^2) + k2*int(K^2)
% Solveur  : fmincon / SQP

%clear; clc; close all;


%% 1. Paramètres
Tf      = 30;       % Durée [s]
n       = 5;        % Degré polynomial
Nc      = 150;      % Points pour contraintes obstacles
Ns      = 500;      % Résolution trajectoire de sortie
lambda1 = 1.0;      % Poids régularité
lambda2 = 0.05;     % Poids courbure
phi0    = pi/4;     % Cap initial [rad]
phi_f   = pi/4;     % Cap final   [rad]
u0      = 2.0;      % Vitesse initiale [m/s]
u_f     = 2.0;      % Vitesse finale   [m/s]
r_safe  = scen.r_obs + scen.ecart;

xs = scen.start_pos(1);  ys = scen.start_pos(2);
xt = scen.target_pos(1); yt = scen.target_pos(2);

%% 2. Base monomiale et dérivées
% B(t)   = [1, t, t^2, …, t^n]        →  p(t)   = B(t)·a
% Bd(t)  = [0, 1, 2t, …, n·t^(n-1)]   →  dp(t)  = Bd(t)·a
% Bdd(t) = [0, 0, 2,  …, n(n-1)t^(n-2)]→  ddp(t)  = Bdd(t)·a
B   = @(t) t .^ (0:n);
Bd  = @(t) [zeros(numel(t),1),  t.^(0:n-1) .* (1:n)           ];
Bdd = @(t) [zeros(numel(t),2),  t.^(0:n-2) .* ((2:n).*(1:n-1))];

%% 3. Matrice de Gram Q_dd
% [Q_dd]_{ij} = (i-1)(i-2)(j-1)(j-2)·Tf^{i+j-5} / (i+j-5),  i,j >= 3
Q_dd = zeros(n+1);
for i = 3:n+1
    for j = 3:n+1
        Q_dd(i,j) = (i-1)*(i-2)*(j-1)*(j-2) * Tf^(i+j-5) / (i+j-5);
    end
end
H_quad = lambda1 * blkdiag(Q_dd, Q_dd);   % (2(n+1) × 2(n+1))

%% 4. Contraintes d'égalité — conditions aux limites
dx0 = u0  * cos(phi0);   dy0 = u0  * sin(phi0);
dxf = u_f * cos(phi_f);  dyf = u_f * sin(phi_f);

B0 = B(0); BTf = B(Tf); Bd0 = Bd(0); BdTf = Bd(Tf);
O  = zeros(1, n+1);

Aeq = [B0, O; BTf, O; O, B0; O, BTf; Bd0, O; O, Bd0; BdTf, O; O, BdTf];
beq = [xs; xt; ys; yt; dx0; dy0; dxf; dyf];

%% 5. Contraintes obstacles et critère
t_col   = linspace(0, Tf, Nc)';
nonlcon = @(w) obs_avoid(w, n, t_col, B, scen.obs, r_safe^2);

t_int = linspace(0, Tf, 300)';
obj   = @(w) eval_cost(w, n, H_quad, lambda2, t_int, Bd, Bdd);

%% 6. Point initial — projection LS sur {Aeq*w = beq}
t_fit  = linspace(0, Tf, 50)';
Bmat   = B(t_fit);
ax0    = Bmat \ (xs + (t_fit/Tf)*(xt - xs));
ay0    = Bmat \ (ys + (t_fit/Tf)*(yt - ys));
w_line = [ax0; ay0];
w0     = w_line + Aeq' * ((Aeq*Aeq') \ (beq - Aeq*w_line));

%% 7. Optimisation fmincon / SQP
opts = optimoptions('fmincon',            ...
    'Algorithm',              'sqp',      ...
    'Display',                'iter',     ...
    'MaxIterations',          500,        ...
    'MaxFunctionEvaluations', 20000,      ...
    'OptimalityTolerance',    1e-6,       ...
    'ConstraintTolerance',    1e-4);

[w_sol, J_sol, exitflag, output] = fmincon(obj, w0, [], [], Aeq, beq, [], [], nonlcon, opts);
fprintf('Exitflag %d | J = %.4f | Iter = %d\n', exitflag, J_sol, output.iterations);

%% 8. Extraction de la référence
ax_sol = w_sol(1:n+1);
ay_sol = w_sol(n+2:end);
p_x    = flip(ax_sol)';
p_y    = flip(ay_sol)';

t_ref = linspace(0, Tf, Ns)';
x_ref = polyval(p_x, t_ref);
y_ref = polyval(p_y, t_ref);

%% Fonctions locales

function J = eval_cost(w, n, H_quad, lambda2, s, Bd, Bdd)
    ax = w(1:n+1);   ay = w(n+2:end);
    J  = w' * H_quad * w;                         % terme quadratique

    Bds = Bd(s); Bdds = Bdd(s);
    xd  = Bds*ax;  yd  = Bds*ay;
    xdd = Bdds*ax; ydd = Bdds*ay;

    kappa2 = ((xd.*ydd - yd.*xdd) ./ ((xd.^2 + yd.^2).^(3/2) + 1e-10)).^2;
    J = J + lambda2 * trapz(s, kappa2);
end

function [c, ceq] = obs_avoid(w, n, s_col, B, obs, r_safe2)
    ax = w(1:n+1);   ay = w(n+2:end);
    Bs = B(s_col);
    xc = Bs*ax;   yc = Bs*ay;

    dx    = xc - obs(:,1)';
    dy    = yc - obs(:,2)';
    c     = r_safe2 - (dx.^2 + dy.^2);
    c     = c(:);
    ceq   = [];
end