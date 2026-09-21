%% plan_poly.m — Planificateur polynomial — sortie coefficients temporels
%  x_eq(t) = polyval(flip(cx), t),  t ∈ [0, Tf_plan]
%  y_eq(t) = polyval(flip(cy), t)
%
%  Dépendances : init_params.m, scenario_random.m (ou scenario.m)
%  Sorties exportées dans le workspace :
%    cx, cy      — vecteurs (n+1 × 1), coefficients de t^0 … t^n
%    Tf_plan     — durée de la trajectoire [s]
%    traj_ref    — struct de commodité pour le correcteur

init_params
scenario_random

%% =========================================================
%% 1. Paramètres du planificateur
%% =========================================================

n       = 7;        % Degré polynomial (≥ 5 pour 4 BC + marge)
Nc      = 150;      % Points de collocation pour les contraintes obstacles
lambda1 = 1.0;      % Poids ∫||p''||² ds  (régularité)
lambda2 = 0.05;     % Poids ∫κ² ds        (proxy courbure)

phi0   = pi/4;      % Cap initial [rad]  — doit correspondre à x0(3)
phi_f  = pi/4;      % Cap final  [rad]

v_mean = 1.5;       % Vitesse moyenne estimée [m/s] — paramètre utilisateur

dist_ref = norm(target_pos - start_pos);
scale    = dist_ref / 2;   % normalisation des tangentes

%% =========================================================
%% 2. Bases polynomiales
%% =========================================================
%  B(s)   : (numel(s) × n+1),  B_k(s)   = s^{k-1}  (1-based)
%  Bd(s)  : dérivée première
%  Bdd(s) : dérivée seconde

B   = @(s) s .^ (0:n);
Bd  = @(s) [ zeros(numel(s),1),  s.^(0:n-1) .* (1:n)            ];
Bdd = @(s) [ zeros(numel(s),2),  s.^(0:n-2) .* ((2:n).*(1:n-1)) ];

%% =========================================================
%% 3. Matrice de Gram — terme quadratique ∫₀¹ ||p''||² ds
%% =========================================================
%  [Q_dd]_{ij} = (i-1)(i-2)(j-1)(j-2)/(i+j-5),  i,j ≥ 3

Q_dd = zeros(n+1);
for i = 3:n+1
    for j = 3:n+1
        Q_dd(i,j) = (i-1)*(i-2)*(j-1)*(j-2) / (i+j-5);
    end
end
H_quad = lambda1 * blkdiag(Q_dd, Q_dd);   % (2(n+1) × 2(n+1))

%% =========================================================
%% 4. Contraintes d'égalité — conditions aux limites
%% =========================================================
xs = start_pos(1);  ys = start_pos(2);
xt = target_pos(1); yt = target_pos(2);

B0  = B(0);  B1  = B(1);
Bd0 = Bd(0); Bd1 = Bd(1);
O   = zeros(1, n+1);

Aeq = [ B0,  O  ;    % x(0) = xs
        B1,  O  ;    % x(1) = xt
        O,   B0 ;    % y(0) = ys
        O,   B1 ;    % y(1) = yt
        Bd0, O  ;    % x'(0) = scale·cos φ₀
        O,   Bd0;    % y'(0) = scale·sin φ₀
        Bd1, O  ;    % x'(1) = scale·cos φ_f
        O,   Bd1 ];  % y'(1) = scale·sin φ_f

beq = [ xs; xt; ys; yt;
        scale*cos(phi0); scale*sin(phi0);
        scale*cos(phi_f); scale*sin(phi_f) ];

%% =========================================================
%% 5. Contraintes obstacles + critère
%% =========================================================
s_col   = linspace(0, 1, Nc)';
r_safe2 = (R_obs + ecart)^2;

s_int   = linspace(0, 1, 300)';
obj     = @(w) eval_cost(w, n, H_quad, lambda2, s_int, Bd, Bdd);
nonlcon = @(w) obs_avoid(w, n, s_col, B, obs, n_obs, r_safe2);

%% =========================================================
%% 6. Point initial — projection sur les contraintes d'égalité
%% =========================================================
s_fit  = linspace(0, 1, 50)';
Bmat   = B(s_fit);
ax0    = Bmat \ (xs + s_fit*(xt - xs));
ay0    = Bmat \ (ys + s_fit*(yt - ys));
w_line = [ax0; ay0];

% Projection orthogonale → satisfait les BC exactement
w0 = w_line + Aeq' * ((Aeq*Aeq') \ (beq - Aeq*w_line));

%% =========================================================
%% 7. Optimisation — SQP
%% =========================================================
opts_fmc = optimoptions('fmincon',            ...
    'Algorithm',              'sqp',          ...
    'Display',                'iter',         ...
    'MaxIterations',          500,            ...
    'MaxFunctionEvaluations', 20000,          ...
    'OptimalityTolerance',    1e-6,           ...
    'ConstraintTolerance',    1e-4);

[w_sol, J_sol, exitflag, output] = fmincon(obj, w0, ...
    [], [], Aeq, beq, [], [], nonlcon, opts_fmc);

fprintf('Exitflag : %d | J = %.4f | Iter = %d\n', ...
    exitflag, J_sol, output.iterations);
if exitflag <= 0
    warning('plan_poly: convergence non garantie.');
end

%% =========================================================
%% 8. Extraction et reparamétrisation temporelle
%% =========================================================
ax_sol = w_sol(1:n+1);       % coefficients en s, axe x
ay_sol = w_sol(n+2:end);     % coefficients en s, axe y

% Longueur d'arc approx. sur la solution (règle des trapèzes)
s_check = linspace(0, 1, 500)';
xd_chk  = Bd(s_check) * ax_sol;
yd_chk  = Bd(s_check) * ay_sol;
L_arc   = trapz(s_check, sqrt(xd_chk.^2 + yd_chk.^2));

% Durée planifiée
Tf_plan = L_arc / v_mean;    % [s]

% Reparamétrisation : s = t/Tf → c_k = a_k / Tf^k
k_vec = (0:n)';              % exposants
cx = ax_sol ./ (Tf_plan .^ k_vec);   % (n+1 × 1)
cy = ay_sol ./ (Tf_plan .^ k_vec);   % (n+1 × 1)

%% =========================================================
%% 9. Export — struct correcteur
%% =========================================================
traj_ref.cx      = cx;         % coefficients x(t) = sum cx_k * t^k
traj_ref.cy      = cy;         % coefficients y(t) = sum cy_k * t^k
traj_ref.Tf      = Tf_plan;    % durée [s]
traj_ref.L_arc   = L_arc;      % longueur arc [m]
traj_ref.n       = n;          % degré

fprintf('\n=== Trajectoire de référence ===\n');
fprintf('Longueur arc  : %.3f m\n',   L_arc);
fprintf('Tf planifié   : %.2f s  (v_mean = %.2f m/s)\n', Tf_plan, v_mean);

% Évaluation à la demande dans le correcteur :
%   x_eq = polyval(flip(cx), t)    pour t ∈ [0, Tf_plan]
%   xd_eq = polyval(flip(cx(2:end) .* (1:n)'), t)   (dérivée)

%% =========================================================
%% Fonctions locales
%% =========================================================

function J = eval_cost(w, n, H_quad, lambda2, s, Bd, Bdd)
    ax = w(1:n+1);   ay = w(n+2:end);
    J  = w' * H_quad * w;
    Bd_s = Bd(s);    Bdd_s = Bdd(s);
    xd   = Bd_s * ax;    yd   = Bd_s * ay;
    xdd  = Bdd_s * ax;   ydd  = Bdd_s * ay;
    num  = xd .* ydd - yd .* xdd;
    den  = (xd.^2 + yd.^2).^(3/2) + 1e-10;
    J    = J + lambda2 * trapz(s, (num ./ den).^2);
end

function [c, ceq] = obs_avoid(w, n, s_col, B, obs, n_obs, r_safe2)
    ax = w(1:n+1);   ay = w(n+2:end);
    Bs = B(s_col);
    xc = Bs * ax;    yc = Bs * ay;
    dx    = xc - obs(:,1)';
    dy    = yc - obs(:,2)';
    dist2 = dx.^2 + dy.^2;
    c   = r_safe2 - dist2(:);
    ceq = [];
end