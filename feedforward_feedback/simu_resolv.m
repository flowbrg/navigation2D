%% Planification de trajectoire polynomiale — Bateau 2D
%  Chemin géométrique p(s) = [x_ref(s), y_ref(s)],  s ∈ [0,1]
%  Optimisation fmincon (SQP)
%  Critère : λ₁·∫||p''||²ds  +  λ₂·∫κ²ds
%
%  Sorties : x_ref, y_ref, + dérivées — entrée de la boucle fermée
%
%  Dépendances : init_params.m, scenario_random.m  (inchangés)

clear; clc; close all;

init_params
scenario_random

%% =========================================================
%% 1. Paramètres du planificateur
%% =========================================================

n       = 5;      % Degré polynomial (min 5 pour 4 BC positions + marge)
Nc      = 150;    % Points de collocation — contraintes obstacles
Ns      = 500;    % Résolution de la trajectoire de sortie

lambda1 = 1.0;    % Poids ∫||p''||² ds  (quadratique, régularité)
lambda2 = 0.05;   % Poids ∫κ² ds        (non linéaire, proxy poussée)

% Caps aux extrémités [rad] — à aligner avec le contrôleur
phi0   = pi/4;    % Cap initial (cohérent avec simu_dyn2 : x0(3) = pi/4)
phi_f  = pi/4;    % Cap final souhaité
impose_tangents = true;

% Normalisation de la tangente imposée :
%   scale = dist/2 donne une bonne approximation arc-longueur
dist_ref = norm(target_pos - start_pos);
scale    = dist_ref / 2;

%% ========================================================
% 2. Base monomiale et ses dérivées
% =========================================================
%
%   p(s)  = B(s)   * a,   B_k(s)   = s^{k-1}
%   p'(s) = Bd(s)  * a,   Bd_k(s)  = (k-1)·s^{k-2}
%   p''(s)= Bdd(s) * a,   Bdd_k(s) = (k-1)(k-2)·s^{k-3}
%
%   Toutes les fonctions acceptent s scalaire ou colonne (Ns×1).
%   Retour : matrice (numel(s) × n+1)

B   = @(s) s .^ (0:n);
Bd  = @(s) [ zeros(numel(s),1),   s.^(0:n-1)  .* (1:n)           ];
Bdd = @(s) [ zeros(numel(s),2),   s.^(0:n-2)  .* ((2:n).*(1:n-1)) ];

%% ========================================================
% 3. Matrice de Gram — terme quadratique ∫₀¹ ||p''||² ds
% =========================================================
%
%   J_quad = a_x'·Q_dd·a_x + a_y'·Q_dd·a_y = w'·H_quad·w
%
%   [Q_dd]_{ij} = (i-1)(i-2)(j-1)(j-2) / (i+j-5),  i,j ≥ 3
%   (indices MATLAB 1-based, degré du monôme = indice - 1)
%   Démonstration : ∫₀¹ (i-1)(i-2)s^{i-3} · (j-1)(j-2)s^{j-3} ds
%                 = (i-1)(i-2)(j-1)(j-2) / (i+j-5)

Q_dd = zeros(n+1);
for i = 3:n+1
    for j = 3:n+1
        Q_dd(i,j) = (i-1)*(i-2)*(j-1)*(j-2) / (i+j-5);
    end
end

% H_quad est symétrique définie positive (les entrées non nulles forment
% un bloc plein dans le coin inférieur-droit de Q_dd)
H_quad = lambda1 * blkdiag(Q_dd, Q_dd);   % (2(n+1) × 2(n+1))

%% ========================================================
% 4. Contraintes d'égalité — conditions aux limites
% =========================================================
%
%   w = [a_x (n+1×1) ; a_y (n+1×1)]
%   Aeq * w = beq

xs = start_pos(1);  ys = start_pos(2);
xt = target_pos(1); yt = target_pos(2);

B0  = B(0);    B1  = B(1);
Bd0 = Bd(0);   Bd1 = Bd(1);
O   = zeros(1, n+1);   % bloc nul pour séparer les axes x et y

if impose_tangents
    % 8 équations : position + tangente aux deux extrémités
    % Tangente imposée : p'(0) = scale·[cos φ₀, sin φ₀]ᵀ
    %                   p'(1) = scale·[cos φ_f, sin φ_f]ᵀ
    Aeq = [ B0,   O  ;    % x(0) = xs
            B1,   O  ;    % x(1) = xt
            O,    B0 ;    % y(0) = ys
            O,    B1 ;    % y(1) = yt
            Bd0,  O  ;    % x'(0) = scale·cos φ₀
            O,    Bd0;    % y'(0) = scale·sin φ₀
            Bd1,  O  ;    % x'(1) = scale·cos φ_f
            O,    Bd1 ];  % y'(1) = scale·sin φ_f

    beq = [ xs; xt; ys; yt;
            scale*cos(phi0); scale*sin(phi0);
            scale*cos(phi_f); scale*sin(phi_f) ];
else
    % 4 équations : positions uniquement
    Aeq = [ B0, O; B1, O; O, B0; O, B1 ];
    beq = [ xs; xt; ys; yt ];
end

%% ========================================================
% 5. Contraintes non linéaires — évitement d'obstacles
% =========================================================
%
%   Condition : dist²(p(sⱼ), obs_i) ≥ r_safe²
%   Forme fmincon : c(w) = r_safe² - dist² ≤ 0
%
%   Nc points × n_obs obstacles = Nc·n_obs inégalités
%   Implémentation vectorisée : dx (Nc × n_obs) par broadcasting

s_col   = linspace(0, 1, Nc)';
r_safe2 = (R_obs + ecart)^2;

nonlcon = @(w) obs_avoid(w, n, s_col, B, obs, n_obs, r_safe2);

%% ========================================================
% 6. Critère d'optimisation (handle)
% =========================================================
s_int = linspace(0, 1, 300)';   % grille d'intégration (règle des trapèzes)

obj = @(w) eval_cost(w, n, H_quad, lambda2, s_int, Bd, Bdd);

%% ========================================================
% 7. Point initial — projection sur les contraintes d'égalité
% =========================================================
%
%   1) Ajustement LS du polynôme sur la ligne droite
%   2) Projection orthogonale sur {w : Aeq·w = beq}
%      → min ||w - w_line||² s.t. Aeq·w = beq
%      → solution : w₀ = w_line + Aeq'·(Aeq·Aeq')⁻¹·(beq - Aeq·w_line)

s_fit = linspace(0, 1, 50)';
Bmat  = B(s_fit);
ax0   = Bmat \ (xs + s_fit*(xt - xs));
ay0   = Bmat \ (ys + s_fit*(yt - ys));
w_line = [ax0; ay0];

% Projection (garantit que le point initial satisfait les BC)
w0 = w_line + Aeq' * ((Aeq*Aeq') \ (beq - Aeq*w_line));

%% ========================================================
% 8. Résolution fmincon — SQP
% =========================================================
%
%   SQP préféré à interior-point ici :
%   - Gère les contraintes actives proprement (obstacles tangents)
%   - Convergence plus rapide sur des NLP de cette taille
%   - Gradient numérique (différences finies, suffisant pour n=5)

opts_fmc = optimoptions('fmincon',             ...
    'Algorithm',               'sqp',          ...
    'Display',                 'iter',         ...
    'MaxIterations',           500,            ...
    'MaxFunctionEvaluations',  20000,          ...
    'OptimalityTolerance',     1e-6,           ...
    'ConstraintTolerance',     1e-4);

[w_sol, J_sol, exitflag, output] = fmincon(obj, w0,   ...
    [], [], Aeq, beq, [], [], nonlcon, opts_fmc);

fprintf('\n=== Planificateur — Résultats ===\n');
fprintf('Exitflag : %d  (%s)\n', exitflag, exit_msg(exitflag));
fprintf('Critère J : %.6f\n',    J_sol);
fprintf('Itérations : %d\n',     output.iterations);

if exitflag <= 0
    warning('Convergence non garantie — vérifier la trajectoire visuellement.');
end

%% ========================================================
% 9. Extraction de la trajectoire de référence
% =========================================================
s_ref  = linspace(0, 1, Ns)';
ax_sol = w_sol(1:n+1);
ay_sol = w_sol(n+2:end);

Bs    = B(s_ref);    Bds   = Bd(s_ref);    Bdds  = Bdd(s_ref);

x_ref   = Bs   * ax_sol;   xd_ref  = Bds  * ax_sol;   xdd_ref = Bdds * ax_sol;
y_ref   = Bs   * ay_sol;   yd_ref  = Bds  * ay_sol;   ydd_ref = Bdds * ay_sol;

% Courbure géométrique κ(s) = (x'y'' - y'x'') / ||p'||³
denom_k   = (xd_ref.^2 + yd_ref.^2).^(3/2) + 1e-10;   % régularisation
kappa_ref = (xd_ref .* ydd_ref - yd_ref .* xdd_ref) ./ denom_k;

% Cap tangentiel φ_ref(s) — utile pour l'initialisation du feedforward
phi_ref = atan2(yd_ref, xd_ref);

%% ========================================================
% 10. Vérification post-optimisation
% =========================================================
min_dist = inf;
for i = 1:n_obs
    d = min(sqrt((x_ref - obs(i,1)).^2 + (y_ref - obs(i,2)).^2));
    min_dist = min(min_dist, d);
end
long_approx = trapz(s_ref, sqrt(xd_ref.^2 + yd_ref.^2));

fprintf('Longueur approx.             : %.2f m\n',   long_approx);
fprintf('Courbure max |κ|             : %.4f m⁻¹\n', max(abs(kappa_ref)));
fprintf('Distance min aux obstacles   : %.4f m\n',   min_dist);
fprintf('Marge requise (R+écart)      : %.4f m\n',   sqrt(r_safe2));

%% ========================================================
% 11. Visualisation
% =========================================================
theta_c = linspace(0, 2*pi, 100);

figure(1); hold on; axis equal; grid on;
xlim([-1, 52]); ylim([-1, 52]);
xlabel('x (m)'); ylabel('y (m)');
title(sprintf('Trajectoire polynomiale  n=%d  |  J=%.3f', n, J_sol));

for i = 1:n_obs
    fill(obs(i,1) + (R_obs+ecart)*cos(theta_c), ...
         obs(i,2) + (R_obs+ecart)*sin(theta_c), ...
         [1.0 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    fill(obs(i,1) + R_obs*cos(theta_c), ...
         obs(i,2) + R_obs*sin(theta_c), ...
         [0.8 0.2 0.2], 'EdgeColor', 'k');
end

plot(x_ref, y_ref, 'b-',  'LineWidth', 2.5, 'DisplayName', 'Trajectoire ref');
plot(xs, ys, 'gs', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', 'Départ');
plot(xt, yt, 'p',  'MarkerSize', 16, 'MarkerFaceColor', 'y', ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Cible');
quiver(xs, ys, cos(phi0), sin(phi0), 0.8, 'g', 'LineWidth', 2, ...
    'DisplayName', 'Cap \phi_0');
legend('Location', 'northwest');

figure(2);
subplot(3,1,1);
plot(s_ref, x_ref, 'b', s_ref, y_ref, 'r', 'LineWidth', 1.5);
legend('x_{ref}', 'y_{ref}'); ylabel('Position (m)'); grid on;
title('Composantes de la trajectoire de référence');

subplot(3,1,2);
plot(s_ref, rad2deg(phi_ref), 'm', 'LineWidth', 1.5);
ylabel('\phi_{ref} (°)'); grid on; title('Cap tangentiel');

subplot(3,1,3);
plot(s_ref, kappa_ref, 'k', 'LineWidth', 1.5);
yline(0, '--', 'Color', [0.5 0.5 0.5]);
ylabel('\kappa (m^{-1})'); xlabel('s'); grid on;
title('Courbure géométrique');

%% ========================================================
% Fonctions locales
% =========================================================

function J = eval_cost(w, n, H_quad, lambda2, s, Bd, Bdd)
%EVAL_COST  Critère J = λ₁·∫||p''||²ds + λ₂·∫κ²ds
%   Terme 1 : analytique via H_quad  (quadratique en w)
%   Terme 2 : numérique par trapèzes (non linéaire en w)
    ax = w(1:n+1);
    ay = w(n+2:end);

    % --- Terme quadratique (analytique) ---
    J = w' * H_quad * w;

    % --- Terme de courbure (numérique) ---
    Bd_s  = Bd(s);    Bdd_s = Bdd(s);
    xd    = Bd_s  * ax;   yd    = Bd_s  * ay;
    xdd   = Bdd_s * ax;   ydd   = Bdd_s * ay;

    num   = xd .* ydd - yd .* xdd;                       % numérateur κ
    denom = (xd.^2 + yd.^2).^(3/2) + 1e-10;             % ||p'||³ + ε

    J = J + lambda2 * trapz(s, (num ./ denom).^2);
end

% ---------------------------------------------------------

function [c, ceq] = obs_avoid(w, n, s_col, B, obs, n_obs, r_safe2)
%OBS_AVOID  Contraintes d'évitement d'obstacles
%   c(k) = r_safe² - dist²(p(sⱼ), obs_i) ≤ 0
%   Vectorisé : dx, dy  de taille (Nc × n_obs)
    ax = w(1:n+1);
    ay = w(n+2:end);

    Bs = B(s_col);
    xc = Bs * ax;   % (Nc × 1)
    yc = Bs * ay;

    % Broadcasting : obs(:,1)' est (1 × n_obs) → dx est (Nc × n_obs)
    dx    = xc - obs(:,1)';
    dy    = yc - obs(:,2)';
    dist2 = dx.^2 + dy.^2;       % (Nc × n_obs)

    c   = r_safe2 - dist2(:);    % (Nc·n_obs × 1), doit être ≤ 0
    ceq = [];
end

% ---------------------------------------------------------

function msg = exit_msg(flag)
%EXIT_MSG  Message lisible pour exitflag fmincon
    switch flag
        case  1;   msg = 'Convergence (critère + contraintes OK)';
        case  2;   msg = 'Convergence (pas trop petit)';
        case  0;   msg = 'Itérations max atteintes';
        case -1;   msg = 'Arrêt par fonction output';
        case -2;   msg = 'Problème infaisable';
        otherwise; msg = sprintf('Code %d', flag);
    end
end