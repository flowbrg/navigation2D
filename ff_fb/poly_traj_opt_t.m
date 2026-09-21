%% Planification de trajectoire polynomiale — Bateau 2D
%  Chemin géométrique p(s) = [x_ref(s), y_ref(s)],  s ∈ [0,1]
%  Optimisation fmincon (SQP)
%  Critère : λ₁·∫||p''||²ds  +  λ₂·∫κ²ds
%
%  Sorties : x_ref, y_ref, + dérivées — entrée de la boucle fermée
%
%  Dépendances : init_params.m, scenario_random.m  (inchangés)

%clear; clc; close all;

init_params
scenario

%% =========================================================
%% 1. Paramètres du planificateur
%% =========================================================

Tf   = 30;    % Durée de la manœuvre [s] — à ajuster
n    = 7;     % Degré polynomial (min 7 pour 4 conditions × 2 axes + dérivées)
              % feedforward3 requiert jusqu'à d³x/dt³ → min degré 5,
              % mais 7 recommandé pour avoir de la liberté d'évitement
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

B   = @(t) t .^ (0:n);                                         % (Nt x n+1)
Bd  = @(t) [zeros(numel(t),1), t.^(0:n-1) .* (1:n)          ]; % dx/dt
Bdd = @(t) [zeros(numel(t),2), t.^(0:n-2) .* ((2:n).*(1:n-1))];% ddx/dt2

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
        Q_dd(i,j) = (i-1)*(i-2)*(j-1)*(j-2) * Tf^(i+j-5) / (i+j-5);
    end
end
H_quad = lambda1 * blkdiag(Q_dd, Q_dd);


% H_quad est symétrique définie positive (les entrées non nulles forment
% un bloc plein dans le coin inférieur-droit de Q_dd)
H_quad = lambda1 * blkdiag(Q_dd, Q_dd);   % (2(n+1) × 2(n+1))

%% ========================================================
% 4. Contraintes d'égalité — conditions aux limites
% =========================================================
%
%  feedforward3 utilise dx/dt, dy/dt → imposer la vitesse initiale
%  cohérente avec u_ref = V = sqrt(dx2+dy2) et phi_ref = atan2(dy,dx)
%
%  Choix : vitesse initiale u0 selon le cap phi0

xs = start_pos(1);  ys = start_pos(2);
xt = target_pos(1); yt = target_pos(2);

u0     = 2.0;    % [m/s] vitesse de croisière initiale — à caler sur le LQR
u_f    = 2.0;    % [m/s] vitesse finale souhaitée

dx0 = u0 * cos(phi0);   dy0 = u0 * sin(phi0);
dxf = u_f * cos(phi_f); dyf = u_f * sin(phi_f);

B0  = B(0);    BTf = B(Tf);
Bd0 = Bd(0);   BdTf = Bd(Tf);
O   = zeros(1, n+1);

Aeq = [ B0,   O  ;
        BTf,  O  ;
        O,    B0 ;
        O,    BTf;
        Bd0,  O  ;
        O,    Bd0;
        BdTf, O  ;
        O,    BdTf];

beq = [xs; xt; ys; yt; dx0; dy0; dxf; dyf];

%% ========================================================
% 5. Contraintes non linéaires — évitement d'obstacles
% =========================================================
%
%   Condition : dist²(p(sⱼ), obs_i) ≥ r_safe²
%   Forme fmincon : c(w) = r_safe² - dist² ≤ 0
%
%   Nc points × n_obs obstacles = Nc·n_obs inégalités
%   Implémentation vectorisée : dx (Nc × n_obs) par broadcasting

t_col   = linspace(0, Tf, Nc)';
r_safe2 = (R_obs + ecart)^2;
nonlcon = @(w) obs_avoid(w, n, t_col, B, obs, n_obs, r_safe2);

%% ========================================================
% 6. Critère d'optimisation (handle)
% =========================================================
t_int = linspace(0, Tf, 300)';  % grille d'intégration (règle des trapèzes)

obj   = @(w) eval_cost(w, n, H_quad, lambda2, t_int, Bd, Bdd);

%% ========================================================
% 7. Point initial — projection sur les contraintes d'égalité
% =========================================================
%
%   1) Ajustement LS du polynôme sur la ligne droite
%   2) Projection orthogonale sur {w : Aeq·w = beq}
%      → min ||w - w_line||² s.t. Aeq·w = beq
%      → solution : w₀ = w_line + Aeq'·(Aeq·Aeq')⁻¹·(beq - Aeq·w_line)

t_fit = linspace(0, Tf, 50)';
Bmat  = B(t_fit);
ax0   = Bmat \ (xs + (t_fit/Tf)*(xt - xs));
ay0   = Bmat \ (ys + (t_fit/Tf)*(yt - ys));
w_line = [ax0; ay0];
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
ax_sol = w_sol(1:n+1);
ay_sol = w_sol(n+2:end);

% flip : croissant → décroissant (convention polyval/polyder MATLAB)
p_x = flip(ax_sol)';   % (1 × n+1), ordre décroissant
p_y = flip(ay_sol)';

p_traj = [p_x; p_y];  % (2 × n+1) — format attendu par trajectory.m

% Vérification
t_ref  = linspace(0, Tf, Ns)';
x_ref  = polyval(p_x, t_ref);
y_ref  = polyval(p_y, t_ref);

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