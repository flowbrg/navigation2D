function dXdt = closed_loop3(X, sigma_r, K, params)
% Fonction de simulation de la boucle fermee
% Args: X       [7x1]   tableau des etats
%       X_ref   [7x1]   tableau de la trajectoire de référence
%       params  [9x1]   tableau des parametres
%       p_traj  [?x2]   tableau des coefficients des polynomes de la
%                       tajectoire ideale

U = controller(X, sigma_r, K, params);
dXdt  = dyn3(X, U, params);
end