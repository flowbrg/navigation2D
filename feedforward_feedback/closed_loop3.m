function dXdt = closed_loop3(t, X, params, K, p_traj)
% Fonction de simulation de la boucle fermee
% Args: t   [s]     temps
%       X   [6x1]   tableau des etats
%       params  [9x1]   tableau des parametres
%       p_traj  [?x2]   tableau des coefficients des polynomes de la
%                       tajectoire ideale

sigma_r = trajectory(t,p_traj); % [x_ref, y_ref, phi_ref,...]
U     = controller(X, sigma_r, K, params);
dXdt  = dyn3(X, U, params);
end