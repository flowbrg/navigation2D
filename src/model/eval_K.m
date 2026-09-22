function [K_lon, K_lat] = eval_K(Klqr,uQuery)
% Interpolate LQR coefficients
% K = [K_lon, K_lat]
% Where K_lon = [K_e_l, K_du] and K_lat = [K_e_n, K_dphi, K_dv, K_dr]

% K_lat is chosen to be constant because the linear model does not make
% sense for u > ~1.2.
% K_e_l is constant because de_n is not a function of u.
% K_du is not constant and needs to be interpolated.

    K_du = interp1(Klqr.u0, Klqr.K_du, uQuery);
    K_lon = [Klqr.K_el, K_du];
    K_lat = Klqr.K_lat;

end