function sigma_r = trajectory(t, p_traj)

% Extrait les coefficients des polynomes pour x et y
p_x = p_traj(1,:);
p_y = p_traj(2,:);

% Calcul des derivees successives
p_dx   = polyder(p_x);
p_dy   = polyder(p_y);
p_ddx  = polyder(p_dx);
p_ddy  = polyder(p_dy);
p_dddx = polyder(p_ddx);
p_dddy = polyder(p_ddy);

% Evaluation des polynomes en t
x_r = polyval(p_x, t);
y_r = polyval(p_y, t);
dx    = polyval(p_dx, t);
dy    = polyval(p_dy, t);
ddx   = polyval(p_ddx, t);
ddy   = polyval(p_ddy, t);
dddx  = polyval(p_dddx, t);
dddy  = polyval(p_dddy, t);

sigma_r = [x_r, y_r, dx, dy, ddx, ddy, dddx, dddy];
end