function traj_fun = make_trajectory(p_x, p_y, Tf)
p_dx   = polyder(p_x);      p_dy   = polyder(p_y);
p_ddx  = polyder(p_dx);     p_ddy  = polyder(p_dy);
p_dddx = polyder(p_ddx);    p_dddy = polyder(p_ddy);

traj_fun = @(t) eval_traj(t, p_x, p_y, p_dx, p_dy, ...
                            p_ddx, p_ddy, p_dddx, p_dddy, Tf);
end