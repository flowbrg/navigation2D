function dxdt = dyn(t, x, T_cmd, theta_cmd, m, I, f, Lg, g, rho, S)
    phi     = x(3);
    alpha   = x(4);
    u       = x(5);
    r       = x(6);

    T       = T_cmd(t);
    theta   = theta_cmd(t);

    beta = alpha - theta;   % direction réelle de la vitesse
    psi = alpha + phi;
    
    Fg = rho*S*sin(2*beta)*u^2; % Force de la gouverne
    
    dx      = u*cos(psi); % x
    dy      = u*sin(psi); % y
    dphi    = r; % phi
    dalpha  = ((-T*sin(alpha) - Fg*cos(beta))/(m*u))-r; % alpha
    du      = (T*cos(alpha) - f*u^2 - Fg*sin(beta))/m; % u
    dr      = (Fg*Lg*cos(theta)-g*r)/I; % r

    dxdt = [dx; dy; dphi; dalpha; du; dr];
end