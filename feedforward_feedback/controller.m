function U = controller(X, sigma_r, K, params)

% Saturations
T_max     = 700;
theta_max = pi/3;

% Feedforward
[U_ff, X_ref] = feedforward3(sigma_r, params);

[dT, dtheta]  = feedback3(X, X_ref, K);

T_ff = U_ff(1);
theta_ff = U_ff(2);

T = max(0, min(T_max, T_ff+dT));
theta = max(-theta_max, min(theta_max, theta_ff+dtheta));

U = [T, theta];
end