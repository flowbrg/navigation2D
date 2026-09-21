function [U, e_n] = controller(X, sigma_r, K, params)

T_max     = params(10);
theta_max = params(11);

% Feedforward
[U_ff, X_ref] = feedforward3(sigma_r, params);

%Feedback
[dT, dtheta, e_n]  = feedback3(X, X_ref, K, params);

T_ff = U_ff(1);
theta_ff = U_ff(2);

T = max(0, min(T_max, T_ff+dT));
theta = max(-theta_max, min(theta_max, theta_ff+dtheta));

U = [T, theta];
end