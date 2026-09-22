function U = controller(X, sigma_r, K, params)

% Feedforward
[U_ff, X_ref] = feedforward3(sigma_r, params);

%Feedback
[dT, dtheta] = feedback3(X, X_ref, K);

T_ff = U_ff(1);
theta_ff = U_ff(2);

T = max(0, min(params.T_max, T_ff+dT));
theta = max(-params.theta_max, min(params.theta_max, theta_ff+dtheta));

U = [T, theta];
end