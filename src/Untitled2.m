clear all
close all
clc

options = optimoptions('fmincon','Display','iter');
samples = 20;
trials = 10;
n=3;
Y = [1+i 0 -1-i; 0 2+i -2-i; -1-i -2-i 3+2i];
x_ode_hat(:, 1) = 0.1 * ones([n*(n-1)/2, 1]);
Z_ode(:, :, 1) = eye(n*(n-1)/2);
lb = 0.9 * ones([n 1]);
ub = 1.1 * ones([n 1]);
Q = duplication_matrix(n) * transformation_matrix(n);
lambda = 0.5;
[~, X] = cordexch(n, samples);
X = X(:, 2:end).*0.1 + 1; 

for i = 2:samples
    for t = 1:trials
        x0(:, t) = lb + (ub-lb) .* rand([length(ub) 1]);
        f = @(x) - log(abs(det(inv(Z_ode(:, :, i-1)) + (kron(x, eye(n)).' * Q)' * (kron(x, eye(n)).' * Q))));
        [trial_res(:, t), trial_f(t)] = fmincon(f, x0(:, t), [], [], [], [], lb, ub, [], options);    
    end
    [min_f, min_index] = min(trial_f);
    V_real_ode(:, i) = trial_res(:, min_index);
    % V_real_ode(:, i) = X(i, :)';
    I_real_ode(:, i) = Y * V_real_ode(:, i) + normrnd(0, 0.1) + sqrt(-1)*normrnd(0, 0.1);
    
    a = kron(V_real_ode(:, i), eye(n)).' * Q;
    Z_ode(:, :, i) = (1/lambda) * (Z_ode(:, :, i-1) - Z_ode(:, :, i-1) * a' * inv(lambda * eye(n) + a * Z_ode(:, :, i-1) * a') * a * Z_ode(:, :, i-1));
    x_ode_hat(:, i) = x_ode_hat(:, i-1) + Z_ode(:, :, i) * a' * (I_real_ode(:, i) - a * x_ode_hat(:, i-1));
    Y_ode_hat(:, :, i) = reshape(Q * x_ode_hat(:, i), [n n]);
end