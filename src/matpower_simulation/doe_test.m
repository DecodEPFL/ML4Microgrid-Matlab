clear all
close all
clc
addpath('./src/');
rng(11);

define_constants;
mpc = loadcase('case4gs');

mpopt_doe = mpoption;
mpopt_doe = mpoption(mpopt_doe, 'fmincon.alg', 2, 'opf.ac.solver', 'FMINCON');

Y = full(makeYbus(mpc));
n = size(Y, 1);
load_mask = mpc.bus(:, BUS_TYPE) == PQ;
controllable_mask = mpc.bus(:, BUS_TYPE) == PV;

Q = duplication_matrix(n) * transformation_matrix(n);
doe_free_params = n * (n - 1) /2; 
x_doe_hat(:, 1) = zeros([doe_free_params, 1]);
Z_doe(:, :, 1) = 100 * eye(doe_free_params);
Y_doe_hat(:, :, 1) = reshape(Q * x_doe_hat(:, 1), [n n]);

line_sample = mpc.branch(1, :);
mpc_orig = mpc;
mpc_iteration = mpc_orig;

lambda = 0.8;
samples = 100;

load_avg_pd = mpc.bus(load_mask, PD);
load_avg_qd = mpc.bus(load_mask, QD);

for i = 2:samples
    %mpc_iteration.bus(load_mask, PD) = normrnd(load_avg_pd, 0.05 * load_avg_pd);
    %mpc_iteration.bus(load_mask, QD) = normrnd(load_avg_qd, 0.05 * load_avg_qd);
    
    mpc = mpc_iteration;
%     mpc.branch = [];
%     for j=1:n
%         for k=j+1:n
%             this_line = line_sample;
%             this_line(1, F_BUS) = j;
%             this_line(1, T_BUS) = k;
%             this_line(1, BR_R) = abs(real(1 / Y_doe_hat(j, k, i-1)));
%             this_line(1, BR_X) = abs(imag(1 / Y_doe_hat(j, k, i-1)));
%             mpc.branch = [mpc.branch; this_line];
%         end
%     end
    res_doe = runopf_doe(mpc, mpopt_doe, Z_doe(:, :, i-1), Q);
    V_star_ode(:, i) = res_doe.bus(:, VM);
    
    mpc = mpc_iteration;
    mpc.gen(1:end-1, VG) = res_doe.bus(controllable_mask, VM);
    res_pf = runpf(mpc);
    
    voltage_phase = res_pf.bus(:, VA) * pi / 180;
    V_real_ode(:, i) = res_pf.bus(:, VM).*(cos(voltage_phase) + sqrt(-1)*sin(voltage_phase));
    I_real_ode(:, i) = Y * V_real_ode(:, i);
    
    a = kron(V_real_ode(:, i), eye(n)).' * Q;
    Z_doe(:, :, i) = (1/lambda) * (Z_doe(:, :, i-1) - Z_doe(:, :, i-1) * a' * inv(lambda * eye(n) + a * Z_doe(:, :, i-1) * a') * a * Z_doe(:, :, i-1));
    x_doe_hat(:, i) = x_doe_hat(:, i-1) + Z_doe(:, :, i) * a' * (I_real_ode(:, i) - a * x_doe_hat(:, i-1));
    Y_doe_hat(:, :, i) = reshape(Q * x_doe_hat(:, i), [n n]);    
end

error_metrics(Y, Y_doe_hat(:, :, end), 'DoE');
figure()
error_evolution(Y, Y_doe_hat, 'DoE');