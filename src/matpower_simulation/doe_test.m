clear all
close all
clc
rng(11);

define_constants;
noise_sd = 1e-4;
mpc = loadcase('case6ww_modified');

mpopt_doe = mpoption('fmincon.alg', 6, 'opf.ac.solver', 'FMINCON', 'opf.start', 3);

Y = full(makeYbus(mpc));
n = size(Y, 1);
load_mask = mpc.bus(:, BUS_TYPE) == PQ;
controllable_mask = mpc.bus(:, BUS_TYPE) == PV;

load_avg_pd = mpc.bus(load_mask, PD);
load_avg_qd = mpc.bus(load_mask, QD);

Q = duplication_matrix(n) * transformation_matrix(n);
doe_free_params = n * (n - 1) / 2; 
% x_hat(:, 1) = [2-4i, 0, 1.1-4.7i, 0.8-3.1i, 0, 0.7-3.8i, 4-8i, 1-3i, 1.5-4.4i, 0, 1.4-3.1i, 1.9-9.6i, 1-2i, 0, 1-3i];
x_hat(:, 1) = 0.1 + rand([doe_free_params, 1]);
Z(:, :, 1) = 100 * eye(doe_free_params);
Y_hat(:, :, 1) = reshape(Q * x_hat(:, 1), [n n]);

line_sample = mpc.branch(1, :);
mpc_orig = mpc;
mpc_iteration = mpc_orig;

lambda = 0.7;
samples = 20;

for i = 2:samples
    mpc_iteration.bus(load_mask, PD) = normrnd(load_avg_pd, 0.05 * load_avg_pd);
    mpc_iteration.bus(load_mask, QD) = normrnd(load_avg_qd, 0.05 * load_avg_qd);
    
    mpc = mpc_iteration;
    mpc.branch = [];
    for j=1:n
        for k=j+1:n
            this_line = line_sample;
            this_line(1, F_BUS) = j;
            this_line(1, T_BUS) = k;
            this_line(1, BR_R) = - real(1 / Y_hat(j, k, i-1));
            this_line(1, BR_X) = - imag(1 / Y_hat(j, k, i-1));
            mpc.branch = [mpc.branch; this_line];
        end
    end
    Y_reconstruct = makeYbus(mpc);
    display(['Difference between estimated and reconstruced matrix ' num2str(max(max(abs(Y_hat(:, :, i-1) - Y_reconstruct))))]);
    res_doe = runopf_doe(mpc, mpopt_doe, Z(:, :, i-1), Q);
    V_star_ode(:, i) = res_doe.bus(:, VM);
    
    mpc = mpc_iteration;
    for j=1:length(controllable_mask)
        if controllable_mask(j)
            mpc.gen(mpc.gen(:, 1) == j, VG) = res_doe.bus(j, VM);
        end
    end
    res_pf = runpf(mpc);
    
    voltage_phase = res_pf.bus(:, VA) * pi / 180;
    V_real_ode(:, i) = res_pf.bus(:, VM).*(cos(voltage_phase) + sqrt(-1)*sin(voltage_phase));
    I_real_ode(:, i) = Y * V_real_ode(:, i) + normrnd(0, noise_sd, size(V_real_ode(:, i)));
    
    [Y_hat(:, :, i), x_hat(:, i), Z(:, :, i)] = rls_step(V_real_ode(:, i), I_real_ode(:, i), Q, lambda, x_hat(:, i-1), Z(:, :, i-1));    
end

error_metrics(Y, Y_hat(:, :, end), 'DoE');
figure()
error_evolution(Y, Y_hat, 'DoE');
