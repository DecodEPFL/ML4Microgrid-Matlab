function [] = progressive_partioning(Y, V, I, T, n_points)
    n = size(V, 1);
    ranks = round(linspace(1, n, n_points));
    
    ols_benchmark = conj(V' \ I');
    [ols_f_norm, ols_max_norm] = errors(Y, ols_benchmark);
    
    % First line: non-symmetric
    % Second line: symmetric
    max_norms = zeros(2, n_points);
    f_norms = zeros(2, n_points);
    y11_max_norms = zeros(2, n_points);
    y11_f_norms = zeros(2, n_points);
    y22_max_norms = zeros(2, n_points);
    y22_f_norms = zeros(2, n_points);
    
    
    
    [matlab_rank, T] = full_rank_partioning(V);
    
    for k = 1:n_points
        r = ranks(k);
        y11_size = n-r;
        y22_size = r;
        
        [Y11_hat, Y22_hat, Y_hat] = partitioned_solver(V, I, T, r, false);
        [Y11_hat_sym, Y22_hat_sym, Y_hat_sym] = partitioned_solver(V, I, T, r, true);
        
        partitioned_Y = T*Y*inv(T);
        y11 = partitioned_Y(1:y11_size, 1:y11_size);
        y22 = partitioned_Y(y11_size+1:end, y11_size+1:end);
        
        if r < n
            [y11_f_norms(1, k), y11_max_norms(1, k)] = errors(y11, Y11_hat);
            [y11_f_norms(2, k), y11_max_norms(2, k)] = errors(y11, Y11_hat_sym);
        end
        [y22_f_norms(1, k), y22_max_norms(1, k)] = errors(y22, Y22_hat);
        [y22_f_norms(2, k), y22_max_norms(2, k)] = errors(y22, Y22_hat_sym);
        [f_norms(1, k), max_norms(1, k)] = errors(Y, Y_hat);
        [f_norms(2, k), max_norms(2, k)] = errors(Y, Y_hat_sym);
        
    end
        figure();
        plot(ranks, y11_max_norms(1, :),'DisplayName','Non-symetric');
        hold on;
        plot(ranks, y11_max_norms(2, :),'DisplayName','Symetric');
        title('Y11 max norm');
        xlabel("Rank");
        ylabel("Max norm error");
        ylim([0 2*ols_max_norm])
        legend
        figure();
        plot(ranks, y11_f_norms(1, :),'DisplayName','Non-symetric');
        hold on;
        plot(ranks, y11_f_norms(2, :),'DisplayName','Symetric');
        title("Y11 frobenius norm")
        xlabel("Rank")
        ylabel("Frobenius norm error")
        ylim([0 2*ols_f_norm])
        legend
        
        figure();
        plot(ranks, y22_max_norms(1, :),'DisplayName','Non-symetric');
        hold on;
        plot(ranks, y22_max_norms(2, :),'DisplayName','Symetric');
        title("Y22 max norm")
        xlabel("Rank")
        ylabel("Max norm error")
        ylim([0 2*ols_max_norm])
        legend
        figure();
        plot(ranks, y22_f_norms(1, :),'DisplayName','Non-symetric');
        hold on;
        plot(ranks, y22_f_norms(2, :),'DisplayName','Symetric');
        title("Y22 frobenius norm")
        xlabel("Rank")
        ylabel("Frobenius norm error")
        ylim([0 2*ols_f_norm])
        legend
       
        figure();
        plot(ranks, max_norms(1, :), 'DisplayName', 'Non-symetric');
        hold on;
        plot(ranks, max_norms(2, :), 'DisplayName', 'Symetric');
        xline(matlab_rank, 'DisplayName', 'Matlab rank');
        yline(ols_max_norm, 'DisplayName', 'OLS benchmark');
        title("Y max norm")
        xlabel("Rank")
        ylabel("Max norm error")
        ylim([0 2*ols_max_norm])
        legend
        figure();
        plot(ranks, f_norms(1, :), 'DisplayName', 'Non-symetric');
        hold on;
        plot(ranks, f_norms(2, :), 'DisplayName', 'Symetric');
        xline(matlab_rank, 'DisplayName', 'Matlab rank');
        yline(ols_f_norm, 'DisplayName', 'OLS benchmark');
        title("Y frobenius norm")
        xlabel("Rank")
        ylabel("Frobenius norm error")
        ylim([0 2*ols_f_norm])
        legend
end

function [f_err, max_err] = errors(M, M_hat)
Me = M - M_hat;
max_err = max(max(abs(Me)));
f_err = norm(Me, 'fro');
end
