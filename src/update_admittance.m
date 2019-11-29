function Y = update_admittance(Y_old, i, j, delta)
% update admittance matrix while respecting Laplacian constraints
Y = Y_old;
Y(i, j) = Y(i, j) + delta;
Y(j, i) = Y(j, i) + delta;
Y(i, i) = Y(i, i) - delta;
Y(j, j) = Y(j, j) - delta;
end

