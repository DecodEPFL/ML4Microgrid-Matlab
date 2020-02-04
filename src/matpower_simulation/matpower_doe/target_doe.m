function z = target_doe(x, n, Z, trans_matrix)
%TARGET_DOE Target for design of experiment
    vm = x(n+1:2*n);
    a = kron(vm, eye(n)).' * trans_matrix;
    z = - log(abs(trace(inv(Z) + a' * a)));
end

