function [Y11_hat, Y22_hat, Y_hat] = partitioned_solver(V, I, T, r, sym)
n = size(V, 1);
y11_size = n-r;
y22_size = r;

partitioned_v = T*V;
partitioned_i = T*I;

v1 = partitioned_v(1:end-r, :);
v2 = partitioned_v(end-r+1:end, :);
i1 = partitioned_i(1:end-r, :);
i2 = partitioned_i(end-r+1:end, :);

X = v1*pinv(v2);

C = i2*pinv(v2) - pinv(v2)'*i1'*X;
vec_c = reshape(C, [], 1);
b = vec_c;

if ~sym
    A11 = -kron(X', X');
    A22 = eye(y22_size^2, y22_size^2);
    A = [A11 A22];
    x = A\b;

    x1 = x(1:y11_size^2);
    x2 = x(y11_size^2+1:end);

    Y11_hat = reshape(x1, [y11_size, y11_size]);
    Y22_hat = reshape(x2, [y22_size, y22_size]);
else
    D1 = duplication_matrix(y11_size);
    D2 = duplication_matrix(y22_size);
    
    x1_size = (y11_size^2-y11_size)/2 + y11_size;
    x2_size = (y22_size^2-y22_size)/2 + y22_size;
    A = [-kron(X', X')*D1, D2];
    
    x = A\b;

    x1 = x(1:x1_size);
    x2 = x(x1_size+1:end);
    Y11_hat = reshape(D1*x1, [y11_size, y11_size]);
    Y22_hat = reshape(D2*x2, [y22_size, y22_size]);
end

Y12_hat = X'\(i2*pinv(v2)-Y22_hat);
Y_hat = inv(T)*[Y11_hat, Y12_hat; Y12_hat', Y22_hat]*T;
end

