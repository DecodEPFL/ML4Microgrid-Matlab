function Q = transformation_matrix(n)
%transformation matrix Create transformation matrix of size n
%   A tranformation matrix is a linear transformation to convert the
%   non-redundant vectorization of a Laplacian matrix into the
%   vectorization of a symmetric matrix
Q = zeros([n*(n+1)/2, n*(n-1)/2]);
row = 1;
for i = 1:n
    for j = 0:n-i
        if j == 0
            Q(row, (row - i +1):(row - i + n - i)) = 1;
            for k = 0:i-2
                s = 0;
                for x = 0:k-1
                    s = s + n - 2 - x;
                end
                Q(row, i - 1 + s) = 1;
            end
        else
            Q(row, row-i) = -1;
        end
        row = row + 1;
    end
end
end

