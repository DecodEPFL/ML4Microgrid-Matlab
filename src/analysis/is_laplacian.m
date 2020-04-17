function [output] = is_laplacian(M)
%IS_LAPLACIAN Summary of this function goes here
%   Detailed explanation goes here

    output = issymmetric(M) & round(real(max(sum(M, 2))), 10) == 0;
end

