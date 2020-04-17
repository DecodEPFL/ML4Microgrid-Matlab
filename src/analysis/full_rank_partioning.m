function [r, T] = full_rank_partioning(M)
    %FULL_RANK_PARTIONING Cumpute the rank of martix M and T matrix such as T*M
    %= [n-r dependant rows; r independant rows]
    %   Detailed explanation goes here
    r = rank(M);
    [~, ~, E] = qr(M');
    T = flip(E');
end

