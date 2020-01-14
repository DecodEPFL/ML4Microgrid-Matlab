function D = duplication_matrix(n)
%duplication matrix Create duplication matrix of size n
%   More details at https://www.mathworks.com/matlabcentral/answers/473737-efficient-algorithm-for-a-duplication-matrix
m   = n * (n + 1) / 2;
nsq = n^2;
r   = 1;
a   = 1;
v   = zeros(1, nsq);
for i = 1:n
   v(r:r + i - 2) = i - n + cumsum(n - (0:i-2));
   r = r + i - 1;
   
   v(r:r + n - i) = a:a + n - i;
   r = r + n - i + 1;
   a = a + n - i + 1;
end
D = sparse(1:nsq, v, 1, nsq, m);
end

