function [KmT, JmT] = f_J_K_kernel(s, ds, oneW, W, n)
%F_J_K_KERNEL Summary of this function goes here
%   Detailed explanation goes here

v1 = -s.^(-n).*ds;
KmT = (oneW-W)*diag(v1);

h1 = s.^(n+1).*ds;
JmT = W*diag(h1);

end