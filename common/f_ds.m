function ds = f_ds(x, R, L)
%COMPUTE_DIFF_INT_ELEMENT Summary of this function goes here
%   Detailed explanation goes here

ds = 2*R.*L./(1-x).^2;

end