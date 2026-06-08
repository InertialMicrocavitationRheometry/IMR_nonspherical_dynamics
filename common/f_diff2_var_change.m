function [dxdr2_coef, dxdr_coef] = f_diff2_var_change(x, R, L)
%F_DIFF2_VAR_CHANGE Summary of this function goes here
%   Detailed explanation goes here

dxdr_coef = -(1-x).^3./(2*R^2*L^2);
dxdr2_coef = (1-x).^4./(4*R^2*L^2);


end