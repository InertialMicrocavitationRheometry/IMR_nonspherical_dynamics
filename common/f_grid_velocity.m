function a = f_grid_velocity(x, R, Rd, L)
%COMPUTE_GRID_VELOCITY Summary of this function goes here
%   Detailed explanation goes here

a = Rd/R.*(1-x).^2./(2*L).*(1+L.*(1+x)./(1-x));

end