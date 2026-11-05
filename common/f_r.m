function r = f_r(x, R, L)
%F_R Map the transformed radial coordinate x to radius r.
%   x = -1 is the bubble wall at R. As x approaches 1, r tends to infinity.
%   L controls the radial stretching. r has the same shape as x.

r = R.*((1+x)./(1-x).*L+1);

end
