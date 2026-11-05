function [A, rhs] = f_boundary_conditions(A, ep, epd, R, Rd, Req, ...
                                  N, xN, Ca, alph, Re, w, rhs, r, dr, ...
                                  forcedep, mod, dxdr, rot, H1)
%F_BOUNDARY_CONDITIONS Apply modal boundary conditions to A and rhs.
%   Replace selected rows of the implicit system with the interface and
%   far-field equations. Prescribed amplitudes enter the right-hand side.

sr = R/Req;

s = r;
ds = dr;

for i = 1:length(N)

    [idxT, idxV, idxep, idxepd, blockSize] = f_get_indicies(forcedep, mod, xN, i, rot);    

    % Evaluate the wall stress and the kappa integral weights for this mode.
    n = N(i);

    zeta = 1/Ca*Req*(2*alph*sr^6-3*alph*sr^4+alph+sr^4);
    fT = (zeta-4*sr^7*1/Re*Rd)/R^2;

    % The time derivative of the kappa integral has a wall-motion term
    % from differentiating its moving integration limit.
    fkappadopt = (n+1)/(2*n+1)*R^(-n)*Rd;

    j2kappa = 2*(2*n+1)/(n+1)*R^(n-3)*(zeta+(n-5)*sr^7/Re*Rd);
    j2kappadot = 2*(2*n+1)/(n+1)*R^(n-2)*sr^7/Re;

    [wkap] = f_kappa_kernel(s, ds, n, w);
    wkapdot = wkap;
    fK = j2kappadot*fkappadopt;

    if mod == "me"
        % Replace the first V row with the tangential stress condition at
        % the moving wall, including the mapped time-derivative term.
        A(idxV(1),(i-1)*blockSize+1:i*blockSize) = 0;

        A(idxV(1), idxV(1)) = A(idxV(1), idxV(1)) + sr^7/(R*Re);

        f_adBC = sr^7/(R*Re)*Rd*dxdr(1).*H1(1,:);
        A(idxV(1), idxT) = A(idxV(1), idxT) + f_adBC;

        % Coefficient for T fReqm BC
        A(idxV(1), idxT(1)) = A(idxV(1), idxT(1)) + fT;

        % Coefficient for T fReqm leibniz rule on kappadot
        A(idxV(1), idxT(1)) = A(idxV(1), idxT(1)) - fK;

        % Coefficient fReqm kappa needs full Reqw replacement on T
        A(idxV(1), idxT) = A(idxV(1), idxT) - j2kappa.*wkap;

        % Coefficient fReqm kappadot needs full Reqw replacement on V
        A(idxV(1), idxV) = A(idxV(1), idxV) - j2kappadot.*wkapdot;
        
        rhs(idxV(1)) = 0;

        % Impose zero T and V at the far-field grid point.

        A(idxT(end), (i-1)*blockSize+1:i*blockSize) = 0;

        A(idxT(end), idxT(end)) = 1;
        rhs(idxT(end)) = 0;

        A(idxV(end), (i-1)*blockSize+1:i*blockSize) = 0;
        A(idxV(end), idxV(end)) = 1;
        rhs(idxV(end)) = 0;


        % Keep the surface terms in the matrix when ep and epd are unknowns.
        if forcedep == 'T'
            rhs(idxV(1)) = rhs(idxV(1)) + 2*(n+2)/(n+1)/R*(ep(i)*(zeta-3*sr^7*1/Re*Rd)+sr^7*1/Re*R*epd(i));
        else
            A(idxV(1), idxep) = A(idxV(1), idxep) - 2*(n+2)/(n+1)/R*(zeta-3*sr^7*1/Re*Rd);
            A(idxV(1), idxepd) = A(idxV(1), idxepd) - 2*(n+2)/(n+1)/R*(sr^7*1/Re*R);
        end
    elseif mod == "Pros"
        % This formulation applies the wall condition in the first T row.
        A(idxT(1), (i-1)*blockSize+1:i*blockSize) = 0;
        A(idxT(1), idxT(1)) = A(idxT(1), idxT(1)) + (n+1)/R*1/Re;
        A(idxT(1), idxT) = A(idxT(1), idxT) - 2*(2*n+1)*R^(n-2)./Re*wkap;

        % Impose zero T at the outer boundary.
        A(idxT(end), :) = 0;

        A(idxT(end), idxT(end)) = 1;
        rhs(idxT(end)) = 0;

        if forcedep == 'T'
            rhs(idxT(1)) = 2*(n+2)./Re*(epd(i)*R+ep(i)*Rd)/R - 2*(n-1)./Re*ep(i)*Rd/R;
        else
            A(idxT(1), idxep) = A(idxT(1), idxep) -2*(n+2)./Re*(Rd)/R + 2*(n-1)./Re*Rd/R - (n+1)/Ca*sr^(-4);
            A(idxT(1), idxepd) =  A(idxT(1), idxepd) - 2*(n+2)./Re;
            rhs(idxT(1)) = 0;
        end
    end
end
end
