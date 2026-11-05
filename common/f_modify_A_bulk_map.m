function A = f_modify_A_bulk_map(A, N, xN, a, mod, blockSize)
%F_MODIFY_A_BULK_MAP Add advection from the moving radial coordinate.
%   The "me" model applies this correction to the V equation. The T
%   equation already receives its advection term in f_A_matrix. The "Pros"
%   model advances T directly, so the correction goes into that block.

a  = a(:);
h = 2/(xN-1);
Aadv = f_advective_D1_operator(a, h);

for k = 1:numel(N)
    if mod == "me"       
        base = (k-1)*blockSize;
        idxV  = base + xN + (1:xN);
        A(idxV, idxV) = A(idxV, idxV) + Aadv;

    elseif mod == "Pros"
        base = (k-1)*blockSize;
        idxT   = base + (1:xN);
        A(idxT, idxT) = A(idxT, idxT) + Aadv;
    end
end
end
