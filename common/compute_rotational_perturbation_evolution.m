function [ep_tot, epd_tot, T_tot, Td_tot, R, Rd, t] = compute_rotational_perturbation_evolution(xN, L, N, ep0, epd0, T0, ...
                    Td0, Req, R, Rd, Rdd, Ca, alph, Re, We, t, TSM, forcedep, mod, rot, varargin)
%COMPUTE_ROTATIONAL_PERTURBATION_EVOLUTION Compute surface perturbation evolution.
%   Use 'Verbose', true to print per-percent condition-number diagnostics.

addpath(fileparts(mfilename('fullpath')));

parser = inputParser;
parser.FunctionName = 'compute_rotational_perturbation_evolution';
addParameter(parser, 'Verbose', false, @(v) islogical(v) && isscalar(v));
parse(parser, varargin{:});
verbose = parser.Results.Verbose;

% xN: number of grid points
% N: a vector containing mode numbers to be simulated
% ep0: vector of epsilon values for each mode at the initial condition
% epd0: vector of time derivative for each mode at initial condition
% T0, T: length(N) x xN, initial condition/ current T field

% T_tot: 3D array,  length(N) x xN x length(t), contains T field for all
% modes for all timesteps

% Td_tot: 3D array,  length(N) x xN x length(t), contains Td field for all
% modes for all timesteps

% ep_tot: length(N) x length(t) all perturbation evolution data
% epd_tot: length(N) x length(t) all time derivative perturbation evolution data
nmode = length(N);
tsteps = length(t);

% Create grid
eps = 1e-4;
x_g = linspace(-1, 1-eps, xN);

% preallocate storage for outputted variables
ep_tot = zeros(nmode, tsteps);
epd_tot = zeros(nmode, tsteps);
T_tot = zeros(nmode, xN, tsteps);
Td_tot = zeros(nmode, xN, tsteps);

% create solution at timestep p
ep = ep0; epd = epd0; T = T0; Td = Td0;

% Assign initial conditions to output arrays
for i = 1:nmode
    ep_tot(i, 1) = ep(i);
    epd_tot(i, 1) = epd(i);
    T_tot(i, :, 1) = T(i,:);
    Td_tot(i, :, 1) = Td(i,:);
end

if forcedep == 'T'
    % Manufacture solution for ep and ep dot for now
    f = 2*pi;
    for i = 1:nmode
        ep_tot(i,:) = ep0(i)*sin(f.*t);
        epd_tot(i,:) = f.*ep0(i)*cos(f.*t);
    end
    P = 2;                 % period in NONDIM time (since t has been scaled by tc)
end

if (isempty(R) && isempty(Rd))
    % Manufacture solution for R and Rd for now
    a = 0.05;
    f = 2*pi/4;

    R = 1 + a*cos(f.*t);
    Rd = (-2*f.*a*sin(f.*t));
    Rdd = (-4*f.^2.*a*cos(f.*t));
    Req = 1;
    % 
    R = 1 + 0*R;
    Rd = 0*Rd;
    Rdd = 0*Rdd;
end

% create state vector 

[~, ~, ~, ~, blockSize] = f_get_indicies(forcedep, mod, xN, 1, rot);
q = zeros(nmode*blockSize,1);
for i = 1:nmode
        [idxT, idxV, idxep, idxepd, blockSize] = f_get_indicies(forcedep, mod, xN, i, rot);
    if rot == "rot"
        q(idxT) = T(i,:);
        if ~isempty(idxV)
            q(idxV) = Td(i,:);
        end
    end
        if ~isempty(idxep)
            q(idxep) = ep(i,:);
        end
        if ~isempty(idxepd)
            q(idxepd) = epd(i,:);
        end
end

q_nm1 = [];

% --- progress tracking ---
Nt = tsteps - 1;      % number of timesteps
nextPct = 1;             % next percentage to print

% ---- quasi-equilibrium detection settings ----
eqTol = 1e-3;               % relative tolerance for period-averaged T change
nConsec = 4;                % require this many consecutive "small change" detections
eqHits = 0;                 % counter

% Precompute difference stencils and integration matrices (xN x xN)
h = x_g(2)-x_g(1);
H1 = f_difference_stencil(xN, 1, h);
H2 = f_difference_stencil(xN, 2, h);
H3 = f_difference_stencil(xN, 3, h);
[W, w, One_wT] = f_trapz_integration_matrices(x_g);


% ---------- cache setup for constant-matrix runs ----------
tolCache = 1e-5;

if numel(t) >= 2
    dt0 = t(2) - t(1);
    uniform_dt = all(abs(diff(t) - dt0) < tolCache*max(1,abs(dt0)));
else
    uniform_dt = true;
end

% Only cache if the radius history is truly constant and dt is uniform
useCache = ...
    all(abs(R  - R(1))  < tolCache) && ...
    all(abs(Rd)         < tolCache) && ...
    all(abs(Rdd)        < tolCache) && ...
    uniform_dt;

% Separate caches for BE and BDF2
dAlhs_BE   = [];
dAlhs_BDF2 = [];

AlhsBase_BE   = [];
AlhsBase_BDF2 = [];

rowScale_BE   = [];
rowScale_BDF2 = [];

Alhs = [];   % keep for diagnostics/condest

% Timestep until end of simulation
for i = 1:tsteps-1
    % compute necessary values for domain mapping:
    dxdr = f_diff_var_change(x_g, R(i+1), L);
    [dxdr2, dxdr21] = f_diff2_var_change(x_g, R(i+1), L);
    [dxdr3, dxdr32, dxdr31] = f_diff3_var_change(x_g, R(i+1), L);
    drdx = f_ds(x_g, R(i+1), L);
    r = f_r(x_g, R(i+1), L);
    a = f_grid_velocity(x_g, R(i+1), Rd(i+1), L);
    
    % commonly used values
    arg = r.^3 - R(i+1).^3 + Req.^3;
    rs  = sign(arg).*abs(arg).^(1/3);

    sr = R(i+1)/Req;
    vs = Req^3*(sr^3-1);

    % construct source term vector 
    s = zeros(nmode*blockSize,1);

    if forcedep == 'T' && mod == "me" && rot == "rot"
        s = f_source_terms(ep_tot(:, i+1), epd_tot(:, i+1), R(i+1), Rd(i+1), ...
                         Rdd(i+1), r, N, xN, Ca, alph, Re, At1, rs, vs, blockSize);
    end

    % -----------------------------------------------------------------------------
    % Need to construct the RHS at t^n+1 though \dot{q} = A*q
    % -----------------------------------------------------------------------------
    dt_n = t(i+1) - t(i);

    if (i == 1 || TSM == 1)
        % -------------------------
        % Step 1: Backward Euler
        % -------------------------
        beta0 = 1.0;
        rhs = q + dt_n*s;

    elseif (i > 1 && TSM == 2)
        % -------------------------
        % Steps >=2: variable-step BDF2
        % -------------------------
        dt_nm1 = t(i) - t(i-1);
        dt_rat = dt_n / dt_nm1;

        beta0 = (1 + 2*dt_rat)/(1 + dt_rat);
        beta1 = 1 + dt_rat;
        beta2 = -dt_rat^2/(1 + dt_rat);

        rhs = beta1*q + beta2*q_nm1 + dt_n*s;
    end

   isBE = (i == 1 || TSM == 1);

   if ~useCache
       buildThisStep = true;
   else
       if isBE
           buildThisStep = isempty(dAlhs_BE);
       else
           buildThisStep = isempty(dAlhs_BDF2);
       end
   end

   if buildThisStep
       % ----- Build A and Alhs for this scheme -----
       A = f_A_matrix(r, R(i+1), Rd(i+1), Rdd(i+1), Req, N, xN, ...
           H1, H2, H3, w, Re, Ca, We, alph, a, dxdr, dxdr2, dxdr21, ...
           dxdr3, dxdr32, dxdr31, drdx, One_wT, W, rs, vs, sr, forcedep, blockSize, mod, rot);

       if forcedep ~= 'T' && rot ~= "irr"
           A = f_ep_source_terms(r, xN, N, A, ...
               R(i+1), Rd(i+1), Rdd(i+1), Req, H1, dxdr, drdx, w, W, One_wT, ...
               Re, Ca, alph, vs, rs, sr, mod, forcedep, rot);
       end

       if rot == "rot"
           A = f_modify_A_bulk_map(A, N, xN, a, mod, blockSize);
       end

       I = eye(size(A));
       AlhsBase = beta0*I - dt_n*A;   % <-- unmodified base matrix
        
       % Apply BCs to a WORKING copy, not to the cached base matrix

       if rot == "rot"
           [AlhsWork, rhs] = f_boundary_conditions(AlhsBase, ep_tot(:, i+1), epd_tot(:, i+1), ...
               R(i+1), Rd(i+1), Req, N, xN, Ca, alph, Re, w, rhs, r, drdx, ...
               forcedep, mod, dxdr, rot, H1);
       else
           AlhsWork = AlhsBase;
       end

       AlhsWork = sparse(AlhsWork);

       % ------------------------------------------------------------
       % Row-scale epsilon rows of AlhsWork and rhs
       % ------------------------------------------------------------

       rowsToScale = zeros(2*nmode,1);
       nRowsToScale = 0;

       for iMode = 1:nmode
           [~, ~, idxep, idxepd, ~] = f_get_indicies(forcedep, mod, xN, iMode, rot);

           if ~isempty(idxep)
               nRowsToScale = nRowsToScale + 2;
               rowsToScale(nRowsToScale-1:nRowsToScale) = [idxep; idxepd];
           end
       end

       rowsToScale = rowsToScale(1:nRowsToScale);
       rowsToScale = unique(rowsToScale(:));  % safety in case of duplicates

       rowScale = ones(size(rhs));

       if ~isempty(rowsToScale)

           % Infinity norm of selected rows:
           rowNorms = full(max(abs(AlhsWork(rowsToScale,:)), [], 2));

           rowScale(rowsToScale) = max(rowNorms, 1);

           % Sparse diagonal left-scaling matrix.
           nEq = size(AlhsWork,1);

           d = ones(nEq,1);
           d(rowsToScale) = 1 ./ rowScale(rowsToScale);

           Dscale = spdiags(d, 0, nEq, nEq);

           AlhsWork = Dscale * AlhsWork;
           rhs      = d .* rhs;
       end

       dThis = AlhsWork;
       dAlhs = dThis;
       Alhs  = AlhsWork;   % keep latest matrix for diagnostics

       % Save cache if allowed
       if useCache
           if isBE
               dAlhs_BE   = dThis;
               AlhsBase_BE = AlhsBase;
               rowScale_BE = rowScale;
           else
               dAlhs_BDF2   = dThis;
               AlhsBase_BDF2 = AlhsBase;
               rowScale_BDF2 = rowScale;
           end
       end

   else
       % ----- Reuse factorization: ONLY update rhs -----
       if isBE
           dAlhs = dAlhs_BE;
           if rot == "rot"
               [~, rhs] = f_boundary_conditions(AlhsBase_BE, ep_tot(:, i+1), epd_tot(:, i+1), ...
                   R(i+1), Rd(i+1), Req, N, xN, Ca, alph, Re, w, rhs, r, drdx, ...
                   forcedep, mod, dxdr, rot, H1);
           end

           rowsToScale = (rowScale_BE ~= 1);
           rhs(rowsToScale) = rhs(rowsToScale) ./ rowScale_BE(rowsToScale);

       else
           dAlhs = dAlhs_BDF2;
           if rot == "rot"
               [~, rhs] = f_boundary_conditions(AlhsBase_BDF2, ep_tot(:, i+1), epd_tot(:, i+1), ...
                   R(i+1), Rd(i+1), Req, N, xN, Ca, alph, Re, w, rhs, r, drdx, ...
                   forcedep, mod, dxdr, rot, H1);
           end
           rowsToScale = (rowScale_BDF2 ~= 1);
           rhs(rowsToScale) = rhs(rowsToScale) ./ rowScale_BDF2(rowsToScale);
       end
   end

   q_pp1 = dAlhs \ rhs;

    % Save fields at i+1
    for j = 1:nmode
        [idxT, idxV, idxep, idxepd, blockSize] = f_get_indicies(forcedep, mod, xN, j, rot);
        if rot == "rot"
            T_tot(j,:,i+1)  = q_pp1(idxT);
            if ~isempty(idxV)
                Td_tot(j,:,i+1) = q_pp1(idxV);
            end
        end
        if ~isempty(idxep)
            ep_tot(j,i+1) = q_pp1(idxep);
        end
        if ~isempty(idxepd)
            epd_tot(j,i+1) = q_pp1(idxepd);
        end
    end

    % shift history:
    q_nm1 = q;
    q = q_pp1;

    if verbose
        pctDone = floor(100 * i / Nt);
        if pctDone >= nextPct
            cdeest = condest(Alhs);
            rdo = rcond(full(Alhs));

            fprintf('%3d%% complete | condest(Alhs) = %.3e | rcond(Alhs) = %.3e\n', ...
                pctDone, cdeest, rdo);

            nextPct = pctDone + 1;
        end
    end

    if forcedep == 'T'
        [isEquil, eqHits, ~, pinc] = check_quasi_equil(r, ...
            t(i+1), t, P, T_tot, Td_tot, nmode, xN, ...
            eqTol, nConsec, eqHits, i+1);

        if isEquil
            T_tot(:,:, i+1:end)  = [];
            Td_tot(:,:, i+1:end) = [];
            ep_tot(:,i+1:end)    = [];
            epd_tot(:,i+1:end)   = [];

            % plot the steady state for all previousa steady values:
            x_g = linspace(-1, 1, xN);

            i_plot  = round(linspace(i-8*pinc, i, 17));
            Nt_plot = numel(i_plot);

            C = magma(Nt_plot+4);
            nPlots = numel(N);
            target = 1.0;

            % --- choose rows/cols for tiledlayout ---
            bestScore = inf; bestRC = [1 nPlots];
            for nRows = 1:ceil(sqrt(nPlots))
                nCols  = ceil(nPlots/nRows);
                tileAR = nCols/nRows;
                score  = abs(tileAR - target) + 0.05*(nRows*nCols - nPlots);
                if score < bestScore
                    bestScore = score;
                    bestRC = [nRows, nCols];
                end
            end
            nRows = bestRC(1); nCols = bestRC(2);

            % USER CONTROLS (edit these)
            axTickFS      = 18;
            axLabelFS     = 28;
            cbTickFS      = 18;
            cbLabelFS     = 50;

            axTickInterp  = 'latex'; axLabelInterp = 'latex'; 
            cbTickInterp  = 'latex'; cbLabelInterp = 'latex';

            % Figure + layout
            figure('Color','w');
            tl = tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');
            title(tl, '$T_n^*(r^*,t^*)$ evolution', 'Interpreter','latex','FontSize',24);

            axesH = gobjects(nPlots,1);

            % Plot each mode in its own tile
            for j = 1:nPlots
                axesH(j) = nexttile(tl);
                hold(axesH(j),'on');

                for ii = 1:Nt_plot
                    s = i_plot(ii);
                    r = f_r(x_g(:), 1, L);
                    plot(axesH(j), r, T_tot(j,:,s), '-', ...
                        'LineWidth', 1.5, 'Color', C(ii+2,:));
                end

                grid(axesH(j),'on'); box(axesH(j),'on');

                % --- axes labels (text) ---
                xlabel(axesH(j), '$r^*$', 'Interpreter', axLabelInterp);
                ylabel(axesH(j), '$T_n^*(r^*,t^*)$', 'Interpreter', axLabelInterp);

                % --- axes ticks (tick labels) ---
                axesH(j).TickLabelInterpreter = axTickInterp;
                axesH(j).FontSize = axTickFS;

                % --- axes label font sizes (independent) ---
                axesH(j).XLabel.FontSize = axLabelFS;
                axesH(j).YLabel.FontSize = axLabelFS;

                xlim(axesH(j), [0.9 25]);
                ylim(axesH(j), [max(-1, min(T_tot(j,1,:))), min(1, max(T_tot(j,1,:)))])
            end
            % Colormap + one shared colorbar (OLDER MATLAB compatible)
            colormap(C)

            % Create colorbar attached to an axes (works on older versions)
            cb = colorbar(axesH(end));
            cb.Layout.Tile = 'east';     % put it on the right of the tiledlayout

            % --- colorbar label (text) ---
            cb.Label.String      = '$t^*$';
            cb.Label.Interpreter = cbLabelInterp;
            cb.Label.FontSize    = cbLabelFS;

            % --- colorbar ticks (tick labels) ---
            cb.Ticks = linspace(0,1,5);
            cb.TickLabels = compose('$%.3g$', ...
                linspace(t(i_plot(1)), t(i_plot(end)), 5));
            cb.TickLabelInterpreter = cbTickInterp;
            cb.FontSize = cbTickFS;
            return
        end
    end

    Nt = length(t);

    % uncomment and change conditions to stop running simulations
    % if any(isnan(q)) || max(abs(ep_tot), [], "all") > 5 || max(abs(epd_tot), [], "all") > 25 || ... 
    %         (Nt > 1500 && i > floor(Nt/3)+1 && max(abs(ep_tot(:, i-floor(Nt/3):i)), [], "all") < 1.1e-3) %|| max(abs(T_tot), [], "all") > 10
    %     T_tot(:,:, i+1:end) = [];
    %     Td_tot(:,:, i+1:end) = [];
    %     ep_tot(:,i+1:end) = [];
    %     epd_tot(:,i+1:end) = [];
    %     t(i+1:end) = [];
    %     R(i+1:end) = [];
    %     Rd(i+1:end) = [];
    %     fprintf('Quit at timestep %d\n', i)
    %     return
    % end
end
end







