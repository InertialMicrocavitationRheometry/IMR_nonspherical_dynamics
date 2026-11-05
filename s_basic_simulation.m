%% Basic simulation of a nonspherical bubble
% Compute the spherical radius history with IMRv2, then use it to evolve
% small spherical-harmonic perturbations with and without rotational flow.
% Run from the repository root. The final two sections export strain fields
% to the current folder and can be skipped when only the histories are needed.
clearvars
clc
setup_paths;
tic

%% Spherical radius history
% Lengths are in meters, stresses and pressures in pascals, and time in
% seconds before nondimensionalization. Rmax is the initial bubble radius.
Rmax = 375e-6;      % Initial (maximum) radius [m]
Req = Rmax/2;    % Equilibrium radius [m]
mu = 0.05;         % Dynamic viscosity [Pa s]
G = 2.77e3;        % Shear modulus [Pa]
alph = 0;          % Nonlinear elastic stress parameter
sig = 0.056;       % Surface tension [N/m]
rho = 1048;        % Density for perturbation scaling [kg/m^3]
p8 = 101325;       % Reference pressure [Pa]

% Ultrasound settings are used only when ultra is true. The IMRv2 wrapper
% converts the nondimensional end time into the requested dimensional grid.
p_a = -50e3;       % Signed acoustic pressure amplitude [Pa]
f_a = 28e3;        % Acoustic frequency [Hz]
ultra = false;
tf_nd = 8;
tsteps = 1000;
[t, R, Rd, Rdd] = f_call_IMRv2(Rmax, Req, mu, G, alph, sig, p_a, f_a, tf_nd, tsteps, ultra);

figure('Name', 'Spherical radius history')
plot(t, R)
xlabel('Nondimensional time (IMRv2)')
ylabel('Radius / initial radius')

%% Perturbation initial conditions
% The solver builds a transformed radial grid from the bubble wall toward
% the far field. N can be a vector to evolve several spherical-harmonic modes.
xN = 256;          % Number of points on the transformed radial grid
N = 5;             % Spherical-harmonic degree, e.g. [2 5 8 11]
L = 1;             % Stretching parameter for the radial coordinate map
ep0 = -0.115 * ones(size(N));  % Initial nondimensional amplitude per mode
epd0 = zeros(size(N));        % Initial nondimensional amplitude rate

% Each mode starts with zero rotational field and zero field derivative.
% Rows index modes, and columns index points on the radial grid.
T0 = zeros(numel(N), xN);
Td0 = T0;

%% Convert to the perturbation solver's reference scales
% The perturbation model uses the equilibrium radius as its length scale.
% The inertial time scale is length times sqrt(density / pressure).
Lmax = Rmax/Req;
tcLIC = Rmax*sqrt(rho/p8);
Lc = Req;
rhoc = rho;
tc = sqrt(rhoc/p8)*Lc;
Uc = Lc/tc;
pc = rhoc*Uc^2;

% Cauchy, Reynolds, and Weber numbers enter the perturbation equations.
% Keep the remaining dimensionless groups available for reporting.
Ca = pc/G;
Re = Lc*sqrt(rhoc*pc)/mu;
We = pc*Lc/(2*sig);
Oh = sqrt(We)/Re;
De = Ca/Re;
Ec = sqrt(We)/Ca;

% Convert the supplied radius, acceleration, and time histories to the
% perturbation scaling. The dimensionless radial velocity is used as supplied.
R = R.*Rmax/Req;
Rdd = Rdd.*Req/Rmax;
t = t.*tcLIC/tc;

%% Optional prescribed perturbation
% Use 'F' to solve for the surface amplitude or 'T' to prescribe a sinusoid.
% The forced arrays have one row per mode and one column per time sample.
% The derivative is with respect to the same nondimensional time used above.
forcedep = 'F';
solverOptions = {'Verbose', true};
if forcedep == 'T'
    forcedEpAngularFrequency = 2*pi;
    forcedEpAmplitude = ep0(:).';
    if isscalar(forcedEpAmplitude)
        forcedEpAmplitude = repmat(forcedEpAmplitude, 1, numel(N));
    end

    if numel(forcedEpAmplitude) ~= numel(N)
        error('The prescribed amplitude must be scalar or have one value per mode in N.');
    end

    % Sample each mode's sinusoid and its analytical time derivative.
    forcedEp = zeros(numel(N), numel(t));
    forcedEpd = zeros(numel(N), numel(t));
    for iMode = 1:numel(N)
        forcedEp(iMode, :) = forcedEpAmplitude(iMode)*sin(forcedEpAngularFrequency.*t);
        forcedEpd(iMode, :) = forcedEpAmplitude(iMode)*forcedEpAngularFrequency.* ...
            cos(forcedEpAngularFrequency.*t);
    end

    solverOptions = [solverOptions, {'ForcedEp', forcedEp, 'ForcedEpd', forcedEpd}];
end

%% Evolve the perturbations
% Compare models on the same radius and time histories. The irrotational
% solve supplies the reference amplitudes. The rotational solve also
% supplies the radial field T used to reconstruct the strain.
equilibriumRadius = 1;  % Req in the perturbation solver's scaling
timeStepping = 2;       % Backward Euler startup followed by variable-step BDF2
rotationalModel = "me"; % The solver also accepts "Pros"
[epirr, ~, ~, ~, ~, ~, ~] = compute_rotational_perturbation_evolution( ...
    xN, L, N, ep0, epd0, T0, Td0, equilibriumRadius, R, Rd, Rdd, ...
    Ca, alph, Re, We, t, timeStepping, forcedep, rotationalModel, "irr", solverOptions{:});

[ep, epd, T, Td, R, Rd, t] = compute_rotational_perturbation_evolution( ...
    xN, L, N, ep0, epd0, T0, Td0, equilibriumRadius, R, Rd, Rdd, ...
    Ca, alph, Re, We, t, timeStepping, forcedep, rotationalModel, "rot", solverOptions{:});
toc

%% Compare surface amplitudes
% Use a separate figure for the perturbation amplitudes. Solid and dashed
% curves distinguish the models, and each mode uses a consistent color.
figure('Name', 'Surface perturbation amplitudes')
hold on
for iMode = 1:numel(N)
    modeCurve = plot(t./Lmax, ep(iMode, :), ...
        'DisplayName', sprintf('Mode %d, rotational', N(iMode)));
    plot(t./Lmax, epirr(iMode, :), '--', 'Color', modeCurve.Color, ...
        'DisplayName', sprintf('Mode %d, irrotational', N(iMode)));
end
hold off
xlabel('Nondimensional time (initial-radius scale)')
ylabel('Perturbation amplitude')
legend('Location', 'best')

%% Export a strain snapshot PDF
% Both exports show the radial-tangential component (e_rt) of the Eulerian
% Almansi strain in the current bubble geometry. Use the custom palette
% when it is available and MATLAB's parula palette in a bare checkout.
strainColormap = @parula;
if exist('brown_centered_diverging', 'file') == 2
    strainColormap = @brown_centered_diverging;
end

strainOptions = { ...
    'Rref', 1, ...
    'StrainMeasure', 'almansi', ...
    'StrainScalar', 'ert', ...
    'StrainColormap', strainColormap, ...
    'SymmetricCLim', true, ...
    'FEM_grid', true, ...
    'GridColor', [1 1 1]};

% Place twelve snapshots on one page and use a fixed strain range for
% comparison across panels. Overlay circles and rays on each snapshot.
make_axisym_displacement_movie_all_fields(T, ep, R, t, N, L, ...
    fullfile(pwd, 'strain_prelim.pdf'), strainOptions{:}, ...
    'OutputMode', 'snapshot_pdf', ...
    'SnapshotLayout', [3 4], ...
    'SnapshotWidthNormalized', 0.75, ...
    'SnapshotHeightNormalized', 0.85, ...
    'SnapshotTimeRange', [0 t(end)], ...
    'GridCircles', 24, ...
    'GridRays', 96, ...
    'tc', 1, ...
    'RLimEq', [1 6], ...
    'GridCircleSpacing', 'current', ...
    'StrainCLim', [-0.1 0.1], ...
    'SnapshotTileSpacing', 'compact', ...
    'SnapshotPadding', 'compact', ...
    'ColorbarLabel', 'Eulerian Almansi shear strain $e_{rt}$');

%% Export a strain movie
% Use the same strain definition with a wider color range and a coarser
% grid overlay so the deformation is easier to see during playback.
make_axisym_displacement_movie_all_fields(T, ep, R, t, N, L, ...
    fullfile(pwd, 'strain_movie.mp4'), strainOptions{:}, ...
    'OutputMode', 'movie', ...
    'MovieType', 'mp4', ...
    'WriteMovie', true, ...
    'FrameRate', 60, ...
    'MovieLength', 3, ...
    'ShowMovieColorbar', false, ...
    'StrainCLim', [-0.2 0.2], ...
    'cbScale', 'linear', ...
    'GridType', 'eulerian', ...
    'GridCircles', 32, ...
    'RadialSpacingPower', 1, ...
    'GridRays', 36, ...
    'RLimEq', [1 5], ...
    'MovieFigurePosition', [100 100 2160 2160]);
