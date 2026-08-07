%% basic simulation input
clear all
clc
% close all

addpath common


%
tic
% -------- Radial Solver ----------------------------------------%
Rmax = 500e-6;
Req = Rmax/7;
mu =  0.1;
G = 0e3;
alph = 0;
sig = 0.056;
p_a = -50e3; f_a = 28e3;
rho = 1048;
p8 = 101325;
tcLIC = Rmax*sqrt(rho/p8);
tf_nd = 2;%%max(texp1)/tcLIC; 
tsteps = 5000; ultra = false;

t = linspace(0, tf_nd, tsteps);
[t, R, Rd, Rdd] = f_call_IMRv2(Rmax, Req, mu, G, alph, sig, p_a, f_a, tf_nd, tsteps, ultra);

figure
plot(t, R, '-')
%%

% -------- perturbation solver initial conditions ---------------%
% define grid in transformed domain
xN = 256;
% Mode numbers
n = 2;
L = 2;

Lmax = Rmax/Req;
N = n;
forcedep = 'F';
ep0 = 0.001.*ones(size(N));
epd0 = zeros(size(N));
T0 = 0.*ones(length(N), xN);
Td0 = T0;
T0(end) = 0;

% Characteristic scales
Lc = Req;
rhoc = rho;
tc = sqrt(rhoc/p8)*Lc;
Uc = Lc/tc;
pc = rhoc*Uc^2;

Ca = pc/G;
Re = Lc*sqrt(rhoc*pc)/mu;
We = pc*Lc/(2*sig);
Oh = sqrt(We)/Re;
De = Ca/Re;
Ec = sqrt(We)/Ca;

% -------- Re-nondimensionalize to perturbation scalings ------- %
R = R.*Rmax/Req;
Rdd = Rdd.*Req/Rmax;
t = t.*tcLIC/tc;

% -------- prescribed perturbation history, if forcedep == 'T' -------- %
solverOptions = {};
if forcedep == 'T'
    forcedEpAngularFrequency = 2*pi;
    forcedEpPeriod = 2*pi/forcedEpAngularFrequency;

    forcedEpAmplitude = ep0(:).';
    if isscalar(forcedEpAmplitude) && numel(N) > 1
        forcedEpAmplitude = repmat(forcedEpAmplitude, 1, numel(N));
    end

    if numel(forcedEpAmplitude) ~= numel(N)
        error('forcedEpAmplitude must be scalar or have one value per mode in N.');
    end

    forcedEp = zeros(numel(N), numel(t));
    forcedEpd = zeros(numel(N), numel(t));
    for iMode = 1:numel(N)
        forcedEp(iMode, :) = forcedEpAmplitude(iMode)*sin(forcedEpAngularFrequency.*t);
        forcedEpd(iMode, :) = forcedEpAmplitude(iMode)*forcedEpAngularFrequency.* ...
            cos(forcedEpAngularFrequency.*t);
    end

    solverOptions = {'ForcedEp', forcedEp, 'ForcedEpd', forcedEpd, ...
        'ForcedEpPeriod', forcedEpPeriod};
end
solverOptions = {'Verbose', true};

mod = "me";
[ep, ~, ~, ~, ~, ~, ~] = compute_rotational_perturbation_evolution(xN, L, N, ep0, epd0, T0, ...
    Td0, 1, R, Rd, Rdd, Ca, alph, Re, We, t, 2, forcedep, mod, "irr", solverOptions{:});
epirr = ep;


[ep, epd, T, Td, R, Rd, t] = compute_rotational_perturbation_evolution(xN, L, N, ep0, epd0, T0, ...
    Td0, 1, R, Rd, Rdd, Ca, alph, Re, We, t, 2, forcedep, mod, "rot", solverOptions{:});

[eppros, epdpros, ~, ~, ~, ~, ~] = compute_rotational_perturbation_evolution(xN, L, N, ep0, epd0, T0, ...
    Td0, 1, R, Rd, Rdd, Ca, alph, Re, We, t, 2, forcedep, "Pros", "rot", solverOptions{:});
toc
%%

Lmax = Rmax/Req;
% figure
hold on
plot(t./Lmax, R./Lmax)

for i = 1:length(n)
figure
plot(t./Lmax, ep(i,:))
hold on
plot(t./Lmax, epirr(i, 1:length(t)), '--')
plot(t./Lmax, eppros(i,:), '-.')
ylim([-.1 .1])
end
%%
make_axisym_displacement_movie_all_fields(T, ep, R, t, N, L, fullfile(pwd,'strain_test.pdf'), ...
    'Req', Req, ...
    'StrainMeasure', 'almansi', ... 
    'StrainScalar', 'ert', ...
    'OutputMode', 'snapshot_pdf', ...
    'SnapshotLayout', [3 4], ...
    'SnapshotWidthNormalized', 0.75, ...
    'SnapshotHeightNormalized', 0.85, ...
    'SnapshotTimeRange', [0 t(end)], ...
    'GridCircles', 24, ...
    'GridRays', 96, ...
    'tc', 1, ...
    'RLimEq', [1 6], ...
    'FEM_grid', true, ...
    'StrainColormap', parula(256), 'StrainCLim', [-0.1 0.1], ...
    'SymmetricCLim', true,  'SnapshotTileSpacing', 'compact', ...
    'SnapshotPadding', 'compact', 'ColorbarLabel', 'Eulerian Almansi shear strain $e_{rt}$');


