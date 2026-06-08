%% basic simulation input
clear all
clc
% close all

addpath common


%
tic
% -------- Radial Solver ----------------------------------------%
Rmax = 150e-6;
Req = Rmax/3;
mu =  0.015;
G = 2.77e3;
alph = 0.25;
sig = 0.056;
p_a = -50e3; f_a = 28e3;
rho = 1048;
p8 = 101325;
tcLIC = Rmax*sqrt(rho/p8);
tf_nd = 4;%%max(texp1)/tcLIC; 
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
n = 8;
L = 5;

Lmax = Rmax/Req;
N = n;
forcedep = 'F';
ep0 = 0.025;
epd0 = 0;
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

mod = "me";
[ep, ~, ~, ~, ~, ~, ~] = compute_rotational_perturbation_evolution(xN, L, N, ep0, epd0, T0, ...
    Td0, 1, R, Rd, Rdd, Ca, alph, Re, We, t, 2, forcedep, mod, "irr");
epirr = ep;


[ep, epd, T, Td, R, Rd, t] = compute_rotational_perturbation_evolution(xN, L, N, ep0, epd0, T0, ...
    Td0, 1, R, Rd, Rdd, Ca, alph, Re, We, t, 2, forcedep, mod, "rot");
toc
%%

Lmax = Rmax/Req;
figure
plot(t./Lmax, R./Lmax)

figure
plot(t./Lmax, ep)
hold on
plot(t./Lmax, epirr(1, 1:length(t)), '--')
ylim([-.1 .1])
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


