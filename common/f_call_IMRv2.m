function [t, R, Rd, Rdd] = f_call_IMRv2(Rmax, Req, mu, G, alph, sig, p_a, f_a, tf_nd, tsteps, ultra, varargin)
%F_CALL_IMRV2 Run the external IMRv2 radial solver and return bubble history.
%
%   f_call_IMRv2(..., ultra) looks for f_imr_fd on the MATLAB path, in the
%   IMRV2_FORWARD_SOLVER environment variable, or in ../IMRv2/src/forward_solver.
%
%   f_call_IMRv2(..., ultra, 'IMRv2Path', pathToForwardSolver) adds an explicit
%   IMRv2 forward-solver folder before calling f_imr_fd.
%   'Kappa' sets the gas polytropic exponent (default 1).
%   Pass radii, material properties, and acoustic inputs in SI units.
%   Outputs are time, radius, velocity, and acceleration in IMRv2's
%   nondimensional scaling, before conversion to the perturbation scaling.

opts = parse_imrv2_options(varargin{:});
imrv2Path = resolve_imrv2_forward_solver_path(opts.IMRv2Path);
if strlength(imrv2Path) > 0
    addpath(char(imrv2Path));
end

if exist('f_imr_fd', 'file') ~= 2
    error('f_call_IMRv2:MissingIMRv2', ['Could not find f_imr_fd. Run setup_paths ' ...
        'with the IMRv2 forward-solver folder, set IMRV2_FORWARD_SOLVER, or add ' ...
        '../IMRv2/src/forward_solver to the MATLAB path.']);
end

% Reference conditions internal to IMRv2. These values also set the radial
% problem's time scale, which is distinct from the perturbation time scale.
kappa = opts.Kappa;
T8 = 298.15;   % Ambient temperature [K]
rho8 = 1000;   % Reference density [kg/m^3]

% IMRv2 model selectors. The 0/1 switches control which physical
% processes enter the spherical solve.
radial = 2;    % Keller-Miksis equation in pressure form
vapor = 1;     % Include vapor
collapse = 0;  % Use the standard IMRv2 initial conditions
bubtherm = 1;  % Solve heat transfer inside the bubble
medtherm = 0;  % Disable the surrounding-medium temperature solve
masstrans = 1; % Include mass transfer
stress = 2;    % Quadratic neo-Hookean Kelvin-Voigt stress model

% IMRv2 expects angular frequency in radians per second. A zero pressure
% amplitude gives the unforced radial calculation when ultra is false.
if ultra
    pa = p_a;
    omega = 2*pi*f_a;
    wavetype = 4;
else
    pa = 0;
    omega = 0;
    wavetype = 0;
end

% IMRv2 accepts the requested time grid in seconds. Method 23 selects
% ode23tb, and Nt sets the thermal grid resolution inside IMRv2.
tc = Rmax*sqrt(rho8/101325);
tfin = tf_nd*tc;
tvector = linspace(0, tfin, tsteps);
radialOptions = { ...
    'progdisplay', 0, 'radial', radial, 'bubtherm', bubtherm, 'tvector', tvector, ...
    'vapor', vapor, 'medtherm', medtherm, 'masstrans', masstrans, 'method', 23, ...
    'stress', stress, 'collapse', collapse, 'mu', mu, 'g', G, 'lambda1', 0e-7, ...
    'lambda2', 0, 'alphax', alph, 'surft', sig, 'r0', Rmax, 'req', Req, 'kappa', kappa, 't8', T8, ...
    'rho8', rho8, 'pa', pa, 'omega', omega, 'wave_type', wavetype};

% Pressure, temperatures, and vapor concentration are also returned by
% IMRv2. The perturbation solver needs only the radial kinematics.
[t, R, Rd, ~, ~, ~, ~, Rdd] = f_imr_fd(radialOptions{:}, 'Nt', 75);
end

function imrv2Path = resolve_imrv2_forward_solver_path(explicitPath)
% Prefer an explicit folder, then the environment setting, then a sibling
% IMRv2 checkout. With none of these, rely on MATLAB's existing path.
explicitPath = string(explicitPath);
if strlength(explicitPath) > 0
    if exist(char(explicitPath), 'dir') ~= 7
        error('f_call_IMRv2:InvalidIMRv2Path', ...
            'IMRv2Path does not exist or is not a folder: %s', char(explicitPath));
    end
    imrv2Path = explicitPath;
    return;
end

envPath = string(getenv('IMRV2_FORWARD_SOLVER'));
if strlength(envPath) > 0
    if exist(char(envPath), 'dir') ~= 7
        error('f_call_IMRv2:InvalidIMRv2EnvironmentPath', ...
            'IMRV2_FORWARD_SOLVER does not exist or is not a folder: %s', char(envPath));
    end
    imrv2Path = envPath;
    return;
end

repoRoot = fileparts(fileparts(mfilename('fullpath')));
defaultPath = string(fullfile(fileparts(repoRoot), 'IMRv2', 'src', 'forward_solver'));
if exist(char(defaultPath), 'dir') == 7
    imrv2Path = defaultPath;
else
    imrv2Path = "";
end
end

function opts = parse_imrv2_options(varargin)
parser = inputParser;
parser.FunctionName = 'f_call_IMRv2';
addParameter(parser, 'IMRv2Path', "", @(p) ischar(p) || isstring(p));
addParameter(parser, 'Kappa', 1, @(v) isnumeric(v) && isscalar(v) && isfinite(v));
parse(parser, varargin{:});
opts = parser.Results;
end
