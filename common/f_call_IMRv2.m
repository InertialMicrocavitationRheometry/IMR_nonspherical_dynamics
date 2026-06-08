function [t, R, Rd, Rdd] = f_call_IMRv2(Rmax, Req, mu, G, alph, sig, p_a, f_a, tf_nd, tsteps, ultra)

addpath ../IMRv2/src/forward_solver/

% equation options
% ------- Material Properties ------------------------%
kappa = 1.4;
T8 = 298.15;
rho8 = 1048;

% ------- Simulation settings ------------------------%
radial = 2;
vapor = 1;
collapse = 0;
bubtherm = 1;
medtherm = 0;
masstrans = 1;
stress = 2;
% vapor = 0;
% collapse = 0;
% bubtherm = 0;
% medtherm = 0;
% masstrans = 0;
% stress = 2;




if ultra
    % --------- Ultrasound settins -----------------------%
    pa = p_a;
    omega = 2*pi*f_a;
    wavetype = 4;
else
    pa = 0;
    omega = 0;
    wavetype = 0;
end

% ------ Simulation time ---------------------------- %
tc = Rmax*sqrt(rho8/101325);
tfin = tf_nd*tc;
tvector = linspace(0,tfin,tsteps);
varin = {'progdisplay',0,'radial',radial,'bubtherm',bubtherm,'tvector',tvector,...
         'vapor',vapor,'medtherm',medtherm,'masstrans',masstrans,'method',23,...
         'stress',stress,'collapse',collapse,'mu',mu,'g',G,'lambda1',0e-7,...
         'lambda2',0,'alphax', alph, 'surft', sig,'r0',Rmax,'req',Req,'kappa',kappa,'t8',T8,...
         'rho8',rho8, 'pa',pa 'omega', omega, 'wave_type', wavetype};

[t,R,Rd,~,~,~,~,Rdd] = f_imr_fd(varin{:},'Nt',75);
end