clear;
clc;
close all;

%% ================================================================
% STANDALONE PLATFORM-BALANCING TEST
%
% Workflow:
% 1) The platform evolves freely with the masses fixed.
% 2) The script automatically detects when the angular velocity remains
% small for to specified time.
% 3) From that instant, X/and balancing is enabled.
% 4) The script detects when phi and theta remain within the deadband
% and the masses are practically stationary.
%
% This script uses the following files without modification:
% platform_parameters.m
% actuator_torque_parameters.m
% platform_dynamics_ode.m
% simulate_platform.m
%% ================================================================

%% ================================================================
% PARAMETERS YOU CAN MODIFY FOR TESTING
% ================================================================

Ktheta_prueba = 0.08;
Kphi_prueba   = 0.08;

vmax_X_prueba = 0.001;
vmax_Y_prueba = 0.001;

amax_X_prueba = 5.0e-4;
amax_Y_prueba = 5.0e-4;

tau_mX_prueba = 1.0;
tau_mY_prueba = 1.0;

error_lim_xy_deg = 0.01;

T_libre_max = 2500;

omega_tol_libre = 1.0e-4;
t_hold_libre = 30;

T_balanceo_max = 2500;

v_tol_balanceo = 1.0e-6;
t_hold_balanceo = 20;

%% ================================================================
% MODEL PARAMETERS
% ================================================================

P = platform_parameters();
A = actuator_torque_parameters();

A.Ktheta = Ktheta_prueba;
A.Kphi   = Kphi_prueba;

A.vmax_X = vmax_X_prueba;
A.vmax_Y = vmax_Y_prueba;

A.amax_X = amax_X_prueba;
A.amax_Y = amax_Y_prueba;

A.tau_mX = tau_mX_prueba;
A.tau_mY = tau_mY_prueba;

A.error_lim_xy = deg2rad(error_lim_xy_deg);

%% ================================================================
% INITIAL CONDITION
% ================================================================

q0 = [0;0;0;1];
w0 = [0;0;0];

d0 = [0;0;0];
v0 = [0;0;0];

xi0 = [0;0];
Omega_RW0 = A.Omega_RW0(:);

x0 = [q0;w0;d0;v0;xi0;Omega_RW0];

cfg = struct;
cfg.modo = 'balanceo_inicial';

%% ================================================================
% HIDE INTERNAL simulate_platform PLOTS
% ================================================================

visibilidad_original = get(groot,'DefaultFigureVisible');
set(groot,'DefaultFigureVisible','off');
cleanup_visibilidad = onCleanup( ...
    @() set(groot,'DefaultFigureVisible',visibilidad_original));

%% ================================================================
% STAGE 1 - FREE EVOLUTION
% ================================================================

A_libre = A;
A_libre.T_balanceo1 = T_libre_max + 10*A.dt;

tspan_libre = 0:A.dt:T_libre_max;

fprintf('\n=============================================================\n');
fprintf(' STAGE 1 - FREE EVOLUTION\n');
fprintf('=============================================================\n');

datos_libre = simulate_platform( ...
    tspan_libre, ...
    x0, ...
    cfg, ...
    P, ...
    A_libre ...
);

close all;

omega_norm = sqrt( ...
    datos_libre.wx_true.^2 + ...
    datos_libre.wy_true.^2 + ...
    datos_libre.wz_true.^2 ...
);

mask_libre = omega_norm < omega_tol_libre;

idx_inicio_balanceo = buscar_condicion_sostenida( ...
    mask_libre, ...
    datos_libre.t, ...
    t_hold_libre ...
);

if isempty(idx_inicio_balanceo)

    warning([ ...
        'No se alcanzo el criterio de estabilizacion natural en %.1f s. ', ...
        'Balancing starts after the free-motion simulation.' ...
    ],T_libre_max);

    idx_inicio_balanceo = length(datos_libre.t);

end

t_inicio_balanceo = datos_libre.t(idx_inicio_balanceo);
x_inicio_balanceo = datos_libre.x(idx_inicio_balanceo,:)';

fprintf('balancing start: t = %.3f s\n',t_inicio_balanceo);
fprintf('phi   = %.6f deg\n',datos_libre.phi_deg(idx_inicio_balanceo));
fprintf('theta = %.6f deg\n',datos_libre.theta_deg(idx_inicio_balanceo));
fprintf('|omega| = %.6e rad/s\n',omega_norm(idx_inicio_balanceo));

%% ================================================================
% STAGE 2 - X/and BALANCING
% ================================================================

A_bal = A;
A_bal.T_balanceo1 = t_inicio_balanceo;

tf_balanceo = t_inicio_balanceo + T_balanceo_max;
tspan_balanceo = t_inicio_balanceo:A.dt:tf_balanceo;

fprintf('\n=============================================================\n');
fprintf(' STAGE 2 - X/and BALANCING\n');
fprintf('=============================================================\n');
fprintf('Ktheta = %.6g m/(s*rad)\n',A_bal.Ktheta);
fprintf('Kphi   = %.6g m/(s*rad)\n',A_bal.Kphi);
fprintf('vmax_X = %.6g m/s\n',A_bal.vmax_X);
fprintf('vmax_Y = %.6g m/s\n',A_bal.vmax_Y);
fprintf('amax_X = %.6g m/s^2\n',A_bal.amax_X);
fprintf('amax_Y = %.6g m/s^2\n',A_bal.amax_Y);

datos_bal = simulate_platform( ...
    tspan_balanceo, ...
    x_inicio_balanceo, ...
    cfg, ...
    P, ...
    A_bal ...
);

close all;

mask_balanceo = ...
    abs(datos_bal.phi)   < A_bal.error_lim_xy & ...
    abs(datos_bal.theta) < A_bal.error_lim_xy & ...
    abs(datos_bal.vX)    < v_tol_balanceo & ...
    abs(datos_bal.vY)    < v_tol_balanceo;

idx_fin_balanceo = buscar_condicion_sostenida( ...
    mask_balanceo, ...
    datos_bal.t, ...
    t_hold_balanceo ...
);

if isempty(idx_fin_balanceo)

    warning('The final balancing criterion was not reached.');
    idx_fin_balanceo = length(datos_bal.t);
    balanceo_completo = false;

else

    balanceo_completo = true;

end

t_fin_balanceo = datos_bal.t(idx_fin_balanceo);
tiempo_balanceo = t_fin_balanceo - t_inicio_balanceo;

fprintf('\n=============================================================\n');
fprintf(' result\n');
fprintf('=============================================================\n');

if balanceo_completo
    fprintf('Balancing completed.\n');
else
    fprintf('The maximum time was reached without satisfying the sustained criterion.\n');
end

fprintf('balancing time = %.3f s\n',tiempo_balanceo);
fprintf('phi final   = %.6f deg\n',datos_bal.phi_deg(idx_fin_balanceo));
fprintf('theta final = %.6f deg\n',datos_bal.theta_deg(idx_fin_balanceo));
fprintf('dX final    = %.3f mm\n',1000*datos_bal.dX(idx_fin_balanceo));
fprintf('dY final    = %.3f mm\n',1000*datos_bal.dY(idx_fin_balanceo));
fprintf('vX final    = %.6e m/s\n',datos_bal.vX(idx_fin_balanceo));
fprintf('vY final    = %.6e m/s\n',datos_bal.vY(idx_fin_balanceo));

%% ================================================================
% TRIM DATA AT THE COMPLETION TIME
% ================================================================

idxL = 1:idx_inicio_balanceo;
idxB = 1:idx_fin_balanceo;

t_total = [ ...
    datos_libre.t(idxL); ...
    datos_bal.t(idxB(2:end)) ...
];

phi_total = [ ...
    datos_libre.phi_deg(idxL); ...
    datos_bal.phi_deg(idxB(2:end)) ...
];

theta_total = [ ...
    datos_libre.theta_deg(idxL); ...
    datos_bal.theta_deg(idxB(2:end)) ...
];

dX_total = [ ...
    datos_libre.dX(idxL); ...
    datos_bal.dX(idxB(2:end)) ...
];

dY_total = [ ...
    datos_libre.dY(idxL); ...
    datos_bal.dY(idxB(2:end)) ...
];

vX_total = [ ...
    datos_libre.vX(idxL); ...
    datos_bal.vX(idxB(2:end)) ...
];

vY_total = [ ...
    datos_libre.vY(idxL); ...
    datos_bal.vY(idxB(2:end)) ...
];

%% ================================================================
% RECONSTRUCT VELOCITY COMMANDS FOR PLOTTING
% ================================================================

vX_cmd_total = zeros(size(t_total));
vY_cmd_total = zeros(size(t_total));

idx_control = t_total >= t_inicio_balanceo;

theta_rad_total = deg2rad(theta_total);
phi_rad_total = deg2rad(phi_total);

indices_control = find(idx_control);

for ii = 1:length(indices_control)

    k = indices_control(ii);

    if abs(theta_rad_total(k)) > A_bal.error_lim_xy
        vX_cmd_total(k) = -A_bal.Ktheta*theta_rad_total(k);
    end

    if abs(phi_rad_total(k)) > A_bal.error_lim_xy
        vY_cmd_total(k) = A_bal.Kphi*phi_rad_total(k);
    end

end

vX_cmd_total = max(-A_bal.vmax_X,min(A_bal.vmax_X,vX_cmd_total));
vY_cmd_total = max(-A_bal.vmax_Y,min(A_bal.vmax_Y,vY_cmd_total));

%% ================================================================
% PLOTS
% ================================================================

set(groot,'DefaultFigureVisible','on');

figure('Name','Prueba de balanceo','NumberTitle','off');

subplot(3,1,1);

plot(t_total,phi_total,'LineWidth',1.3);
hold on;
plot(t_total,theta_total,'LineWidth',1.3);

yline(+error_lim_xy_deg,'--');
yline(-error_lim_xy_deg,'--');
xline(t_inicio_balanceo,':','Inicio balanceo');

grid on;
ylabel('Angulo [deg]');
legend('\phi','\theta','Location','best');
title('Evolucion libre y balanceo X/Y');

subplot(3,1,2);

plot(t_total,1000*dX_total,'LineWidth',1.3);
hold on;
plot(t_total,1000*dY_total,'LineWidth',1.3);
xline(t_inicio_balanceo,':','Inicio balanceo');

grid on;
ylabel('Posicion [mm]');
legend('d_X','d_Y','Location','best');

subplot(3,1,3);

plot(t_total,1000*vX_total,'LineWidth',1.3);
hold on;
plot(t_total,1000*vX_cmd_total,'--','LineWidth',1.1);
plot(t_total,1000*vY_total,'LineWidth',1.3);
plot(t_total,1000*vY_cmd_total,'--','LineWidth',1.1);
xline(t_inicio_balanceo,':','Inicio balanceo');

grid on;
xlabel('Tiempo [s]');
ylabel('Velocidad [mm/s]');
legend('v_X','v_{X,cmd}','v_Y','v_{Y,cmd}','Location','best');

%% ================================================================
% RESULTS FOR COMPARISON ACROSS TESTS
% ================================================================

RESULTADO_BALANCEO = struct;

RESULTADO_BALANCEO.Ktheta = A_bal.Ktheta;
RESULTADO_BALANCEO.Kphi = A_bal.Kphi;

RESULTADO_BALANCEO.vmax_X = A_bal.vmax_X;
RESULTADO_BALANCEO.vmax_Y = A_bal.vmax_Y;

RESULTADO_BALANCEO.amax_X = A_bal.amax_X;
RESULTADO_BALANCEO.amax_Y = A_bal.amax_Y;

RESULTADO_BALANCEO.t_inicio_balanceo = t_inicio_balanceo;
RESULTADO_BALANCEO.t_fin_balanceo = t_fin_balanceo;
RESULTADO_BALANCEO.tiempo_balanceo = tiempo_balanceo;
RESULTADO_BALANCEO.balanceo_completo = balanceo_completo;

RESULTADO_BALANCEO.phi_final_deg = datos_bal.phi_deg(idx_fin_balanceo);
RESULTADO_BALANCEO.theta_final_deg = datos_bal.theta_deg(idx_fin_balanceo);

RESULTADO_BALANCEO.dX_final = datos_bal.dX(idx_fin_balanceo);
RESULTADO_BALANCEO.dY_final = datos_bal.dY(idx_fin_balanceo);

RESULTADO_BALANCEO.datos_libre = datos_libre;
RESULTADO_BALANCEO.datos_balanceo = datos_bal;

%% ================================================================
% LOCAL FUNCTION
% ================================================================

function idx = buscar_condicion_sostenida(mask,t,t_hold)

    mask = logical(mask(:));
    t = t(:);

    idx = [];

    if isempty(mask) || length(mask) ~= length(t)
        return;
    end

    i_ini = 1;

    while i_ini <= length(t)

        if ~mask(i_ini)
            i_ini = i_ini + 1;
            continue;
        end

        i_fin = i_ini;

        while i_fin < length(t) && mask(i_fin+1)
            i_fin = i_fin + 1;
        end

        if t(i_fin)-t(i_ini) >= t_hold
            idx = i_ini;
            return;
        end

        i_ini = i_fin + 1;

    end

end
