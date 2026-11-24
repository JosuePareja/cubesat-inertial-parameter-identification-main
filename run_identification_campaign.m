function secuencia = run_identification_campaign(P,A,etiqueta)
%RUN_IDENTIFICATION_CAMPAIGN
% Runs to complete identification campaign on plant P.
%
% The same function is used for:
% 1) the BASE platform;
% 2) the BASE + CubeSat platform.
%
% The difference between both campaigns is contained in P and in
% to.rCMz_obj. All final outputs are reconstructed in the
% canonical configuration dX=dY=dZ=0 so the campaigns can be subtracted.

    if nargin < 3 || isempty(etiqueta)
        etiqueta = 'CAMPAIGN';
    end

    if nargin < 1 || isempty(P)
        P = platform_parameters();
    end

    if nargin < 2 || isempty(A)
        A = actuator_torque_parameters();
    end

    etiqueta = char(string(etiqueta));
    A = configurar_archivos_campana(A,etiqueta);
    P.etiqueta_campana = etiqueta;

    fprintf('\n\n');
    fprintf('#############################################################\n');
    fprintf(' IDENTIFICATION CAMPAIGN: %s\n',etiqueta);
    fprintf(' M_total = %.6f kg\n',P.M_total);
    fprintf(' rCMz target before 7P estimation = %.3f mm\n', ...
        1000*A.rCMz_obj);
    fprintf('#############################################################\n');

%% ================================================================
% INITIAL CONDITIONS
% ================================================================

q0 = [0;0;0;1];
w0 = [0;0;0];
d0 = [0;0;0];
v0 = [0;0;0];

% X/and integral states. They are reused by stage and reset
% explicitly when switching between angular PI and position PI control.
xi0 = [0;0];

% Initial relative reaction-wheel speeds.
Omega_RW0 = A.Omega_RW0(:);

x0 = [q0;w0;d0;v0;xi0;Omega_RW0];


%% ================================================================
% STEP 1A - FREE OSCILLATION + INITIAL X/and BALANCING
% ================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' STEP 1A - INITIAL X/and BALANCING\n');
fprintf('=============================================================\n');

cfg = struct;
cfg.modo = 'balanceo_inicial';

tspan = 0:A.dt:A.T_exc1;

datos_balanceo1 = ...
    simulate_platform(tspan,x0,cfg,P,A);

x_bal1 = datos_balanceo1.x(end,:)';
x_bal1(14:15) = 0;

dX_bal1 = x_bal1(8);
dY_bal1 = x_bal1(9);
dZ_bal1 = x_bal1(10);

fprintf('balanced dX = %.3f mm\n',1000*dX_bal1);
fprintf('balanced dY = %.3f mm\n',1000*dY_bal1);
fprintf('dZ            = %.3f mm\n',1000*dZ_bal1);


%% ================================================================
% STEP 1B - FIRST EXCITATION
% ================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' STEP 1B - FIRST EXCITATION\n');
fprintf('=============================================================\n');

cfg = struct;

cfg.modo = 'excitacion_1';

cfg.dX_bal = dX_bal1;
cfg.dY_bal = dY_bal1;
cfg.dZ_bal = dZ_bal1;

tspan = A.T_exc1:A.dt:A.T_fin_exc1;

datos1 = ...
    simulate_platform(tspan,x_bal1,cfg,P,A);

save(A.archivo_etapa1,'datos1');

fprintf('Data saved to %s\n',A.archivo_etapa1);


%% ================================================================
% STEP 1C - FIRST 10-PARAMETER ESTIMATION
% ================================================================

R1 = estimate_platform_parameters(datos1,P,A);

% In to true test, only the noisy estimate would be available.
% Therefore, the physical sequence is determined from this estimate.
rCM_est = R1.ruido.rCM_hat_ini;

fprintf('\nCenter of mass used to continue the sequence:\n');
fprintf('rx = %.4f mm\n',1000*rCM_est(1));
fprintf('ry = %.4f mm\n',1000*rCM_est(2));
fprintf('rz = %.4f mm\n',1000*rCM_est(3));


%% ================================================================
% STEP 2 - X/and REBALANCING
% ================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' STEP 2 - X/and REBALANCING\n');
fprintf('=============================================================\n');

x_inicio_reb = datos1.x(end,:)';
x_inicio_reb(14:15) = 0;

t0_reb = datos1.t(end);
tf_reb = t0_reb + A.T_rebalanceo_max;

cfg = struct;

cfg.modo = 'rebalanceo_xy';
cfg.dZ_hold = x_inicio_reb(10);

tspan = t0_reb:A.dt:tf_reb;

datos_reb = ...
    simulate_platform(tspan,x_inicio_reb,cfg,P,A);

mask_xy = ...
    abs(datos_reb.phi) < A.error_lim_xy ...
    & abs(datos_reb.theta) < A.error_lim_xy ...
    & abs(datos_reb.vX) < A.v_tol ...
    & abs(datos_reb.vY) < A.v_tol;

idx_reb = buscar_estabilidad( ...
    mask_xy,datos_reb.t,A.t_hold_xy);

if isempty(idx_reb)

    warning('The rebalancing criterion was not reached within the maximum time.');
    idx_reb = length(datos_reb.t);

end

x_reb = datos_reb.x(idx_reb,:)';
x_reb(14:15) = 0;
t_reb = datos_reb.t(idx_reb);

fprintf('Rebalancing completed at t = %.2f s\n',t_reb);
fprintf('phi   = %.5f deg\n',datos_reb.phi_deg(idx_reb));
fprintf('theta = %.5f deg\n',datos_reb.theta_deg(idx_reb));
fprintf('dX    = %.3f mm\n',1000*x_reb(8));
fprintf('dY    = %.3f mm\n',1000*x_reb(9));


%% ================================================================
% STEP 3 - Z-AXIS CENTER-OF-MASS COMPENSATION
% ================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' STEP 3 - Z-AXIS COMPENSATION\n');
fprintf('=============================================================\n');

% The estimated rCM,z corresponds to the original configuration.
% X/and displacements do not change the CM Z coordinate.
%
% NOTE:
% The target Z-axis CM position depends on the campaign.
% BASE uses to.rCMz_obj_base; BASE+CubeSat uses to.rCMz_obj_cubesat.
% The current travel is limited by to.dZ_min/to.dZ_max.
%
% rCMz_nuevo =
% rCMz_estimado + (mZ/M_total)*dZ
%
% Therefore:
%
% dZ_target = (M_total/mZ)*(rCMz_obj-rCMz_estimado)

dZ_target_raw = ...
    (P.M_total/P.mZ)*(A.rCMz_obj-rCM_est(3));

dZ_target = max( ...
    A.dZ_min, ...
    min(A.dZ_max,dZ_target_raw) ...
);

fprintf('initial estimated rCMz = %.4f mm\n',1000*rCM_est(3));
fprintf('target rCMz         = %.4f mm\n',1000*A.rCMz_obj);
fprintf('computed dZ          = %.3f mm\n',1000*dZ_target_raw);
fprintf('applied dZ           = %.3f mm\n',1000*dZ_target);
fprintf('Configured maximum Z travel = %.3f mm\n',1000*P.dZ_max_abs);
if abs(dZ_target) > P.dZ_max_abs
    fprintf(['NOTE: this compensation requires to travel larger than the ', ...
        'recorrido maximo configurado. Se aplicara saturacion.\n']);
end

if abs(dZ_target-dZ_target_raw) > 1e-12

    warning(['The required Z displacement exceeds the allowed travel. ' ...
        'Se aplicara el limite del actuador.']);

end

cfg = struct;

cfg.modo = 'compensacion_z';
cfg.dZ_target = dZ_target;

tf_comp = t_reb + A.T_compZ_max;

tspan = t_reb:A.dt:tf_comp;

datos_comp = ...
    simulate_platform(tspan,x_reb,cfg,P,A);

% CM,z predicted by the FIRST ESTIMATION during Z-axis motion
rCMz_est_series = ...
    rCM_est(3) ...
    + (P.mZ/P.M_total).*datos_comp.dZ;

mask_z = ...
    abs(rCMz_est_series) < A.rCMz_lim ...
    & abs(datos_comp.dZ-dZ_target) < A.pos_tol_Z ...
    & abs(datos_comp.vZ) < A.v_tol ...
    & abs(datos_comp.phi) < A.error_lim_xy ...
    & abs(datos_comp.theta) < A.error_lim_xy;

idx_comp = buscar_estabilidad( ...
    mask_z,datos_comp.t,A.t_hold_Z);

if isempty(idx_comp)

    warning('The Z-compensation criterion was not fully reached.');
    idx_comp = length(datos_comp.t);

end

x_comp = datos_comp.x(idx_comp,:)';
x_comp(14:15) = 0;
t_comp = datos_comp.t(idx_comp);

rCMz_est_final = ...
    rCM_est(3) ...
    + (P.mZ/P.M_total)*x_comp(10);

fprintf('\nCompensation completed at t = %.2f s\n',t_comp);
fprintf('final dZ               = %.3f mm\n',1000*x_comp(10));
fprintf('final estimated rCMz    = %.4f mm\n',1000*rCMz_est_final);
fprintf('TRUE simulated rCMz   = %.4f mm\n', ...
    1000*datos_comp.rCM_real(idx_comp,3));

if abs(rCMz_est_final) < A.rCMz_lim

    fprintf('CRITERION |rCMz| < %.2f mm: MET\n', ...
        1000*A.rCMz_lim);

else

    fprintf('CRITERION |rCMz| < %.2f mm: NOT MET\n', ...
        1000*A.rCMz_lim);

end


%% ================================================================
% STEP 4 - WAIT UNTIL THE MASSES ARE STATIONARY
% ================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' STEP 4 - WAIT / SETTLING\n');
fprintf('=============================================================\n');

cfg = struct;

cfg.modo = 'espera';

cfg.dX_hold = x_comp(8);
cfg.dY_hold = x_comp(9);
cfg.dZ_hold = x_comp(10);

tspan = t_comp:A.dt:(t_comp+A.T_espera);

datos_espera = ...
    simulate_platform(tspan,x_comp,cfg,P,A);

x_quieto = datos_espera.x(end,:)';
x_quieto(14:15) = 0;

fprintf('Mass velocities at the end of the settling interval:\n');
fprintf('vX = %.6e m/s\n',x_quieto(11));
fprintf('vY = %.6e m/s\n',x_quieto(12));
fprintf('vZ = %.6e m/s\n',x_quieto(13));

% Represents the mechanical lock once the masses are considered
% fully stationary.
if max(abs(x_quieto(11:13))) <= 10*A.v_tol

    x_quieto(11:13) = 0;

else

    warning(['The mass velocities are still not sufficiently small. ' ...
        'Considere aumentar A.T_espera.']);

end


%% ================================================================
% KNOWN Mr FOR THE SECOND ESTIMATION
% ================================================================

% The second excitation starts with the PI integrators reset.
% This avoids retaining numerical memory from previous stages.
x_quieto(14:15) = 0;

d_ref2 = x_quieto(8:10);

% Alias retained for compatibility with previous results.
d_final = d_ref2;

DeltaMr_final = ...
    P.mX*P.eX*d_final(1) ...
    + P.mY*P.eY*d_final(2) ...
    + P.mZ*P.eZ*d_final(3);

Mr_conocido = struct;

Mr_conocido.ideal = ...
    R1.ideal.Mr_hat_ini + DeltaMr_final;

Mr_conocido.ruido = ...
    R1.ruido.Mr_hat_ini + DeltaMr_final;

Mr_conocido.real = ...
    P.Mr0 + DeltaMr_final;

fprintf('\nMr of the compensated configuration [kg*m]:\n');

fprintf('TRUE:\n');
disp(Mr_conocido.real);

fprintf('Used by the ideal estimator:\n');
disp(Mr_conocido.ideal);

fprintf('Used by the NOISY estimator:\n');
disp(Mr_conocido.ruido);


%% ================================================================
% STEP 5 - SECOND EXCITATION
%
% omega_x : fast and mass
% omega_y : fast X mass
% omega_z : Z reaction wheel
%
% The estimator knows d(t), v(t), Mr(t), DeltaJ(t), h_rel(t), and
% the angular momentum stored in the wheels, h_RW(t).
% ================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' STEP 5 - FAST X/and MASSES + Z WHEEL WITH h_RW\n');
fprintf('=============================================================\n');

t0_2 = datos_espera.t(end);

cfg = struct;

cfg.modo = 'excitacion_2';
cfg.t0 = t0_2;

% Reference configuration for the second estimation.
cfg.dX_ref2 = d_ref2(1);
cfg.dY_ref2 = d_ref2(2);
cfg.dZ_ref2 = d_ref2(3);

tspan = ...
    t0_2:A.dt:(t0_2+A.T_exc2_total);

datos2 = ...
    simulate_platform(tspan,x_quieto,cfg,P,A);

fprintf('\nWheel speed at the start of stage 2 [rpm]:\n');
disp(datos2.Omega_RW_rpm(1,:));
fprintf('Wheel speed at the end of stage 2 [rpm]:\n');
disp(datos2.Omega_RW_rpm(end,:));

datos2.Mr_conocido = Mr_conocido;
datos2.d_ref2 = d_ref2;
datos2.d_final = d_ref2;  % alias

save(A.archivo_etapa2,'datos2');

fprintf('\nZ wheel:\n');
fprintf('  +%.3f mNm: t_rel = %.1f to %.1f s\n', ...
    1e3*A.tauZ2,A.tZ2,A.tZ2+A.DT_Z2);
fprintf('  -%.3f mNm: t_rel = %.1f to %.1f s\n', ...
    1e3*A.tauZ2,A.tZ2_opuesto,A.tZ2_opuesto+A.DT_Z2);

fprintf('\nY mass -> omega_x excitation:\n');
fprintf('  amplitude = +/- %.1f mm\n',1000*A.deltaY2);
fprintf('  +delta: %.1f to %.1f s\n',A.tY2_pos_ini,A.tY2_pos_fin);
fprintf('  -delta: %.1f to %.1f s\n',A.tY2_neg_ini,A.tY2_neg_fin);

fprintf('\nX mass -> omega_y excitation:\n');
fprintf('  amplitude = +/- %.1f mm\n',1000*A.deltaX2);
fprintf('  +delta: %.1f to %.1f s\n',A.tX2_pos_ini,A.tX2_pos_fin);
fprintf('  -delta: %.1f to %.1f s\n',A.tX2_neg_ini,A.tX2_neg_fin);

fprintf('\nMaximum mass speed in stage 2:\n');
fprintf('  X = %.1f mm/s\n',1000*A.vmax_X2);
fprintf('  Y = %.1f mm/s\n',1000*A.vmax_Y2);


%% ================================================================
% STEP 5B - SAVITZKY-GOLAY FILTERING
%
% IMPORTANT:
% datos2 retains the original noisy measurements.
% datos2_filtrados retains the same measurements after the
% Savitzky-Golay filter.
%
% Which of the two structures is passed to the second estimator depends on:
%
% to.modo_datos_est2
%
% defined in actuator_torque_parameters.m.
%
% Filtered signals:
% wx, wy, wz
% gx, gy, gz
%
% The ideal signals, Omega_RW, and diagnostic torques remain
% unchanged. Omega_RW is treated as to known measurement from the encoder/model.
% ================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' STEP 5B - SAVITZKY-GOLAY FILTER\n');
fprintf('=============================================================\n');

datos2_filtrados = filtrar_datos_sg(datos2,A);

save(A.archivo_etapa2_filtrado,'datos2_filtrados');

fprintf('SG order             = %d\n',A.SG_orden);
fprintf('SG window           = %d samples\n',A.SG_ventana);
fprintf('Window duration  = %.3f s\n',A.SG_ventana*A.dt);
fprintf('Filtered data saved to %s\n',A.archivo_etapa2_filtrado);

graficar_comparacion_sg(datos2,datos2_filtrados);


%% ================================================================
% STEP 6 - DATA SELECTION AND SECOND ESTIMATION
%
% To switch between noisy and filtered data, only modify:
%
% to.modo_datos_est2 = 'noisy';
%
% or:
%
% to.modo_datos_est2 = 'filtered';
%
% in actuator_torque_parameters.m.
% ================================================================

modo_est2 = lower(strtrim(A.modo_datos_est2));

switch modo_est2

    case 'ruidosa'

        datos2_estimador = datos2;

 % Metadata explicitly recording that not SG filter was applied.
        datos2_estimador.filtrado_sg = struct;
        datos2_estimador.filtrado_sg.aplicado = false;
        datos2_estimador.filtrado_sg.orden = [];
        datos2_estimador.filtrado_sg.ventana = [];
        datos2_estimador.filtrado_sg.dt = A.dt;
        datos2_estimador.filtrado_sg.ventana_s = [];
        datos2_estimador.filtrado_sg.gravedad_filtrada = false;

        fprintf('\nSECOND ESTIMATION: UNFILTERED NOISY data will be used.\n');

    case 'filtrada'

        datos2_estimador = datos2_filtrados;

        fprintf('\nSECOND ESTIMATION: NOISY + SAVITZKY-GOLAY data will be used.\n');

    otherwise

        error('to.modo_datos_est2 must be ''noisy'' or ''filtered''.');

end

R2 = estimate_inertia_parameters( ...
    datos2_estimador,P,A,Mr_conocido,d_ref2);


%% ================================================================
% final SUMMARY
% ================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' final CAMPAIGN SUMMARY: %s\n',etiqueta);
fprintf('=============================================================\n');

fprintf('\nData mode used in the second estimation: %s\n', ...
    upper(modo_est2));

fprintf('\nCANONICAL CONFIGURATION FOR CAMPAIGN COMPARISON:\n');
fprintf('  dX = dY = dZ = 0\n');

fprintf('\nCanonical center of mass estimated with noise - 10P [mm]:\n');
disp(1000*R1.canonica.ruido.rCM);

fprintf('TRUE compensated center of mass at the start of stage 2 [mm]:\n');
disp(1000*datos2.rCM_real(1,:)');

fprintf('\nTRUE CANONICAL J [kg*m^2]:\n');
disp(R2.canonica.real.J);

fprintf('ESTIMATED CANONICAL J - ideal second stage [kg*m^2]:\n');
disp(R2.canonica.ideal.J);

fprintf('ESTIMATED CANONICAL J - MEASURED second stage [kg*m^2]:\n');
disp(R2.canonica.medida.J);

fprintf('\nESTIMATED CANONICAL S - NOISY 10P [kg*m]:\n');
disp(R1.canonica.ruido.S);

fprintf('CANONICAL S PROPAGATED TO THE 7P OUTPUT [kg*m]:\n');
disp(R2.canonica.medida.S);

omega_pico_ideal = [ ...
    max(abs(datos2.wx_true)), ...
    max(abs(datos2.wy_true)), ...
    max(abs(datos2.wz_true)) ...
];

omega_pico_medida = [ ...
    max(abs(datos2.wx)), ...
    max(abs(datos2.wy)), ...
    max(abs(datos2.wz)) ...
];

fprintf('\nAngular-velocity peaks - ideal [rad/s]:\n');
fprintf('  |wx|max = %.6f\n',omega_pico_ideal(1));
fprintf('  |wy|max = %.6f\n',omega_pico_ideal(2));
fprintf('  |wz|max = %.6f\n',omega_pico_ideal(3));

fprintf('\nAngular-velocity peaks - NOISY [rad/s]:\n');
fprintf('  |wx|max = %.6f\n',omega_pico_medida(1));
fprintf('  |wy|max = %.6f\n',omega_pico_medida(2));
fprintf('  |wz|max = %.6f\n',omega_pico_medida(3));

fprintf('\nApproximate X/Y target = %.3f rad/s\n',A.omega_obj_exc2);

fprintf('=============================================================\n');


%% ================================================================
% SECOND-EXCITATION PLOTS
% ================================================================

trel2 = datos2.t-datos2.t(1);

% ------------------------------------------------------------
% velocity angular
% ------------------------------------------------------------

figure( ...
    'Name','Segunda excitacion - velocidad angular', ...
    'NumberTitle','off' ...
);

subplot(3,1,1);
plot(trel2,datos2.wx_true,'LineWidth',1.3);
hold on;
plot(trel2,datos2.wx,'LineWidth',0.5);
grid on;
ylabel('\omega_x [rad/s]');
title('Segunda etapa: respuesta angular');
legend('Ideal','Ruidosa','Location','best');

subplot(3,1,2);
plot(trel2,datos2.wy_true,'LineWidth',1.3);
hold on;
plot(trel2,datos2.wy,'LineWidth',0.5);
grid on;
ylabel('\omega_y [rad/s]');

subplot(3,1,3);
plot(trel2,datos2.wz_true,'LineWidth',1.3);
hold on;
plot(trel2,datos2.wz,'LineWidth',0.5);
grid on;
ylabel('\omega_z [rad/s]');
xlabel('Tiempo relativo [s]');

% ------------------------------------------------------------
% position and velocity of the masses X/and
% ------------------------------------------------------------

figure( ...
    'Name','Segunda excitacion - masas X Y', ...
    'NumberTitle','off' ...
);

subplot(2,1,1);
plot(trel2,1000*(datos2.dX-d_ref2(1)),'LineWidth',1.3);
hold on;
plot(trel2,1000*(datos2.dY-d_ref2(2)),'LineWidth',1.3);
grid on;
ylabel('\Delta d [mm]');
title('Movimiento relativo de las masas durante la segunda etapa');
legend('\Delta d_X','\Delta d_Y','Location','best');

subplot(2,1,2);
plot(trel2,1000*datos2.vX,'LineWidth',1.3);
hold on;
plot(trel2,1000*datos2.vY,'LineWidth',1.3);
grid on;
ylabel('v [mm/s]');
xlabel('Tiempo relativo [s]');
legend('v_X','v_Y','Location','best');

% ------------------------------------------------------------
% wheel Z and torques gravitational incremental of the masses
% ------------------------------------------------------------

figure( ...
    'Name','Second excitation - equivalent torques', ...
    'NumberTitle','off' ...
);

subplot(3,1,1);
plot(trel2,1e3*datos2.tau_masaY(:,1),'LineWidth',1.3);
grid on;
ylabel('\tau_X [mNm]');
title('Torque gravitacional incremental: masa Y -> eje X');

subplot(3,1,2);
plot(trel2,1e3*datos2.tau_masaX(:,2),'LineWidth',1.3);
grid on;
ylabel('\tau_Y [mNm]');
title('Torque gravitacional incremental: masa X -> eje Y');

subplot(3,1,3);
plot(trel2,1e3*datos2.tau_RW(:,3),'LineWidth',1.3);
grid on;
ylabel('\tau_Z [mNm]');
xlabel('Tiempo relativo [s]');
title('Rueda de reaccion Z');


%% ================================================================
% complete history of TODA the sequence
%
% Concatenate only the segments that were actually used.
% for continuar the sequence. this es importante because the stages
% of rebalancing and compensation simulan until to time maximum,
% but the sequence continues from the instant in that reaches the
% criterion of stability.
% ================================================================

historico = construir_historico_completo( ...
    datos_balanceo1, ...
    datos1, ...
    datos_reb,idx_reb, ...
    datos_comp,idx_comp, ...
    datos_espera, ...
    datos2);


%% ================================================================
% save TODO
% ================================================================

secuencia = struct;

secuencia.R1 = R1;
secuencia.R2 = R2;

% results estandarizados for futuras comparisons BASE / BASE+CubeSat.
% both estructuras estan expressed in the same configuration:
%
% dX = dY = dZ = 0
%
% the output 10P provides principalmente S and rCM.
% the output 7P provides principalmente J; its S es the value of 10P
% transported algebraically to the configuration canonical.
secuencia.configuracion_canonica = R2.configuracion_canonica;
secuencia.canonica_10P = R1.canonica;
secuencia.canonica_7P = R2.canonica;

% Recommended canonical result for future CubeSat identification:
% - S and rCM come of the estimation 10P.
% - J comes of the estimation 7P.
% - b is retained as a parameter of the complete test system.
secuencia.canonica_final = struct;

secuencia.canonica_final.ideal = struct;
secuencia.canonica_final.ideal.J = R2.canonica.ideal.J;
secuencia.canonica_final.ideal.S = R1.canonica.ideal.S;
secuencia.canonica_final.ideal.Mr = R1.canonica.ideal.S;
secuencia.canonica_final.ideal.rCM = R1.canonica.ideal.rCM;
secuencia.canonica_final.ideal.b = R2.canonica.ideal.b;

secuencia.canonica_final.medida = struct;
secuencia.canonica_final.medida.J = R2.canonica.medida.J;
secuencia.canonica_final.medida.S = R1.canonica.ruido.S;
secuencia.canonica_final.medida.Mr = R1.canonica.ruido.S;
secuencia.canonica_final.medida.rCM = R1.canonica.ruido.rCM;
secuencia.canonica_final.medida.b = R2.canonica.medida.b;

secuencia.canonica_final.real = struct;
secuencia.canonica_final.real.J = R2.canonica.real.J;
secuencia.canonica_final.real.S = R1.canonica.real.S;
secuencia.canonica_final.real.Mr = R1.canonica.real.S;
secuencia.canonica_final.real.rCM = R1.canonica.real.rCM;
secuencia.canonica_final.real.b = R2.canonica.real.b;

secuencia.canonica_final.fuente_J = 'Estimacion 7P';
secuencia.canonica_final.fuente_S = 'Estimacion 10P';
secuencia.canonica_final.d = [0;0;0];

secuencia.datos2_filtrados = datos2_filtrados;
secuencia.datos2_estimador = datos2_estimador;
secuencia.modo_datos_est2 = modo_est2;
secuencia.historico = historico;

secuencia.d_ref2 = d_ref2;
secuencia.d_final = d_ref2;  % alias
secuencia.Mr_conocido = Mr_conocido;

secuencia.omega_pico_ideal = omega_pico_ideal;
secuencia.omega_pico_medida = omega_pico_medida;

secuencia.rCM_est_inicial = rCM_est;
secuencia.rCMz_est_final = rCMz_est_final;

secuencia.dZ_target = dZ_target;

secuencia.t_rebalanceo = t_reb;
secuencia.t_compensacion = t_comp;
secuencia.t_inicio_etapa2 = t0_2;

secuencia.etiqueta = etiqueta;
secuencia.P = P;
secuencia.A = A;
secuencia.rCMz_objetivo = A.rCMz_obj;

save(A.archivo_secuencia,'secuencia');

fprintf('\nGlobal result saved to:\n%s\n',A.archivo_secuencia);


%% ================================================================
% PLOTEO final - historical complete
% ================================================================

graficar_historico_completo(historico);

% Last figure of the sequence: relative rotor speeds
% of the three wheels of reaction, expressed in rpm.
graficar_velocidades_ruedas(historico);



end


function A = configurar_archivos_campana(A,etiqueta)
%CONFIGURAR_ARCHIVOS_CAMPANA avoids that to campaign sobreescriba other.

    prefijo = upper(regexprep(etiqueta,'[^A-Za-z0-9]+','_'));
    prefijo = regexprep(prefijo,'^_+|_+$','');

    if isempty(prefijo)
        prefijo = 'CAMPAIGN';
    end

    campos = { ...
        'archivo_etapa1', ...
        'archivo_resultado1', ...
        'archivo_etapa2', ...
        'archivo_etapa2_filtrado', ...
        'archivo_resultado2', ...
        'archivo_secuencia' ...
    };

    for k = 1:numel(campos)
        c = campos{k};
        if isfield(A,c) && ~isempty(A.(c))
            [ruta,nombre,ext] = fileparts(A.(c));
            nombre_nuevo = [prefijo '_' nombre ext];
            if isempty(ruta)
                A.(c) = nombre_nuevo;
            else
                A.(c) = fullfile(ruta,nombre_nuevo);
            end
        end
    end
end

%% ========================================================================
% LOCAL FUNCTION
% ========================================================================

function idx = buscar_estabilidad(mask,t,hold_time)
% Devuelve the indice final of the first interval continuo in the that
% mask remains true during to the menos hold_time.

    idx = [];

    if isempty(mask)
        return;
    end

    inicio = [];

    for k = 1:length(mask)

        if mask(k)

            if isempty(inicio)
                inicio = k;
            end

            if t(k)-t(inicio) >= hold_time
                idx = k;
                return;
            end

        else

            inicio = [];

        end

    end

end


function datos_f = filtrar_datos_sg(datos,A)
%FILTRAR_DATA_SG
% Applies Savitzky-Golay filtering to the second-stage measurements.
%
% the target es that estimate_inertia_parameters.m not receive directamente the
% unprocessed noisy measurements.
%
% the signals ideales keep without modify for conservar
% the ideal simulation reference.

    if exist('sgolayfilt','file') ~= 2
        error([ ...
            'No se encontro sgolayfilt. ' ...
            'Signal Processing Toolbox is required to apply ' ...
            'el filtro Savitzky-Golay.' ...
        ]);
    end

    N = length(datos.t);

    if A.SG_ventana > N
        error('The Savitzky-Golay window is larger than the number of samples.');
    end

    if mod(A.SG_ventana,2) == 0
        error('The Savitzky-Golay window must be odd.');
    end

    if A.SG_ventana <= A.SG_orden
        error('The Savitzky-Golay window must be larger than the order.');
    end

    datos_f = datos;

 % Save a copy of the original measurements inside the
 % estructura filtered, for poder compararlas posteriormente.
    datos_f.wx_raw = datos.wx;
    datos_f.wy_raw = datos.wy;
    datos_f.wz_raw = datos.wz;

    datos_f.gx_raw = datos.gx;
    datos_f.gy_raw = datos.gy;
    datos_f.gz_raw = datos.gz;

 % ------------------------------------------------------------
 % GIROSCOPIO
 % ------------------------------------------------------------

    datos_f.wx = sgolayfilt( ...
        datos.wx,A.SG_orden,A.SG_ventana);

    datos_f.wy = sgolayfilt( ...
        datos.wy,A.SG_orden,A.SG_ventana);

    datos_f.wz = sgolayfilt( ...
        datos.wz,A.SG_orden,A.SG_ventana);


 % ------------------------------------------------------------
 % GRAVITY VECTOR
 % ------------------------------------------------------------

    if A.SG_filtrar_gravedad

        datos_f.gx = sgolayfilt( ...
            datos.gx,A.SG_orden,A.SG_ventana);

        datos_f.gy = sgolayfilt( ...
            datos.gy,A.SG_orden,A.SG_ventana);

        datos_f.gz = sgolayfilt( ...
            datos.gz,A.SG_orden,A.SG_ventana);

    else

        datos_f.gx = datos.gx;
        datos_f.gy = datos.gy;
        datos_f.gz = datos.gz;

    end

    datos_f.gB = [ ...
        datos_f.gx, ...
        datos_f.gy, ...
        datos_f.gz ...
    ];


 % ------------------------------------------------------------
 % METADATA of the filtered
 % ------------------------------------------------------------

    datos_f.filtrado_sg = struct;

    datos_f.filtrado_sg.aplicado = true;
    datos_f.filtrado_sg.orden = A.SG_orden;
    datos_f.filtrado_sg.ventana = A.SG_ventana;
    datos_f.filtrado_sg.dt = A.dt;
    datos_f.filtrado_sg.ventana_s = A.SG_ventana*A.dt;
    datos_f.filtrado_sg.gravedad_filtrada = ...
        A.SG_filtrar_gravedad;

end


function graficar_comparacion_sg(datos,datos_f)
%PLOTR_COMPARACION_SG
% Compara:
% 1) signal ideal
% 2) medicion noisy
% 3) medicion filtered with Savitzky-Golay

    trel = datos.t-datos.t(1);

 % ------------------------------------------------------------
 % angular VELOCITY
 % ------------------------------------------------------------

    figure( ...
        'Name','Savitzky-Golay - velocidad angular', ...
        'NumberTitle','off' ...
    );

    subplot(3,1,1);

    plot(trel,datos.wx_true,'LineWidth',1.5);
    hold on;
    plot(trel,datos.wx,'LineWidth',0.5);
    plot(trel,datos_f.wx,'LineWidth',1.2);

    grid on;

    ylabel('\omega_x [rad/s]');

    title([ ...
        'Segunda etapa: ideal, ruidosa y filtrada ', ...
        'con Savitzky-Golay' ...
    ]);

    legend( ...
        'Ideal', ...
        'Ruidosa', ...
        'Filtrada SG', ...
        'Location','best' ...
    );


    subplot(3,1,2);

    plot(trel,datos.wy_true,'LineWidth',1.5);
    hold on;
    plot(trel,datos.wy,'LineWidth',0.5);
    plot(trel,datos_f.wy,'LineWidth',1.2);

    grid on;

    ylabel('\omega_y [rad/s]');


    subplot(3,1,3);

    plot(trel,datos.wz_true,'LineWidth',1.5);
    hold on;
    plot(trel,datos.wz,'LineWidth',0.5);
    plot(trel,datos_f.wz,'LineWidth',1.2);

    grid on;

    ylabel('\omega_z [rad/s]');
    xlabel('Tiempo relativo [s]');


 % ------------------------------------------------------------
 % GRAVITY VECTOR
 % ------------------------------------------------------------

    figure( ...
        'Name','Savitzky-Golay - vector de gravedad', ...
        'NumberTitle','off' ...
    );

    subplot(3,1,1);

    plot(trel,datos.gx_true,'LineWidth',1.5);
    hold on;
    plot(trel,datos.gx,'LineWidth',0.5);
    plot(trel,datos_f.gx,'LineWidth',1.2);

    grid on;

    ylabel('g_x [m/s^2]');

    title([ ...
        'Vector de gravedad: ideal, ruidoso y filtrado ', ...
        'con Savitzky-Golay' ...
    ]);

    legend( ...
        'Ideal', ...
        'Ruidosa', ...
        'Filtrada SG', ...
        'Location','best' ...
    );


    subplot(3,1,2);

    plot(trel,datos.gy_true,'LineWidth',1.5);
    hold on;
    plot(trel,datos.gy,'LineWidth',0.5);
    plot(trel,datos_f.gy,'LineWidth',1.2);

    grid on;

    ylabel('g_y [m/s^2]');


    subplot(3,1,3);

    plot(trel,datos.gz_true,'LineWidth',1.5);
    hold on;
    plot(trel,datos.gz,'LineWidth',0.5);
    plot(trel,datos_f.gz,'LineWidth',1.2);

    grid on;

    ylabel('g_z [m/s^2]');
    xlabel('Tiempo relativo [s]');

end


function H = construir_historico_completo( ...
    datos_balanceo1,datos1, ...
    datos_reb,idx_reb, ...
    datos_comp,idx_comp, ...
    datos_espera,datos2)
%CONSTRUIR_HISTORICO_COMPLETO
% concatenates cronologicamente all the stages used actually.
%
% removes the first shows of each stage subsequent to the first
% for not duplicate the instant of junction between simulations.

    H = struct;

    campos_vector = { ...
        't', ...
        'wx_true','wy_true','wz_true', ...
        'wx','wy','wz', ...
        'gx_true','gy_true','gz_true', ...
        'gx','gy','gz', ...
        'dX','dY','dZ', ...
        'vX','vY','vZ' ...
    };

    for k = 1:numel(campos_vector)
        H.(campos_vector{k}) = [];
    end

    H.tau_RW = zeros(0,3);
    H.tau_masaX = zeros(0,3);
    H.tau_masaY = zeros(0,3);
    H.tau_exc = zeros(0,3);
    H.Omega_RW = zeros(0,3);
    H.h_RW = zeros(0,3);

 % stage 1A complete.
    H = anexar_etapa(H,datos_balanceo1,1,length(datos_balanceo1.t));

 % stage 1B complete.
    H = anexar_etapa(H,datos1,2,length(datos1.t));

 % Rebalancing: only up to the instant actually used.
    H = anexar_etapa(H,datos_reb,2,idx_reb);

 % Z compensation: only up to the instant actually used.
    H = anexar_etapa(H,datos_comp,2,idx_comp);

 % Espera complete.
    H = anexar_etapa(H,datos_espera,2,length(datos_espera.t));

 % second excitation complete.
    H = anexar_etapa(H,datos2,2,length(datos2.t));

 % torque gravitacional incremental equivalente due to the masses.
    H.tau_masas = H.tau_masaX + H.tau_masaY;

end


function H = anexar_etapa(H,D,i1,i2)
%ANEXAR_ETAPA adds to segment of to stage to the historical.

    if i2 < i1
        return;
    end

    idx = i1:i2;

    campos_vector = { ...
        't', ...
        'wx_true','wy_true','wz_true', ...
        'wx','wy','wz', ...
        'gx_true','gy_true','gz_true', ...
        'gx','gy','gz', ...
        'dX','dY','dZ', ...
        'vX','vY','vZ' ...
    };

    for k = 1:numel(campos_vector)

        campo = campos_vector{k};

        H.(campo) = [ ...
            H.(campo); ...
            D.(campo)(idx,:) ...
        ];

    end

    H.tau_RW = [H.tau_RW;D.tau_RW(idx,:)];
    H.tau_masaX = [H.tau_masaX;D.tau_masaX(idx,:)];
    H.tau_masaY = [H.tau_masaY;D.tau_masaY(idx,:)];
    H.tau_exc = [H.tau_exc;D.tau_exc(idx,:)];
    H.Omega_RW = [H.Omega_RW;D.Omega_RW(idx,:)];
    H.h_RW = [H.h_RW;D.h_RW(idx,:)];

end


function graficar_historico_completo(H)
%PLOTR_HISTORICO_COMPLETO
% Final dashboard for the entire sequence:
%
% row 1: ideal and noisy angular velocity
% row 2: ideal and noisy gravity vector
% row 3: position of the three masses moving
% row 4: torques of the wheels and torque of excitation total
%
% Torques are shown in mNm and positions in mm.

    t = H.t;

    figure( ...
        'Name','Historico completo de la secuencia', ...
        'NumberTitle','off', ...
        'Units','normalized', ...
        'OuterPosition',[0 0 1 1] ...
    );

 % ============================================================
 % row 1 - angular VELOCITY
 % ============================================================

    ax(1) = subplot(4,3,1);
    plot(t,H.wx_true,'LineWidth',1.25);
    hold on;
    plot(t,H.wx,'LineWidth',0.45);
    grid on;
    ylabel('\omega_x [rad/s]');
    title('Velocidad angular X');
    legend('Ideal','Ruidosa','Location','best');

    ax(2) = subplot(4,3,2);
    plot(t,H.wy_true,'LineWidth',1.25);
    hold on;
    plot(t,H.wy,'LineWidth',0.45);
    grid on;
    ylabel('\omega_y [rad/s]');
    title('Velocidad angular Y');

    ax(3) = subplot(4,3,3);
    plot(t,H.wz_true,'LineWidth',1.25);
    hold on;
    plot(t,H.wz,'LineWidth',0.45);
    grid on;
    ylabel('\omega_z [rad/s]');
    title('Velocidad angular Z');


 % ============================================================
 % row 2 - GRAVITY VECTOR
 % ============================================================

    ax(4) = subplot(4,3,4);
    plot(t,H.gx_true,'LineWidth',1.25);
    hold on;
    plot(t,H.gx,'LineWidth',0.45);
    grid on;
    ylabel('g_x [m/s^2]');
    title('Gravedad X');
    legend('Ideal','Ruidosa','Location','best');

    ax(5) = subplot(4,3,5);
    plot(t,H.gy_true,'LineWidth',1.25);
    hold on;
    plot(t,H.gy,'LineWidth',0.45);
    grid on;
    ylabel('g_y [m/s^2]');
    title('Gravedad Y');

    ax(6) = subplot(4,3,6);
    plot(t,H.gz_true,'LineWidth',1.25);
    hold on;
    plot(t,H.gz,'LineWidth',0.45);
    grid on;
    ylabel('g_z [m/s^2]');
    title('Gravedad Z');


 % ============================================================
 % row 3 - position of masses
 % ============================================================

    ax(7) = subplot(4,3,7);
    plot(t,1e3*H.dX,'LineWidth',1.2);
    grid on;
    ylabel('dX [mm]');
    title('Masa movil X');

    ax(8) = subplot(4,3,8);
    plot(t,1e3*H.dY,'LineWidth',1.2);
    grid on;
    ylabel('dY [mm]');
    title('Masa movil Y');

    ax(9) = subplot(4,3,9);
    plot(t,1e3*H.dZ,'LineWidth',1.2);
    grid on;
    ylabel('dZ [mm]');
    title('Masa movil Z');


 % ============================================================
 % row 4 - torques
 %
 % RW: torque generated by the reaction wheels
 % masses : torque gravitacional incremental equivalente by
 % displacements X/and in the stages of excitation
 % Total: RW + incremental contribution of the masses.
 % Es to magnitud of diagnostic; in the second
 % estimation, the masses are modeled using Mr(t), DeltaJ(t),
 % and h_rel(t), not as to torque externo independiente.
 % ============================================================

    tau_masas = H.tau_masas;

    ax(10) = subplot(4,3,10);
    plot(t,1e3*H.tau_RW(:,1),'LineWidth',1.15);
    hold on;
    plot(t,1e3*tau_masas(:,1),'--','LineWidth',1.0);
    plot(t,1e3*H.tau_exc(:,1),':','LineWidth',1.25);
    grid on;
    ylabel('\tau_X [mNm]');
    xlabel('Tiempo [s]');
    title('Torques X');
    legend('RW','Masas','Excitacion total','Location','best');

    ax(11) = subplot(4,3,11);
    plot(t,1e3*H.tau_RW(:,2),'LineWidth',1.15);
    hold on;
    plot(t,1e3*tau_masas(:,2),'--','LineWidth',1.0);
    plot(t,1e3*H.tau_exc(:,2),':','LineWidth',1.25);
    grid on;
    ylabel('\tau_Y [mNm]');
    xlabel('Tiempo [s]');
    title('Torques Y');

    ax(12) = subplot(4,3,12);
    plot(t,1e3*H.tau_RW(:,3),'LineWidth',1.15);
    hold on;
    plot(t,1e3*tau_masas(:,3),'--','LineWidth',1.0);
    plot(t,1e3*H.tau_exc(:,3),':','LineWidth',1.25);
    grid on;
    ylabel('\tau_Z [mNm]');
    xlabel('Tiempo [s]');
    title('Torques Z');

    linkaxes(ax,'x');

    sgtitle('HISTORICO COMPLETO DE TODA LA SECUENCIA');

end

function graficar_velocidades_ruedas(H)
%PLOTR_VELOCIDADES_RUEDAS
% Plot the complete history of the relative rotor speeds of the
% three reaction wheels.
%
% H.Omega_RW almacena internamente in rad/s. for compare of way
% directly with the RW400-50 specification, it is converted to rpm.

    t = H.t;
    Omega_rpm = H.Omega_RW*(60/(2*pi));

    figure( ...
        'Name','Velocidades angulares de las ruedas de reaccion', ...
        'NumberTitle','off' ...
    );

    ax(1) = subplot(3,1,1);
    plot(t,Omega_rpm(:,1),'LineWidth',1.3);
    grid on;
    ylabel('\Omega_{RW,X} [rpm]');
    title('Rueda de reaccion X');

    ax(2) = subplot(3,1,2);
    plot(t,Omega_rpm(:,2),'LineWidth',1.3);
    grid on;
    ylabel('\Omega_{RW,Y} [rpm]');
    title('Rueda de reaccion Y');

    ax(3) = subplot(3,1,3);
    plot(t,Omega_rpm(:,3),'LineWidth',1.3);
    grid on;
    ylabel('\Omega_{RW,Z} [rpm]');
    xlabel('Tiempo [s]');
    title('Rueda de reaccion Z');

    linkaxes(ax,'x');

    sgtitle('VELOCIDADES ANGULARES DE LAS RUEDAS DE REACCION');

    pico_rpm = max(abs(Omega_rpm),[],1);

    fprintf('\nMaximum absolute reaction-wheel speed:\n');
    fprintf('  RW-X: %.3f rpm\n',pico_rpm(1));
    fprintf('  RW-and: %.3f rpm\n',pico_rpm(2));
    fprintf('  RW-Z: %.3f rpm\n',pico_rpm(3));

end
