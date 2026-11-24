clear;
clc;
close all;

%% ================================================================
% BASE + CUBESAT EXPERIMENT + FIGURE GENERATION FOR THE REPORT
%
% This script reuses the existing project functions:
% platform_parameters
% actuator_torque_parameters
% cubesat_parameters
% add_cubesat_to_platform
% run_identification_campaign
% extract_cubesat_parameters
%
% It automatically generates PNG figures (300 dpi) and CSV tables for
% Sections 7--13 of the LaTeX report.
% ================================================================

%% configuration of the experiment
nombre_cubesat = 'INNOCUBE';

usar_semilla_fija = true;
semilla_ruido = 20261001;

carpeta_salida = fullfile(pwd,'result_figures');
resolucion_png = 300;

if ~exist(carpeta_salida,'dir')
    mkdir(carpeta_salida);
end

if usar_semilla_fija
    rng(semilla_ruido,'twister');
else
    rng('shuffle');
end

fprintf('\n=============================================================\n');
fprintf(' EXPERIMENT + FIGURES FOR REPORT\n');
fprintf(' CubeSat: %s\n',nombre_cubesat);
fprintf(' Output folder: %s\n',carpeta_salida);
if usar_semilla_fija
    fprintf(' noise seed: %d\n',semilla_ruido);
else
    fprintf(' noise seed: random\n');
end
fprintf('=============================================================\n');

%% parameters
P_base = platform_parameters();
A0 = actuator_torque_parameters();
C = cubesat_parameters(nombre_cubesat);

%% EXECUTION OF THE TWO CAMPAIGNS
% Diagnostic figures created by run_identification_campaign are hidden so that
% the report script explicitly controls the exported figures.
visibilidad_original = get(groot,'DefaultFigureVisible');
set(groot,'DefaultFigureVisible','off');
cleanup_visibilidad = onCleanup(@() set(groot,'DefaultFigureVisible',visibilidad_original));

A_base = A0;
A_base.rCMz_obj = A0.rCMz_obj_base;
BASE = run_identification_campaign(P_base,A_base,'BASE');
close all;

P_con_cubesat = add_cubesat_to_platform(P_base,C);
A_cubesat = A0;
A_cubesat.rCMz_obj = A0.rCMz_obj_cubesat;
CON_CUBESAT = run_identification_campaign(P_con_cubesat,A_cubesat,['BASE_' C.id]);
close all;

RESULTADO_CUBESAT = extract_cubesat_parameters(BASE,CON_CUBESAT,C);
close all;

%% save result NUMERICO complete
archivo_resultados = fullfile( ...
    carpeta_salida, ...
    ['resultados_reporte_BASE_' C.id '.mat'] ...
);

save(archivo_resultados, ...
    'nombre_cubesat','semilla_ruido','usar_semilla_fija', ...
    'P_base','P_con_cubesat','A_base','A_cubesat','C', ...
    'BASE','CON_CUBESAT','RESULTADO_CUBESAT');

%% figures of each campaign
fprintf('\nGenerating figures for the BASE campaign...\n');
generar_figuras_campana(BASE,'BASE',carpeta_salida,resolucion_png);

fprintf('Generating figures for the BASE + CubeSat campaign...\n');
generar_figuras_campana( ...
    CON_CUBESAT,['BASE_' C.id],carpeta_salida,resolucion_png);

%% figures of comparison between campaigns
fprintf('Generando comparisons BASE / BASE + CubeSat...\n');
generar_figuras_comparacion_campanas( ...
    BASE,CON_CUBESAT,C,carpeta_salida,resolucion_png);

%% CUBESAT EXTRACTION AND METRIC FIGURES
fprintf('Generating CubeSat extraction and metric figures...\n');
METRICAS_REPORTE = generar_figuras_cubesat( ...
    RESULTADO_CUBESAT,C,carpeta_salida,resolucion_png);

%% TABLES CSV for INCORPORAR results in LATEX
generar_tablas_csv_campana(BASE,'BASE',carpeta_salida);
generar_tablas_csv_campana(CON_CUBESAT,['BASE_' C.id],carpeta_salida);
generar_tablas_csv_cubesat( ...
    RESULTADO_CUBESAT,C,METRICAS_REPORTE,carpeta_salida);

%% file of TEXTO with results CLAVE
archivo_resumen = fullfile( ...
    carpeta_salida,['resumen_resultados_' C.id '.txt']);

fid = fopen(archivo_resumen,'w');
if fid < 0
    warning('Could not create %s.',archivo_resumen);
else
    fprintf(fid,'RESULTS FOR REPORT - %s\n',C.nombre);
    fprintf(fid,'========================================\n\n');
    fprintf(fid,'noise seed: %d\n\n',semilla_ruido);

    escribir_resumen_campana(fid,BASE,'BASE');
    escribir_resumen_campana(fid,CON_CUBESAT,['BASE + ' C.nombre]);

    fprintf(fid,'\nCUBESAT EXTRACTION\n');
    fprintf(fid,'-----------------------\n');
    fprintf(fid,'error norma CM ideal [mm]   = %.9g\n', ...
        METRICAS_REPORTE.error_CM_norma_ideal_mm);
    fprintf(fid,'error norma CM measured [mm]  = %.9g\n', ...
        METRICAS_REPORTE.error_CM_norma_medida_mm);
    fprintf(fid,'error Frobenius J_CR ideal [%%]  = %.9g\n', ...
        METRICAS_REPORTE.error_JCR_F_ideal_pct);
    fprintf(fid,'error Frobenius J_CR measured [%%] = %.9g\n', ...
        METRICAS_REPORTE.error_JCR_F_medida_pct);

    fclose(fid);
end

save(fullfile(carpeta_salida,['report_metrics_' C.id '.mat']), ...
    'METRICAS_REPORTE');

set(groot,'DefaultFigureVisible',visibilidad_original);

fprintf('\n=============================================================\n');
fprintf(' GENERATION COMPLETED\n');
fprintf(' results MAT: %s\n',archivo_resultados);
fprintf(' Resumen TXT:    %s\n',archivo_resumen);
fprintf(' Figures/CSV:    %s\n',carpeta_salida);
fprintf('=============================================================\n');


%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function generar_figuras_campana(S,tag,carpeta,resolucion)

    H = S.historico;
    P = S.P;
    A = S.A;

    tag_archivo = sanitizar_nombre(tag);

 % Time limits of the stages actually used.
    t_bal_0 = H.t(1);
    t_bal_f = A.T_exc1;

    t_exc1_0 = A.T_exc1;
    t_exc1_f = A.T_fin_exc1;

    t_reb_0 = A.T_fin_exc1;
    t_reb_f = S.t_rebalanceo;

    t_comp_0 = S.t_rebalanceo;
    t_comp_f = S.t_compensacion;

    t_exc2_0 = S.t_inicio_etapa2;
    t_exc2_f = S.t_inicio_etapa2 + A.T_exc2_total;

 %% --------------------------------------------------------------------
 % section 7.2 - balancing X/and
 % ---------------------------------------------------------------------
    idx = intervalo(H.t,t_bal_0,t_bal_f);
    [phi_deg,theta_deg] = actitud_desde_gravedad( ...
        H.gx_true(idx),H.gy_true(idx),P.g);

    f = nueva_figura();
    tl = tiledlayout(f,2,1,'TileSpacing','compact','Padding','compact');

    ax = nexttile(tl);
    plot(ax,H.t(idx),phi_deg,'LineWidth',1.2); hold(ax,'on');
    plot(ax,H.t(idx),theta_deg,'LineWidth',1.2);
    yline(ax,+rad2deg(A.error_lim_xy),'--');
    yline(ax,-rad2deg(A.error_lim_xy),'--');
    xline(ax,A.T_balanceo1,':','Inicio control','LabelVerticalAlignment','bottom');
    grid(ax,'on');
    ylabel(ax,'Angulo [deg]');
    legend(ax,{'\phi','\theta','+\epsilon_{xy}','-\epsilon_{xy}'}, ...
        'Location','best');
    title(ax,['Balanceo X/Y - ' tag],'Interpreter','none');

    ax = nexttile(tl);
    plot(ax,H.t(idx),1e3*H.dX(idx),'LineWidth',1.2); hold(ax,'on');
    plot(ax,H.t(idx),1e3*H.dY(idx),'LineWidth',1.2);
    grid(ax,'on');
    xlabel(ax,'Tiempo [s]');
    ylabel(ax,'Desplazamiento [mm]');
    legend(ax,{'d_X','d_Y'},'Location','best');

    guardar_figura(f,carpeta, ...
        ['res_07_02_balanceo_' tag_archivo '.png'],resolucion);

 %% --------------------------------------------------------------------
 % section 7.3 - first excitation
 % ---------------------------------------------------------------------
    idx = intervalo(H.t,t_exc1_0,t_exc1_f);
    tref = S.R1.referencia.d(:);

    targetX = tref(1)*ones(nnz(idx),1);
    targetY = tref(2)*ones(nnz(idx),1);
    tloc = H.t(idx);
    targetX(tloc >= A.T_X1) = tref(1) + A.deltaX1;
    targetY(tloc >= A.T_Y1) = tref(2) + A.deltaY1;

    f = nueva_figura();
    tl = tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');

    ax = nexttile(tl);
    plot(ax,tloc,1e3*H.dX(idx),'LineWidth',1.2); hold(ax,'on');
    plot(ax,tloc,1e3*targetX,'--','LineWidth',1.0);
    plot(ax,tloc,1e3*H.dY(idx),'LineWidth',1.2);
    plot(ax,tloc,1e3*targetY,'--','LineWidth',1.0);
    grid(ax,'on');
    ylabel(ax,'Posicion [mm]');
    legend(ax,{'d_X','d_{X,target}','d_Y','d_{Y,target}'},'Location','best');
    title(ax,['Primera excitacion - ' tag],'Interpreter','none');

    ax = nexttile(tl);
    plot(ax,tloc,H.wx_true(idx),'LineWidth',1.2); hold(ax,'on');
    plot(ax,tloc,H.wy_true(idx),'LineWidth',1.2);
    plot(ax,tloc,H.wz_true(idx),'LineWidth',1.2);
    grid(ax,'on');
    ylabel(ax,'\omega [rad/s]');
    legend(ax,{'\omega_x','\omega_y','\omega_z'},'Location','best');

    ax = nexttile(tl);
    tau_masas = H.tau_masaX(idx,:) + H.tau_masaY(idx,:);
    plot(ax,tloc,1e3*tau_masas(:,1),'LineWidth',1.2); hold(ax,'on');
    plot(ax,tloc,1e3*tau_masas(:,2),'LineWidth',1.2);
    plot(ax,tloc,1e3*tau_masas(:,3),'LineWidth',1.2);
    grid(ax,'on');
    xlabel(ax,'Tiempo [s]');
    ylabel(ax,'\tau masas [mN m]');
    legend(ax,{'\tau_x','\tau_y','\tau_z'},'Location','best');

    guardar_figura(f,carpeta, ...
        ['res_07_03_primera_excitacion_' tag_archivo '.png'],resolucion);

 %% --------------------------------------------------------------------
 % section 7.4 - rebalancing
 % ---------------------------------------------------------------------
    idx = intervalo(H.t,t_reb_0,t_reb_f);
    [phi_deg,theta_deg] = actitud_desde_gravedad( ...
        H.gx_true(idx),H.gy_true(idx),P.g);

    f = nueva_figura();
    tl = tiledlayout(f,2,1,'TileSpacing','compact','Padding','compact');

    ax = nexttile(tl);
    plot(ax,H.t(idx),phi_deg,'LineWidth',1.2); hold(ax,'on');
    plot(ax,H.t(idx),theta_deg,'LineWidth',1.2);
    yline(ax,+rad2deg(A.error_lim_xy),'--');
    yline(ax,-rad2deg(A.error_lim_xy),'--');
    grid(ax,'on');
    ylabel(ax,'Angulo [deg]');
    legend(ax,{'\phi','\theta','+\epsilon_{xy}','-\epsilon_{xy}'}, ...
        'Location','best');
    title(ax,['Rebalanceo X/Y - ' tag],'Interpreter','none');

    ax = nexttile(tl);
    plot(ax,H.t(idx),1e3*H.dX(idx),'LineWidth',1.2); hold(ax,'on');
    plot(ax,H.t(idx),1e3*H.dY(idx),'LineWidth',1.2);
    grid(ax,'on');
    xlabel(ax,'Tiempo [s]');
    ylabel(ax,'Desplazamiento [mm]');
    legend(ax,{'d_X','d_Y'},'Location','best');

    guardar_figura(f,carpeta, ...
        ['res_07_04_rebalanceo_' tag_archivo '.png'],resolucion);

 %% --------------------------------------------------------------------
 % section 7.5 - compensation Z
 % ---------------------------------------------------------------------
    idx = intervalo(H.t,t_comp_0,t_comp_f);
    [phi_deg,theta_deg] = actitud_desde_gravedad( ...
        H.gx_true(idx),H.gy_true(idx),P.g);

    rCM_real = reconstruir_rCM(P,H.dX(idx),H.dY(idx),H.dZ(idx));
    rCMz_est = S.rCM_est_inicial(3) + (P.mZ/P.M_total).*H.dZ(idx);

    f = nueva_figura();
    tl = tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');

    ax = nexttile(tl);
    plot(ax,H.t(idx),1e3*H.dZ(idx),'LineWidth',1.2); hold(ax,'on');
    yline(ax,1e3*S.dZ_target,'--','d_{Z,target}');
    grid(ax,'on');
    ylabel(ax,'d_Z [mm]');
    title(ax,['Compensacion vertical - ' tag],'Interpreter','none');

    ax = nexttile(tl);
    plot(ax,H.t(idx),1e3*rCM_real(:,3),'LineWidth',1.2); hold(ax,'on');
    plot(ax,H.t(idx),1e3*rCMz_est,'--','LineWidth',1.2);
    yline(ax,1e3*S.rCMz_objetivo,':','Objetivo');
    grid(ax,'on');
    ylabel(ax,'r_{CM,z} [mm]');
    legend(ax,{'Real','Predicho por 10P','Objetivo'},'Location','best');

    ax = nexttile(tl);
    plot(ax,H.t(idx),phi_deg,'LineWidth',1.2); hold(ax,'on');
    plot(ax,H.t(idx),theta_deg,'LineWidth',1.2);
    grid(ax,'on');
    xlabel(ax,'Tiempo [s]');
    ylabel(ax,'Angulo [deg]');
    legend(ax,{'\phi','\theta'},'Location','best');

    guardar_figura(f,carpeta, ...
        ['res_07_05_compensacion_Z_' tag_archivo '.png'],resolucion);

 %% --------------------------------------------------------------------
 % sections 7.6 and 7.7 - PI and second excitation
 % ---------------------------------------------------------------------
    idx = intervalo(H.t,t_exc2_0,t_exc2_f);
    trel = H.t(idx)-t_exc2_0;
    dref = S.d_ref2(:);

    dX_target = perfil_bipolar( ...
        trel,dref(1),A.deltaX2, ...
        A.tX2_pos_ini,A.tX2_pos_fin,A.tX2_neg_ini,A.tX2_neg_fin);

    dY_target = perfil_bipolar( ...
        trel,dref(2),A.deltaY2, ...
        A.tY2_pos_ini,A.tY2_pos_fin,A.tY2_neg_ini,A.tY2_neg_fin);

    f = nueva_figura();
    tl = tiledlayout(f,2,1,'TileSpacing','compact','Padding','compact');

    ax = nexttile(tl);
    plot(ax,trel,1e3*H.dX(idx),'LineWidth',1.2); hold(ax,'on');
    plot(ax,trel,1e3*dX_target,'--','LineWidth',1.0);
    grid(ax,'on');
    ylabel(ax,'d_X [mm]');
    legend(ax,{'Real','Referencia'},'Location','best');
    title(ax,['Seguimiento PI durante segunda excitacion - ' tag], ...
        'Interpreter','none');

    ax = nexttile(tl);
    plot(ax,trel,1e3*H.dY(idx),'LineWidth',1.2); hold(ax,'on');
    plot(ax,trel,1e3*dY_target,'--','LineWidth',1.0);
    grid(ax,'on');
    xlabel(ax,'Tiempo relativo [s]');
    ylabel(ax,'d_Y [mm]');
    legend(ax,{'Real','Referencia'},'Location','best');

    guardar_figura(f,carpeta, ...
        ['res_07_06_PI_segunda_excitacion_' tag_archivo '.png'],resolucion);

    f = nueva_figura();
    tl = tiledlayout(f,4,1,'TileSpacing','compact','Padding','compact');

    ax = nexttile(tl);
    plot(ax,trel,H.wx_true(idx),'LineWidth',1.2); hold(ax,'on');
    plot(ax,trel,H.wy_true(idx),'LineWidth',1.2);
    plot(ax,trel,H.wz_true(idx),'LineWidth',1.2);
    grid(ax,'on');
    ylabel(ax,'\omega [rad/s]');
    legend(ax,{'\omega_x','\omega_y','\omega_z'},'Location','best');
    title(ax,['Segunda excitacion multieje - ' tag],'Interpreter','none');

    ax = nexttile(tl);
    plot(ax,trel,1e3*H.tau_masaY(idx,1),'LineWidth',1.2); hold(ax,'on');
    plot(ax,trel,1e3*H.tau_masaX(idx,2),'LineWidth',1.2);
    plot(ax,trel,1e3*H.tau_RW(idx,3),'LineWidth',1.2);
    grid(ax,'on');
    ylabel(ax,'Torque [mN m]');
    legend(ax,{'Masa Y \rightarrow X','Masa X \rightarrow Y','RW-Z'}, ...
        'Location','best');

    ax = nexttile(tl);
    plot(ax,trel,1e3*(H.dX(idx)-dref(1)),'LineWidth',1.2); hold(ax,'on');
    plot(ax,trel,1e3*(H.dY(idx)-dref(2)),'LineWidth',1.2);
    grid(ax,'on');
    ylabel(ax,'\Delta d [mm]');
    legend(ax,{'\Delta d_X','\Delta d_Y'},'Location','best');

    ax = nexttile(tl);
    omega_rpm = H.Omega_RW(idx,3)*60/(2*pi);
    plot(ax,trel,omega_rpm,'LineWidth',1.2);
    grid(ax,'on');
    xlabel(ax,'Tiempo relativo [s]');
    ylabel(ax,'\Omega_{RW,Z} [rpm]');

    guardar_figura(f,carpeta, ...
        ['res_07_07_perfil_segunda_excitacion_' tag_archivo '.png'],resolucion);

 %% --------------------------------------------------------------------
 % section 8 - noise of MEDICION during second excitation
 % ---------------------------------------------------------------------
    nwx = H.wx(idx)-H.wx_true(idx);
    nwy = H.wy(idx)-H.wy_true(idx);
    nwz = H.wz(idx)-H.wz_true(idx);

    f = nueva_figura();
    tl = tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');
    ruido = {nwx,nwy,nwz};
    nombres = {'n_{\omega,x}','n_{\omega,y}','n_{\omega,z}'};

    for k = 1:3
        ax = nexttile(tl);
        plot(ax,trel,ruido{k},'LineWidth',0.6); hold(ax,'on');
        yline(ax,+3*A.sigma_gyro,'--');
        yline(ax,-3*A.sigma_gyro,'--');
        grid(ax,'on');
        ylabel(ax,[nombres{k} ' [rad/s]']);
        if k == 1
            title(ax,['Ruido de giroscopo - ' tag],'Interpreter','none');
        end
    end
    xlabel(ax,'Tiempo relativo [s]');

    guardar_figura(f,carpeta, ...
        ['res_08_ruido_giroscopo_' tag_archivo '.png'],resolucion);

    ngx = H.gx(idx)-H.gx_true(idx);
    ngy = H.gy(idx)-H.gy_true(idx);
    ngz = H.gz(idx)-H.gz_true(idx);

    f = nueva_figura();
    tl = tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');
    ruido = {ngx,ngy,ngz};
    nombres = {'n_{g,x}','n_{g,y}','n_{g,z}'};

    for k = 1:3
        ax = nexttile(tl);
        plot(ax,trel,ruido{k},'LineWidth',0.6); hold(ax,'on');
        yline(ax,+3*A.sigma_g,'--');
        yline(ax,-3*A.sigma_g,'--');
        grid(ax,'on');
        ylabel(ax,[nombres{k} ' [m/s^2]']);
        if k == 1
            title(ax,['Ruido del vector gravedad - ' tag],'Interpreter','none');
        end
    end
    xlabel(ax,'Tiempo relativo [s]');

    guardar_figura(f,carpeta, ...
        ['res_08_ruido_gravedad_' tag_archivo '.png'],resolucion);

 %% --------------------------------------------------------------------
 % section 9 - results of the estimator 10P
 % ---------------------------------------------------------------------
    R1 = S.R1;

    etiquetasJ = {'J_{xx}','J_{yy}','J_{zz}','J_{xy}','J_{xz}','J_{yz}'};
    datosJ10 = 1e3*[ ...
        R1.theta_real(1:6), ...
        R1.ideal.theta_hat_ini(1:6), ...
        R1.ruido.theta_hat_ini(1:6) ...
    ];

    f = nueva_figura();
    ax = axes(f);
    bar(ax,datosJ10,'grouped');
    grid(ax,'on');
    xticks(ax,1:6);
    xticklabels(ax,etiquetasJ);
    ylabel(ax,'Inercia [10^{-3} kg m^2]');
    legend(ax,{'Real','Ideal','Ruidosa'},'Location','best');
    title(ax,['Estimacion 10P: componentes de J - ' tag], ...
        'Interpreter','none');
    guardar_figura(f,carpeta, ...
        ['res_09_10P_J_' tag_archivo '.png'],resolucion);

    datosCM = 1e3*[ ...
        R1.rCM_real(:), ...
        R1.ideal.rCM_hat_canonica(:), ...
        R1.ruido.rCM_hat_canonica(:) ...
    ];

    f = nueva_figura();
    ax = axes(f);
    bar(ax,datosCM,'grouped');
    grid(ax,'on');
    xticks(ax,1:3);
    xticklabels(ax,{'x','y','z'});
    ylabel(ax,'r_{CM} [mm]');
    legend(ax,{'Real','Ideal','Ruidosa'},'Location','best');
    title(ax,['Estimacion 10P: centro de masa canonico - ' tag], ...
        'Interpreter','none');
    guardar_figura(f,carpeta, ...
        ['res_09_10P_CM_' tag_archivo '.png'],resolucion);

    graficar_valores_singulares( ...
        R1.ideal.singular_values,R1.ruido.singular_values, ...
        ['Estimador 10P - ' tag], ...
        fullfile(carpeta,['res_09_10P_singular_values_' tag_archivo '.png']), ...
        resolucion);

 %% --------------------------------------------------------------------
 % section 10 - results of the estimator 7P
 % ---------------------------------------------------------------------
    R2 = S.R2;
    Jreal6 = matriz_a_6_local(R2.J_real_original);
    Jideal6 = matriz_a_6_local(R2.ideal.J_hat_canonica);
    Jmed6 = matriz_a_6_local(R2.medida.J_hat_canonica);

    f = nueva_figura();
    ax = axes(f);
    bar(ax,1e3*[Jreal6,Jideal6,Jmed6],'grouped');
    grid(ax,'on');
    xticks(ax,1:6);
    xticklabels(ax,etiquetasJ);
    ylabel(ax,'Inercia [10^{-3} kg m^2]');
    legend(ax,{'Real','Ideal','Medida'},'Location','best');
    title(ax,['Estimacion 7P: tensor canonico - ' tag], ...
        'Interpreter','none');
    guardar_figura(f,carpeta, ...
        ['res_10_7P_J_' tag_archivo '.png'],resolucion);

    graficar_valores_singulares( ...
        R2.ideal.singular_values,R2.medida.singular_values, ...
        ['Estimador 7P - ' tag], ...
        fullfile(carpeta,['res_10_7P_singular_values_' tag_archivo '.png']), ...
        resolucion);

 %% --------------------------------------------------------------------
 % section 11 - historical complete of to campaign
 % ---------------------------------------------------------------------
    f = nueva_figura();
    tl = tiledlayout(f,4,1,'TileSpacing','compact','Padding','compact');

    limites = [A.T_exc1,A.T_fin_exc1,S.t_rebalanceo, ...
        S.t_compensacion,S.t_inicio_etapa2];

    ax = nexttile(tl);
    plot(ax,H.t,H.wx_true,'LineWidth',1.0); hold(ax,'on');
    plot(ax,H.t,H.wy_true,'LineWidth',1.0);
    plot(ax,H.t,H.wz_true,'LineWidth',1.0);
    agregar_limites_etapas(ax,limites);
    grid(ax,'on');
    ylabel(ax,'\omega [rad/s]');
    legend(ax,{'\omega_x','\omega_y','\omega_z'},'Location','best');
    title(ax,['Historico completo - ' tag],'Interpreter','none');

    ax = nexttile(tl);
    plot(ax,H.t,H.gx_true,'LineWidth',1.0); hold(ax,'on');
    plot(ax,H.t,H.gy_true,'LineWidth',1.0);
    plot(ax,H.t,H.gz_true,'LineWidth',1.0);
    agregar_limites_etapas(ax,limites);
    grid(ax,'on');
    ylabel(ax,'g_B [m/s^2]');
    legend(ax,{'g_x','g_y','g_z'},'Location','best');

    ax = nexttile(tl);
    plot(ax,H.t,1e3*H.dX,'LineWidth',1.0); hold(ax,'on');
    plot(ax,H.t,1e3*H.dY,'LineWidth',1.0);
    plot(ax,H.t,1e3*H.dZ,'LineWidth',1.0);
    agregar_limites_etapas(ax,limites);
    grid(ax,'on');
    ylabel(ax,'d [mm]');
    legend(ax,{'d_X','d_Y','d_Z'},'Location','best');

    ax = nexttile(tl);
    plot(ax,H.t,1e3*H.tau_exc(:,1),'LineWidth',1.0); hold(ax,'on');
    plot(ax,H.t,1e3*H.tau_exc(:,2),'LineWidth',1.0);
    plot(ax,H.t,1e3*H.tau_exc(:,3),'LineWidth',1.0);
    agregar_limites_etapas(ax,limites);
    grid(ax,'on');
    xlabel(ax,'Tiempo [s]');
    ylabel(ax,'\tau_{exc} [mN m]');
    legend(ax,{'\tau_x','\tau_y','\tau_z'},'Location','best');

    guardar_figura(f,carpeta, ...
        ['res_11_historico_' tag_archivo '.png'],resolucion);
end


function generar_figuras_comparacion_campanas(BASE,BS,C,carpeta,resolucion)

 %% center of mass canonico final
    datos = 1e3*[ ...
        BASE.canonica_final.real.rCM(:), ...
        BASE.canonica_final.medida.rCM(:), ...
        BS.canonica_final.real.rCM(:), ...
        BS.canonica_final.medida.rCM(:) ...
    ];

    f = nueva_figura();
    ax = axes(f);
    bar(ax,datos,'grouped');
    grid(ax,'on');
    xticks(ax,1:3);
    xticklabels(ax,{'x','y','z'});
    ylabel(ax,'r_{CM} [mm]');
    legend(ax,{'BASE real','BASE medida', ...
        ['BASE+' C.id ' real'],['BASE+' C.id ' medida']}, ...
        'Location','best','Interpreter','none');
    title(ax,'Comparacion del centro de masa canonico entre campanas');
    guardar_figura(f,carpeta,'res_11_canonical_CM_comparison.png',resolucion);

 %% tensor canonico final
    J = 1e3*[ ...
        matriz_a_6_local(BASE.canonica_final.real.J), ...
        matriz_a_6_local(BASE.canonica_final.medida.J), ...
        matriz_a_6_local(BS.canonica_final.real.J), ...
        matriz_a_6_local(BS.canonica_final.medida.J) ...
    ];

    f = nueva_figura();
    ax = axes(f);
    bar(ax,J,'grouped');
    grid(ax,'on');
    xticks(ax,1:6);
    xticklabels(ax,{'J_{xx}','J_{yy}','J_{zz}','J_{xy}','J_{xz}','J_{yz}'});
    ylabel(ax,'Inercia [10^{-3} kg m^2]');
    legend(ax,{'BASE real','BASE medida', ...
        ['BASE+' C.id ' real'],['BASE+' C.id ' medida']}, ...
        'Location','best','Interpreter','none');
    title(ax,'Comparacion de tensores canonicos entre campanas');
    guardar_figura(f,carpeta,'res_11_canonical_J_comparison.png',resolucion);

 %% Estimated friction as a diagnostic quantity
    b = [ ...
        BASE.R1.canonica.ruido.b, BS.R1.canonica.ruido.b; ...
        BASE.R2.canonica.medida.b, BS.R2.canonica.medida.b ...
    ];

    f = nueva_figura();
    ax = axes(f);
    bar(ax,b,'grouped');
    grid(ax,'on');
    xticks(ax,1:2);
    xticklabels(ax,{'10P','7P'});
    ylabel(ax,'b [N m s/rad]');
    legend(ax,{'BASE',['BASE+' C.id]},'Location','best','Interpreter','none');
    title(ax,'Coeficiente de friccion estimado en ambas campanas');
    guardar_figura(f,carpeta, ...
        ['res_12_friccion_' sanitizar_nombre(C.id) '.png'],resolucion);
end


function M = generar_figuras_cubesat(R,C,carpeta,resolucion)

    id = sanitizar_nombre(C.id);

 %% Local center of mass in C axes
    cm_gt = 1e3*C.r_G_CM_C(:);
    cm_ideal = 1e3*R.ideal.r_G_CM_C(:);
    cm_med = 1e3*R.medida.r_G_CM_C(:);

    f = nueva_figura();
    ax = axes(f);
    bar(ax,[cm_gt,cm_ideal,cm_med],'grouped');
    grid(ax,'on');
    xticks(ax,1:3);
    xticklabels(ax,{'x_C','y_C','z_C'});
    ylabel(ax,'r_{G\rightarrow CM} [mm]');
    legend(ax,{'Ground truth','Ideal','Medida'},'Location','best');
    title(ax,['Centro de masa extraido - ' C.nombre],'Interpreter','none');
    guardar_figura(f,carpeta,['res_12_cubesat_CM_' id '.png'],resolucion);

 %% tensor with respect to to CR in axes C
 % uses as ground truth the tensor that actually entra to the plant.
 % this keeps the figure consistente with the simulation incluso if 
 % the interpretation of the input tensor is modified later.
    Jgt_C = C.R_BC'*C.J_CR_B*C.R_BC;
    Jideal_C = C.R_BC'*R.ideal.J_CR_B*C.R_BC;
    Jmed_C = C.R_BC'*R.medida.J_CR_B*C.R_BC;

    Jgt6 = matriz_a_6_local(Jgt_C);
    Jideal6 = matriz_a_6_local(Jideal_C);
    Jmed6 = matriz_a_6_local(Jmed_C);

    f = nueva_figura();
    ax = axes(f);
    bar(ax,1e3*[Jgt6,Jideal6,Jmed6],'grouped');
    grid(ax,'on');
    xticks(ax,1:6);
    xticklabels(ax,{'J_{xx}','J_{yy}','J_{zz}','J_{xy}','J_{xz}','J_{yz}'});
    ylabel(ax,'J^{CR} [10^{-3} kg m^2]');
    legend(ax,{'Ground truth','Ideal','Medida'},'Location','best');
    title(ax,['Tensor del CubeSat respecto a CR - ' C.nombre], ...
        'Interpreter','none');
    guardar_figura(f,carpeta,['res_12_cubesat_J_CR_' id '.png'],resolucion);

 %% Center-of-mass error by component
    err_cm_ideal = cm_ideal-cm_gt;
    err_cm_med = cm_med-cm_gt;

    f = nueva_figura();
    ax = axes(f);
    bar(ax,[err_cm_ideal,err_cm_med],'grouped');
    grid(ax,'on');
    xticks(ax,1:3);
    xticklabels(ax,{'x_C','y_C','z_C'});
    ylabel(ax,'Error [mm]');
    legend(ax,{'Ideal','Medida'},'Location','best');
    title(ax,['Error de centro de masa - ' C.nombre],'Interpreter','none');
    guardar_figura(f,carpeta,['res_13_error_CM_' id '.png'],resolucion);

 %% error of tensor with respect to to CR
    errJ_ideal = 1e3*(Jideal6-Jgt6);
    errJ_med = 1e3*(Jmed6-Jgt6);

    f = nueva_figura();
    ax = axes(f);
    bar(ax,[errJ_ideal,errJ_med],'grouped');
    grid(ax,'on');
    xticks(ax,1:6);
    xticklabels(ax,{'J_{xx}','J_{yy}','J_{zz}','J_{xy}','J_{xz}','J_{yz}'});
    ylabel(ax,'Error [10^{-3} kg m^2]');
    legend(ax,{'Ideal','Medida'},'Location','best');
    title(ax,['Error del tensor respecto a CR - ' C.nombre], ...
        'Interpreter','none');
    guardar_figura(f,carpeta,['res_13_error_J_CR_' id '.png'],resolucion);

 %% Global metrics
    M = struct;
    M.J_CR_GT_C = Jgt_C;
    M.J_CR_ideal_C = Jideal_C;
    M.J_CR_medida_C = Jmed_C;

    M.error_CM_ideal_mm = err_cm_ideal;
    M.error_CM_medida_mm = err_cm_med;
    M.error_CM_norma_ideal_mm = norm(err_cm_ideal);
    M.error_CM_norma_medida_mm = norm(err_cm_med);

    M.error_JCR_ideal = Jideal6-Jgt6;
    M.error_JCR_medida = Jmed6-Jgt6;
    M.error_JCR_F_ideal_pct = ...
        100*norm(Jideal_C-Jgt_C,'fro')/norm(Jgt_C,'fro');
    M.error_JCR_F_medida_pct = ...
        100*norm(Jmed_C-Jgt_C,'fro')/norm(Jgt_C,'fro');

    datos_metricas = [ ...
        M.error_CM_norma_ideal_mm, M.error_CM_norma_medida_mm; ...
        M.error_JCR_F_ideal_pct, M.error_JCR_F_medida_pct ...
    ];

    f = nueva_figura();
    tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

    ax = nexttile(tl);
    bar(ax,datos_metricas(1,:));
    grid(ax,'on');
    xticks(ax,1:2);
    xticklabels(ax,{'Ideal','Medida'});
    ylabel(ax,'Norma error CM [mm]');

    ax = nexttile(tl);
    bar(ax,datos_metricas(2,:));
    grid(ax,'on');
    xticks(ax,1:2);
    xticklabels(ax,{'Ideal','Medida'});
    ylabel(ax,'Error Frobenius J^{CR} [%]');

    sgtitle(tl,['Metricas globales - ' C.nombre],'Interpreter','none');
    guardar_figura(f,carpeta,['res_13_metricas_globales_' id '.png'],resolucion);
end


function generar_tablas_csv_campana(S,tag,carpeta)

    tag = sanitizar_nombre(tag);

 %% table 10P
    nombres10 = {'Jxx';'Jyy';'Jzz';'Jxy';'Jxz';'Jyz';'Sx';'Sy';'Sz';'b'};
    T10 = table( ...
        nombres10, ...
        S.R1.theta_real(:), ...
        S.R1.ideal.theta_hat_ini(:), ...
        S.R1.ruido.theta_hat_ini(:), ...
        'VariableNames',{'Parametro','Real','Ideal','Ruidosa'});
    writetable(T10,fullfile(carpeta,['table_09_10P_' tag '.csv']));

 %% table 7P in configuration canonical
    Jreal6 = matriz_a_6_local(S.R2.J_real_original);
    Jideal6 = matriz_a_6_local(S.R2.ideal.J_hat_canonica);
    Jmed6 = matriz_a_6_local(S.R2.medida.J_hat_canonica);

    nombres7 = {'Jxx';'Jyy';'Jzz';'Jxy';'Jxz';'Jyz';'b'};
    real7 = [Jreal6;S.P.b];
    ideal7 = [Jideal6;S.R2.ideal.b_hat];
    med7 = [Jmed6;S.R2.medida.b_hat];

    T7 = table(nombres7,real7,ideal7,med7, ...
        'VariableNames',{'Parametro','Real','Ideal','Medida'});
    writetable(T7,fullfile(carpeta,['table_10_7P_' tag '.csv']));

 %% Regressor diagnostics
    estimador = {'10P ideal';'10P ruido';'7P ideal';'7P medida'};
    rango = [S.R1.ideal.rank_Phi;S.R1.ruido.rank_Phi; ...
        S.R2.ideal.rank_Phi;S.R2.medida.rank_Phi];
    condicion = [S.R1.ideal.cond_Phi;S.R1.ruido.cond_Phi; ...
        S.R2.ideal.cond_Phi;S.R2.medida.cond_Phi];
    RMSE = [S.R1.ideal.RMSE_batch;S.R1.ruido.RMSE_batch; ...
        S.R2.ideal.RMSE_batch;S.R2.medida.RMSE_batch];

    TD = table(estimador,rango,condicion,RMSE);
    writetable(TD,fullfile(carpeta,['table_diagnostics_' tag '.csv']));
end


function generar_tablas_csv_cubesat(R,C,M,carpeta)

    id = sanitizar_nombre(C.id);

 %% CM
    eje = {'x';'y';'z'};
    gt_mm = 1e3*C.r_G_CM_C(:);
    ideal_mm = 1e3*R.ideal.r_G_CM_C(:);
    medida_mm = 1e3*R.medida.r_G_CM_C(:);
    error_ideal_mm = ideal_mm-gt_mm;
    error_medida_mm = medida_mm-gt_mm;

    TCM = table(eje,gt_mm,ideal_mm,medida_mm,error_ideal_mm,error_medida_mm, ...
        'VariableNames',{'Eje','GroundTruthmm','Idealmm','Medidamm', ...
        'ErrorIdealmm','ErrorMedidamm'});
    writetable(TCM,fullfile(carpeta,['table_12_CM_' id '.csv']));

 %% tensor CR
    parametro = {'Jxx';'Jyy';'Jzz';'Jxy';'Jxz';'Jyz'};
    gt = matriz_a_6_local(M.J_CR_GT_C);
    ideal = matriz_a_6_local(M.J_CR_ideal_C);
    medida = matriz_a_6_local(M.J_CR_medida_C);
    error_ideal = ideal-gt;
    error_medida = medida-gt;

    TJ = table(parametro,gt,ideal,medida,error_ideal,error_medida, ...
        'VariableNames',{'Parametro','GroundTruthkgm2','Idealkgm2', ...
        'Medidakgm2','ErrorIdealkgm2','ErrorMedidakgm2'});
    writetable(TJ,fullfile(carpeta,['table_12_J_CR_' id '.csv']));

 %% Global metrics
    metrica = {'Norma error CM';'Error Frobenius JCR'};
    ideal_m = [M.error_CM_norma_ideal_mm;M.error_JCR_F_ideal_pct];
    medida_m = [M.error_CM_norma_medida_mm;M.error_JCR_F_medida_pct];
    TM = table(metrica,ideal_m,medida_m, ...
        'VariableNames',{'Metrica','Ideal','Medida'});
    writetable(TM,fullfile(carpeta,['table_13_metrics_' id '.csv']));
end


function escribir_resumen_campana(fid,S,nombre)

    fprintf(fid,'\n%s\n',nombre);
    fprintf(fid,'%s\n',repmat('-',1,length(nombre)));
    fprintf(fid,'rango 10P ideal = %d / 10\n',S.R1.ideal.rank_Phi);
    fprintf(fid,'cond 10P ideal  = %.9g\n',S.R1.ideal.cond_Phi);
    fprintf(fid,'RMSE 10P ideal  = %.9g\n',S.R1.ideal.RMSE_batch);
    fprintf(fid,'rango 10P noise = %d / 10\n',S.R1.ruido.rank_Phi);
    fprintf(fid,'cond 10P noise  = %.9g\n',S.R1.ruido.cond_Phi);
    fprintf(fid,'RMSE 10P noise  = %.9g\n',S.R1.ruido.RMSE_batch);
    fprintf(fid,'rango 7P ideal  = %d / 7\n',S.R2.ideal.rank_Phi);
    fprintf(fid,'cond 7P ideal   = %.9g\n',S.R2.ideal.cond_Phi);
    fprintf(fid,'RMSE 7P ideal   = %.9g\n',S.R2.ideal.RMSE_batch);
    fprintf(fid,'rango 7P measured = %d / 7\n',S.R2.medida.rank_Phi);
    fprintf(fid,'cond 7P measured  = %.9g\n',S.R2.medida.cond_Phi);
    fprintf(fid,'RMSE 7P measured  = %.9g\n',S.R2.medida.RMSE_batch);
    fprintf(fid,'rCM canonico true [mm] = %.9g %.9g %.9g\n', ...
        1e3*S.canonica_final.real.rCM);
    fprintf(fid,'rCM canonico med  [mm] = %.9g %.9g %.9g\n', ...
        1e3*S.canonica_final.medida.rCM);
    fprintf(fid,'dZ target [mm] = %.9g\n',1e3*S.dZ_target);
    fprintf(fid,'target rCMz [mm] = %.9g\n',1e3*S.rCMz_objetivo);
end


function graficar_valores_singulares(sideal,smed,titulo,ruta,resolucion)

    sideal = sideal(:);
    smed = smed(:);

    if isempty(sideal) || isempty(smed)
        return;
    end

    sideal_n = sideal/max(sideal(1),eps);
    smed_n = smed/max(smed(1),eps);

    f = nueva_figura();
    ax = axes(f);
    semilogy(ax,1:numel(sideal_n),max(sideal_n,eps),'o-','LineWidth',1.2); hold(ax,'on');
    semilogy(ax,1:numel(smed_n),max(smed_n,eps),'s-','LineWidth',1.2);
    grid(ax,'on');
    xlabel(ax,'Indice del valor singular');
    ylabel(ax,'\sigma_i/\sigma_{max}');
    legend(ax,{'Ideal','Con ruido'},'Location','best');
    title(ax,['Valores singulares normalizados - ' titulo],'Interpreter','none');
    guardar_figura_ruta(f,ruta,resolucion);
end


function [phi_deg,theta_deg] = actitud_desde_gravedad(gx,gy,g)

    s_theta = gx./g;
    s_theta = max(-1,min(1,s_theta));
    theta = asin(s_theta);

    c_theta = cos(theta);
    c_theta(abs(c_theta) < 1e-9) = 1e-9;

    s_phi = -gy./(g.*c_theta);
    s_phi = max(-1,min(1,s_phi));
    phi = asin(s_phi);

    phi_deg = rad2deg(phi);
    theta_deg = rad2deg(theta);
end


function rCM = reconstruir_rCM(P,dX,dY,dZ)

    dX = dX(:);
    dY = dY(:);
    dZ = dZ(:);
    N = numel(dX);
    rCM = zeros(N,3);

    for k = 1:N
        rX = P.rX0 + P.eX*dX(k);
        rY = P.rY0 + P.eY*dY(k);
        rZ = P.rZ0 + P.eZ*dZ(k);

        S = P.M_p*P.r_P + P.mX*rX + P.mY*rY + P.mZ*rZ;
        rCM(k,:) = (S/P.M_total)';
    end
end


function dtarget = perfil_bipolar(trel,dref,delta,tpi,tpf,tni,tnf)

    trel = trel(:);
    dtarget = dref*ones(size(trel));
    dtarget(trel >= tpi & trel < tpf) = dref + delta;
    dtarget(trel >= tni & trel < tnf) = dref - delta;
end


function idx = intervalo(t,t0,tf)
    idx = (t >= t0-1e-9) & (t <= tf+1e-9);
end


function v = matriz_a_6_local(J)
    v = [J(1,1);J(2,2);J(3,3);J(1,2);J(1,3);J(2,3)];
end


function agregar_limites_etapas(ax,limites)
    for k = 1:numel(limites)
        xline(ax,limites(k),':');
    end
end


function f = nueva_figura()
    f = figure( ...
        'Visible','off', ...
        'Color','w', ...
        'Units','pixels', ...
        'Position',[100 100 1100 760]);
end


function guardar_figura(f,carpeta,nombre,resolucion)
    guardar_figura_ruta(f,fullfile(carpeta,nombre),resolucion);
end


function guardar_figura_ruta(f,ruta,resolucion)
    drawnow;
    if exist('exportgraphics','file') == 2
        exportgraphics(f,ruta,'Resolution',resolucion);
    else
        print(f,ruta,'-dpng',['-r' num2str(resolucion)]);
    end
    close(f);
end


function s = sanitizar_nombre(s)
    s = char(string(s));
    s = regexprep(s,'[^A-Za-z0-9_-]','_');
end
