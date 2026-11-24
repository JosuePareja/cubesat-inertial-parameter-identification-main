clear;
clc;
close all;

%% ================================================================
% GENERATE_FINAL_ARTICLE_RESULTS
%
% Master script for the article.
%
% Runs one BASE campaign and five BASE+CubeSat campaigns:
% InnoCube, TPA-1, SOURCE, M3Sat, and InflateSail.
%
% Requires the .m files included in this package.
%
% Outputs in ./article_results:
% fig_results_01_campaigns.png
% fig_results_02_vertical_compensation.png
% fig_results_03_second_excitation.png
% fig_results_04_innocube_final.png
% table_innocube_final.csv / .tex
% table_friction.csv / .tex
% table_multicubesat.csv / .tex
% results_macros.tex
% complete_results.mat
%
% Article convention:
% Reference | ideal | Noisy
% ================================================================

%% configuration
carpeta_script = fileparts(mfilename('fullpath'));
if isempty(carpeta_script)
    carpeta_script = pwd;
end
addpath(carpeta_script);

carpeta_salida = fullfile(carpeta_script,'article_results');
if ~exist(carpeta_salida,'dir')
    mkdir(carpeta_salida);
end

semilla_base = 20261001;
resolucion_png = 300;

% Article figure style
set(groot,'DefaultAxesFontSize',14);
set(groot,'DefaultAxesLabelFontSizeMultiplier',16/14);
set(groot,'DefaultLegendFontSize',16);
set(groot,'DefaultTextFontSize',14);

vis_original = get(groot,'DefaultFigureVisible');
set(groot,'DefaultFigureVisible','off');
cleanup_vis = onCleanup(@() set(groot,'DefaultFigureVisible',vis_original));

fprintf('\n=============================================================\n');
fprintf(' final ARTICLE RESULTS\n');
fprintf(' PI balancing: Kp = 0.08, Ki = 0.002\n');
fprintf(' CubeSats: InnoCube, TPA-1, SOURCE, M3Sat, InflateSail\n');
fprintf(' Output: %s\n',carpeta_salida);
fprintf('=============================================================\n');

%% parameters BASE
P_base = platform_parameters();
A0 = actuator_torque_parameters();

if ~isfield(A0,'Kp_bal') || ~isfield(A0,'Ki_bal')
    error(['This script requires the PI version of the package. ', ...
           'Fields A.Kp_bal and A.Ki_bal were not found.']);
end

%% campaign BASE
rng(semilla_base,'twister');
A_base = A0;
A_base.rCMz_obj = A0.rCMz_obj_base;

fprintf('\n--- campaign BASE ---\n');
BASE = run_identification_campaign(P_base,A_base,'BASE');
close all;

%% campaigns CUBESAT
ids = {'INNOCUBE','TPA1','SOURCE','M3SAT','INFLATESAIL'};
RESULTADOS = struct;
COMBINADAS = struct;
CUBESATS = struct;

for k = 1:numel(ids)

    id = ids{k};
    C = cubesat_parameters(id);
    P_comb = add_cubesat_to_platform(P_base,C);

    A_comb = A0;
    A_comb.rCMz_obj = A0.rCMz_obj_cubesat;

 % Fixed and reproducible seed, different for each spacecraft.
    rng(semilla_base+k,'twister');

    fprintf('\n--- campaign BASE + %s ---\n',C.nombre);
    S = run_identification_campaign(P_comb,A_comb,['BASE_' C.id]);
    close all;

    R = extract_cubesat_parameters(BASE,S,C);
    close all;

    CUBESATS.(C.id) = C;
    COMBINADAS.(C.id) = S;
    RESULTADOS.(C.id) = R;
end

%% figures principal: INNOCUBE
INNO = RESULTADOS.INNOCUBE;
S_INNO = COMBINADAS.INNOCUBE;
C_INNO = CUBESATS.INNOCUBE;

generar_figura_campanas(BASE,S_INNO,carpeta_salida,resolucion_png);
generar_figura_compensacion(BASE,S_INNO,carpeta_salida,resolucion_png);
generar_figura_exc2(BASE,S_INNO,carpeta_salida,resolucion_png);
generar_figura_innocube_final(INNO,C_INNO,carpeta_salida,resolucion_png);

%% TABLES
T_inno = crear_tabla_innocube(INNO,C_INNO);
writetable(T_inno,fullfile(carpeta_salida,'table_innocube_final.csv'));
escribir_tabla_innocube_tex(T_inno, ...
    fullfile(carpeta_salida,'table_innocube_final.tex'));

T_b = crear_table_friction(BASE,S_INNO);
writetable(T_b,fullfile(carpeta_salida,'table_friction.csv'));
escribir_table_friction_tex(T_b, ...
    fullfile(carpeta_salida,'table_friction.tex'));

T_multi = crear_table_multicubesat(RESULTADOS,CUBESATS);
writetable(T_multi,fullfile(carpeta_salida,'table_multicubesat.csv'));
escribir_table_multicubesat_tex(T_multi, ...
    fullfile(carpeta_salida,'table_multicubesat.tex'));

%% MACROS for ABSTRACT, RESULTS and CONCLUSIONS
crear_macros_articulo(T_multi,INNO,C_INNO, ...
    fullfile(carpeta_salida,'results_macros.tex'));

%% save TODO
save(fullfile(carpeta_salida,'complete_results.mat'), ...
    'BASE','RESULTADOS','COMBINADAS','CUBESATS', ...
    'T_inno','T_b','T_multi','A0','semilla_base');

set(groot,'DefaultFigureVisible',vis_original);

fprintf('\n=============================================================\n');
fprintf(' COMPLETED\n');
fprintf(' Figures: 4 PNG to 300 dpi\n');
fprintf(' Tables: 3 CSV + 3 LaTeX\n');
fprintf(' Macros: results_macros.tex\n');
fprintf(' MAT: complete_results.mat\n');
fprintf(' Folder: %s\n',carpeta_salida);
fprintf('=============================================================\n');


%% ========================================================================
% figure 1 - campaigns COMPLETAS
% ========================================================================

function generar_figura_campanas(BASE,COMB,carpeta,resolucion)

    f = figure('Units','centimeters','Position',[2 2 28 21], ...
        'PaperPositionMode','auto');
    tl = tiledlayout(f,3,2,'TileSpacing','compact','Padding','compact');

    dibujar_campana_columna(tl,BASE,1);
    dibujar_campana_columna(tl,COMB,2);

    exportgraphics(f,fullfile(carpeta,'fig_results_01_campaigns.png'), ...
        'Resolution',resolucion);
    close(f);
end

function dibujar_campana_columna(tl,S,col)

    H = S.historico;

    ax = nexttile(tl,col);
    plot(ax,H.t,H.wx_true,'LineWidth',1.0); hold(ax,'on');
    plot(ax,H.t,H.wy_true,'LineWidth',1.0);
    plot(ax,H.t,H.wz_true,'LineWidth',1.0);
    grid(ax,'on');
    ylabel(ax,'Angular rate [rad/s]');
    legend(ax,{'\omega_x','\omega_y','\omega_z'},'Location','best');

    ax = nexttile(tl,2+col);
    plot(ax,H.t,1e3*H.dX,'LineWidth',1.0); hold(ax,'on');
    plot(ax,H.t,1e3*H.dY,'LineWidth',1.0);
    plot(ax,H.t,1e3*H.dZ,'LineWidth',1.0);
    grid(ax,'on');
    ylabel(ax,'Mass position [mm]');
    legend(ax,{'d_X','d_Y','d_Z'},'Location','best');

    ax = nexttile(tl,4+col);
    tau_total = H.tau_RW + H.tau_masaX + H.tau_masaY;
    plot(ax,H.t,1e3*tau_total(:,1),'LineWidth',1.0); hold(ax,'on');
    plot(ax,H.t,1e3*tau_total(:,2),'LineWidth',1.0);
    plot(ax,H.t,1e3*tau_total(:,3),'LineWidth',1.0);
    grid(ax,'on');
    xlabel(ax,'Time [s]');
    ylabel(ax,'Known excitation [mN m]');
    legend(ax,{'\tau_x','\tau_y','\tau_z'},'Location','best');
end


%% ========================================================================
% figure 2 - compensation VERTICAL
% ========================================================================

function generar_figura_compensacion(BASE,COMB,carpeta,resolucion)

    f = figure('Units','centimeters','Position',[2 2 28 16], ...
        'PaperPositionMode','auto');
    tl = tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');

    dibujar_comp_columna(tl,BASE,1);
    dibujar_comp_columna(tl,COMB,2);

    exportgraphics(f,fullfile(carpeta, ...
        'fig_results_02_vertical_compensation.png'), ...
        'Resolution',resolucion);
    close(f);
end

function dibujar_comp_columna(tl,S,col)

    H = S.historico;
    idx = H.t >= S.t_rebalanceo & H.t <= S.t_compensacion;

    trel = H.t(idx)-S.t_rebalanceo;

    ax = nexttile(tl,col);
    plot(ax,trel,1e3*H.dZ(idx),'LineWidth',1.2);
    grid(ax,'on');
    ylabel(ax,'d_Z [mm]');

    [phi_deg,theta_deg] = actitud_desde_gravedad( ...
        H.gx_true(idx),H.gy_true(idx),S.P.g);

    ax = nexttile(tl,2+col);
    plot(ax,trel,phi_deg,'LineWidth',1.1); hold(ax,'on');
    plot(ax,trel,theta_deg,'LineWidth',1.1);
    grid(ax,'on');
    xlabel(ax,'Time from compensation start [s]');
    ylabel(ax,'Angle [deg]');
    legend(ax,{'\phi','\theta'},'Location','best');
end


%% ========================================================================
% figure 3 - second excitation
% ========================================================================

function generar_figura_exc2(BASE,COMB,carpeta,resolucion)

    f = figure('Units','centimeters','Position',[2 2 28 16], ...
        'PaperPositionMode','auto');
    tl = tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');

    dibujar_exc2_columna(tl,BASE,1);
    dibujar_exc2_columna(tl,COMB,2);

    exportgraphics(f,fullfile(carpeta, ...
        'fig_results_03_second_excitation.png'), ...
        'Resolution',resolucion);
    close(f);
end

function dibujar_exc2_columna(tl,S,col)

    H = S.historico;
    t0 = S.t_inicio_etapa2;
    tf = t0 + S.A.T_exc2_total;
    idx = H.t >= t0 & H.t <= tf;
    trel = H.t(idx)-t0;

    dX0 = H.dX(find(idx,1,'first'));
    dY0 = H.dY(find(idx,1,'first'));

    ax = nexttile(tl,col);
    plot(ax,trel,1e3*(H.dX(idx)-dX0),'LineWidth',1.1); hold(ax,'on');
    plot(ax,trel,1e3*(H.dY(idx)-dY0),'LineWidth',1.1);
    grid(ax,'on');
    ylabel(ax,'Relative position [mm]');
    legend(ax,{'d_X-d_{X,ref}','d_Y-d_{Y,ref}'},'Location','best');

    ax = nexttile(tl,2+col);
    plot(ax,trel,H.wx_true(idx),'LineWidth',1.1); hold(ax,'on');
    plot(ax,trel,H.wy_true(idx),'LineWidth',1.1);
    plot(ax,trel,H.wz_true(idx),'LineWidth',1.1);
    grid(ax,'on');
    xlabel(ax,'Time from second excitation start [s]');
    ylabel(ax,'Angular rate [rad/s]');
    legend(ax,{'\omega_x','\omega_y','\omega_z'},'Location','best');
end


%% ========================================================================
% FIGURE 4 - FINAL INNOCUBE RESULT
% ========================================================================

function generar_figura_innocube_final(R,C,carpeta,resolucion)

    cm_ref = 1e3*C.r_G_CM_C(:);
    cm_ideal = 1e3*R.ideal.r_G_CM_C(:);
    cm_noisy = 1e3*R.medida.r_G_CM_C(:);

    J_ref = matriz_a_6(C.J_reportado_C);
    J_ideal = matriz_a_6(R.ideal.J_CM_C);
    J_noisy = matriz_a_6(R.medida.J_CM_C);

    f = figure('Units','centimeters','Position',[2 2 28 11], ...
        'PaperPositionMode','auto');
    tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

    ax = nexttile(tl);
    bar(ax,[cm_ref,cm_ideal,cm_noisy]);
    grid(ax,'on');
    xticks(ax,1:3);
    xticklabels(ax,{'x_C','y_C','z_C'});
    ylabel(ax,'Center-of-mass offset [mm]');
    legend(ax,{'Reference','Ideal','Noisy'},'Location','best');

    ax = nexttile(tl);
    bar(ax,[J_ref,J_ideal,J_noisy]);
    grid(ax,'on');
    xticks(ax,1:6);
    xticklabels(ax,{'J_{xx}','J_{yy}','J_{zz}','J_{xy}','J_{xz}','J_{yz}'});
    ylabel(ax,'Inertia [kg m^2]');
    legend(ax,{'Reference','Ideal','Noisy'},'Location','best');

    exportgraphics(f,fullfile(carpeta, ...
        'fig_results_04_innocube_final.png'), ...
        'Resolution',resolucion);
    close(f);
end


%% ========================================================================
% TABLES
% ========================================================================

function T = crear_tabla_innocube(R,C)

    nombres = { ...
        'r_GCM_x';'r_GCM_y';'r_GCM_z'; ...
        'Jxx';'Jyy';'Jzz';'Jxy';'Jxz';'Jyz'};
    unidades = { ...
        'mm';'mm';'mm'; ...
        'kg m^2';'kg m^2';'kg m^2';'kg m^2';'kg m^2';'kg m^2'};

    ref = [1e3*C.r_G_CM_C(:); matriz_a_6(C.J_reportado_C)];
    ideal = [1e3*R.ideal.r_G_CM_C(:); matriz_a_6(R.ideal.J_CM_C)];
    noisy = [1e3*R.medida.r_G_CM_C(:); matriz_a_6(R.medida.J_CM_C)];

    abs_err = abs(noisy-ref);
    rel_err = 100*abs_err./abs(ref);

    T = table(nombres,unidades,ref,ideal,noisy,abs_err,rel_err, ...
        'VariableNames',{'Parameter','Unit','Reference','Ideal','Noisy', ...
        'AbsoluteError','RelativeError_pct'});
end

function T = crear_table_friction(BASE,COMB)

    P = BASE.P;

    Stage = {'BASE 10P';'BASE 7P';'BASE+InnoCube 10P';'BASE+InnoCube 7P'};
    Reference = repmat(P.b,4,1);

    Ideal = [ ...
        BASE.R1.canonica.ideal.b; ...
        BASE.R2.canonica.ideal.b; ...
        COMB.R1.canonica.ideal.b; ...
        COMB.R2.canonica.ideal.b];

    Noisy = [ ...
        BASE.R1.canonica.ruido.b; ...
        BASE.R2.canonica.medida.b; ...
        COMB.R1.canonica.ruido.b; ...
        COMB.R2.canonica.medida.b];

    T = table(Stage,Reference,Ideal,Noisy);
end

function T = crear_table_multicubesat(Rs,Cs)

    ids = {'INNOCUBE','TPA1','SOURCE','M3SAT','INFLATESAIL'};
    Spacecraft = cell(numel(ids),1);
    CM_Ideal_mm = zeros(numel(ids),1);
    CM_Noisy_mm = zeros(numel(ids),1);
    Inertia_Ideal_pct = zeros(numel(ids),1);
    Inertia_Noisy_pct = zeros(numel(ids),1);

    for k = 1:numel(ids)

        id = ids{k};
        R = Rs.(id);
        C = Cs.(id);

        Spacecraft{k} = C.nombre;
        CM_Ideal_mm(k) = norm(1e3*(R.ideal.r_G_CM_C-C.r_G_CM_C));
        CM_Noisy_mm(k) = norm(1e3*(R.medida.r_G_CM_C-C.r_G_CM_C));

        if C.productos_inercia_reportados
            Jref = C.J_reportado_C;
            Jid = R.ideal.J_CM_C;
            Jno = R.medida.J_CM_C;

            Inertia_Ideal_pct(k) = 100*norm(Jid-Jref,'fro')/norm(Jref,'fro');
            Inertia_Noisy_pct(k) = 100*norm(Jno-Jref,'fro')/norm(Jref,'fro');
        else
            r = diag(C.J_reportado_C);
            i = diag(R.ideal.J_CM_C);
            n = diag(R.medida.J_CM_C);

            Inertia_Ideal_pct(k) = 100*norm(i-r)/norm(r);
            Inertia_Noisy_pct(k) = 100*norm(n-r)/norm(r);
        end
    end

    T = table(Spacecraft,CM_Ideal_mm,CM_Noisy_mm, ...
        Inertia_Ideal_pct,Inertia_Noisy_pct);
end


%% ========================================================================
% EXPORTACION LATEX
% ========================================================================

function escribir_tabla_innocube_tex(T,archivo)

    fid = abrir_archivo(archivo);
    c = onCleanup(@() fclose(fid));

    fprintf(fid,'\\begin{table*}[t]\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,['\\caption{Reference and estimated InnoCube mass properties. ', ...
        'The Ideal column is obtained from noise-free simulated measurements, ', ...
        'whereas Noisy uses the sensor-noise levels defined for the numerical ', ...
        'campaign. The absolute and relative errors correspond to the Noisy ', ...
        'estimate.}\n']);
    fprintf(fid,'\\label{tab:innocube_final_results}\n');
    fprintf(fid,'\\resizebox{\\textwidth}{!}{%%\n');
    fprintf(fid,'\\begin{tabular}{llrrrrr}\n\\toprule\n');
    fprintf(fid,['Parameter & Unit & Reference & ideal & Noisy & ', ...
        'Absolute error & Relative error [\\%%] \\\\\\\\\n\\midrule\n']);

    latex_names = { ...
        '$r_{G\\rightarrow CM,x}^{C}$'; ...
        '$r_{G\\rightarrow CM,y}^{C}$'; ...
        '$r_{G\\rightarrow CM,z}^{C}$'; ...
        '$J_{xx}^{CM,C}$';'$J_{yy}^{CM,C}$';'$J_{zz}^{CM,C}$'; ...
        '$J_{xy}^{CM,C}$';'$J_{xz}^{CM,C}$';'$J_{yz}^{CM,C}$'};

    latex_units = {'mm';'mm';'mm'; ...
        '$\\mathrm{kg\\,m^2}$';'$\\mathrm{kg\\,m^2}$';'$\\mathrm{kg\\,m^2}$'; ...
        '$\\mathrm{kg\\,m^2}$';'$\\mathrm{kg\\,m^2}$';'$\\mathrm{kg\\,m^2}$'};

    for k = 1:height(T)
        fprintf(fid,'%s & %s & %.9g & %.9g & %.9g & %.9g & %.6g \\\\\\\\\n', ...
            latex_names{k},latex_units{k},T.Reference(k),T.Ideal(k), ...
            T.Noisy(k),T.AbsoluteError(k),T.RelativeError_pct(k));
    end

    fprintf(fid,'\\bottomrule\n\\end{tabular}}\n\\end{table*}\n');
end

function escribir_table_friction_tex(T,archivo)

    fid = abrir_archivo(archivo);
    c = onCleanup(@() fclose(fid));

    fprintf(fid,'\\begin{table}[t]\n\\centering\n');
    fprintf(fid,['\\caption{Equivalent viscous-friction coefficient estimated ', ...
        'during the two identification stages of the detailed InnoCube ', ...
        'validation.}\n']);
    fprintf(fid,'\\label{tab:friction_results}\n');
    fprintf(fid,'\\begin{tabular}{lrrr}\n\\toprule\n');
    fprintf(fid,'Stage & Reference & ideal & Noisy \\\\\\\\\n\\midrule\n');

    for k = 1:height(T)
        fprintf(fid,'%s & %.7f & %.7f & %.7f \\\\\\\\\n', ...
            T.Stage{k},T.Reference(k),T.Ideal(k),T.Noisy(k));
    end

    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n');
end

function escribir_table_multicubesat_tex(T,archivo)
    fid = fopen(archivo,'w');
    if fid < 0
        error('Could not create %s.',archivo);
    end

    fprintf(fid,'\\begin{table*}[t]\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,['\\caption{Cross-validation of the identification method for ', ...
        'five nanosatellite mass-property distributions. Center-of-mass ', ...
        'error is the Euclidean norm of the local CM-vector error. Inertia ', ...
        'error uses the complete centroidal tensor for InnoCube, TPA-1, and ', ...
        'InflateSail. For SOURCE and M3Sat, only the available diagonal ', ...
        'moments are included because source products of inertia are ', ...
        'unavailable.}\n']);
    fprintf(fid,'\\label{tab:multicubesat_validation}\n');
    fprintf(fid,'\\resizebox{\\textwidth}{!}{%%\n');
    fprintf(fid,'\\begin{tabular}{lrrrr}\n');
    fprintf(fid,'\\toprule\n');
    fprintf(fid,['Spacecraft & CM ideal [mm] & CM Noisy [mm] & ', ...
        'Inertia Ideal [\\%%] & Inertia Noisy [\\%%] \\\\\n']);
    fprintf(fid,'\\midrule\n');

    for i = 1:height(T)
        fprintf(fid,'%s & %.6g & %.6g & %.6g & %.6g \\\\\n', ...
            T.Spacecraft{i},T.CM_Ideal_mm(i),T.CM_Noisy_mm(i), ...
            T.Inertia_Ideal_pct(i),T.Inertia_Noisy_pct(i));
    end

    fprintf(fid,'\\bottomrule\n');
    fprintf(fid,'\\end{tabular}}\n');
    fprintf(fid,'\\end{table*}\n');
    fclose(fid);
end

function crear_macros_articulo(T,R,C,archivo)

    cm_noisy = norm(1e3*(R.medida.r_G_CM_C-C.r_G_CM_C));

    Jref = C.J_reportado_C;
    j_noisy = 100*norm(R.medida.J_CM_C-Jref,'fro')/norm(Jref,'fro');

    max_cm = max(T.CM_Noisy_mm);
    max_j = max(T.Inertia_Noisy_pct);

    fid = abrir_archivo(archivo);
    c = onCleanup(@() fclose(fid));

    fprintf(fid,'%% Generated automatically by MATLAB.\n');
    fprintf(fid,'\\renewcommand{\\InnoCMNoisy}{%.4g}\n',cm_noisy);
    fprintf(fid,'\\renewcommand{\\InnoJNoisy}{%.4g}\n',j_noisy);
    fprintf(fid,'\\renewcommand{\\MaxCMNoisy}{%.4g}\n',max_cm);
    fprintf(fid,'\\renewcommand{\\MaxJNoisy}{%.4g}\n',max_j);
    fprintf(fid,'\\renewcommand{\\NumValidationCases}{%d}\n',height(T));
end


%% ========================================================================
% UTILIDADES
% ========================================================================

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

function v = matriz_a_6(J)
    v = [J(1,1);J(2,2);J(3,3);J(1,2);J(1,3);J(2,3)];
end

function fid = abrir_archivo(nombre)
    fid = fopen(nombre,'w');
    if fid < 0
        error('Could not create the file: %s',nombre);
    end
end
