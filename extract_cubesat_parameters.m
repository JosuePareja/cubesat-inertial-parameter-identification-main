function R = extract_cubesat_parameters(BASE,CON_CUBESAT,C)
%EXTRACT_CUBESAT_PARAMETERS
% Recovers CubeSat parameters by subtracting two campaigns expressed in
% the same canonical configuration dX=dY=dZ=0:
%
% CubeSat = (BASE + CubeSat) - BASE
%
% The CubeSat mass is assumed known.
%
% Main outputs:
% R.measured.r_CR_CM_B
% R.measured.r_G_CM_C
% R.measured.J_CR_B
% R.measured.J_CM_B
% R.measured.J_CM_C
%
% The final comparison is performed in the original CubeSat C axes.

    validar_campana(BASE,'BASE');
    validar_campana(CON_CUBESAT,'BASE+CubeSat');

    modos = {'ideal','medida','real'};

    R = struct;
    R.nombre = C.nombre;
    R.id = C.id;
    R.M = C.M;
    R.C = C;

    for k = 1:numel(modos)

        modo = modos{k};

        B = BASE.canonica_final.(modo);
        BC = CON_CUBESAT.canonica_final.(modo);

        S_sat_B = BC.S-B.S;
        r_CR_CM_B = S_sat_B/C.M;

        J_CR_B = BC.J-B.J;

        I3 = eye(3);
        J_CM_B = J_CR_B - ...
            C.M*((r_CR_CM_B'*r_CR_CM_B)*I3-r_CR_CM_B*r_CR_CM_B');

        r_G_CM_B = r_CR_CM_B-C.r_CR_G_B;
        r_G_CM_C = C.R_BC'*r_G_CM_B;
        J_CM_C = C.R_BC'*J_CM_B*C.R_BC;

        R.(modo) = struct;
        R.(modo).S_CR_B = S_sat_B;
        R.(modo).r_CR_CM_B = r_CR_CM_B;
        R.(modo).r_G_CM_B = r_G_CM_B;
        R.(modo).r_G_CM_C = r_G_CM_C;
        R.(modo).J_CR_B = J_CR_B;
        R.(modo).J_CM_B = J_CM_B;
        R.(modo).J_CM_C = J_CM_C;

    end

 % =============================================================
 % Center-of-mass comparison
 % =============================================================

    ejes = {'x_C';'y_C';'z_C'};

    cm_original_mm = 1000*C.r_G_CM_C;
    cm_ideal_mm = 1000*R.ideal.r_G_CM_C;
    cm_medida_mm = 1000*R.medida.r_G_CM_C;

    error_cm_ideal_mm = cm_ideal_mm-cm_original_mm;
    error_cm_medida_mm = cm_medida_mm-cm_original_mm;

    R.TablaCM = table( ...
        ejes, ...
        cm_original_mm, ...
        cm_ideal_mm, ...
        cm_medida_mm, ...
        error_cm_ideal_mm, ...
        error_cm_medida_mm, ...
        'VariableNames',{ ...
            'Eje','Original_mm','Ideal_mm','Estimada_mm', ...
            'ErrorIdeal_mm','ErrorEstimada_mm' ...
        } ...
    );

    R.error_CM_norma_ideal_mm = norm(error_cm_ideal_mm);
    R.error_CM_norma_medida_mm = norm(error_cm_medida_mm);

 % =============================================================
 % Inertia-tensor comparison in C axes
 % =============================================================

    nombresJ = {'Jxx';'Jyy';'Jzz';'Jxy';'Jxz';'Jyz'};

    J_original_6 = matriz_a_6_con_nan(C.J_reportado_C);
    J_ideal_6 = matriz_a_6(R.ideal.J_CM_C);
    J_medida_6 = matriz_a_6(R.medida.J_CM_C);

    error_J_ideal = J_ideal_6-J_original_6;
    error_J_medida = J_medida_6-J_original_6;

    error_J_ideal_pct = porcentaje_error(J_ideal_6,J_original_6);
    error_J_medida_pct = porcentaje_error(J_medida_6,J_original_6);

    R.TablaJ = table( ...
        nombresJ, ...
        J_original_6, ...
        J_ideal_6, ...
        J_medida_6, ...
        error_J_ideal, ...
        error_J_medida, ...
        error_J_ideal_pct, ...
        error_J_medida_pct, ...
        'VariableNames',{ ...
            'Parametro','Original_kgm2','Ideal_kgm2','Estimada_kgm2', ...
            'ErrorIdeal_kgm2','ErrorEstimada_kgm2', ...
            'ErrorIdeal_pct','ErrorEstimada_pct' ...
        } ...
    );

 % Global inertia error is evaluated only over source-supported reference
 % components. This prevents unreported products of inertia from entering the
 % cross-validation metric for SOURCE and M3Sat.
    maskJ = isfinite(C.J_reportado_C);
    Jref_available = C.J_reportado_C(maskJ);
    R.error_J_F_ideal_pct = ...
        100*norm(R.ideal.J_CM_C(maskJ)-Jref_available)/norm(Jref_available);

    R.error_J_F_medida_pct = ...
        100*norm(R.medida.J_CM_C(maskJ)-Jref_available)/norm(Jref_available);

 % =============================================================
 % Friction coefficient - diagnostic quantity, not a CubeSat parameter
 % =============================================================

    R.friccion = struct;
    R.friccion.BASE_10P = BASE.R1.canonica.ruido.b;
    R.friccion.CON_CUBESAT_10P = CON_CUBESAT.R1.canonica.ruido.b;
    R.friccion.BASE_7P = BASE.R2.canonica.medida.b;
    R.friccion.CON_CUBESAT_7P = CON_CUBESAT.R2.canonica.medida.b;

 % =============================================================
 % DISPLAY RESULTS
 % =============================================================

    fprintf('\n');
    fprintf('=============================================================\n');
    fprintf(' CUBESAT IDENTIFICATION: %s\n',C.nombre);
    fprintf('=============================================================\n');

    fprintf('\nKnown mass = %.6f kg\n',C.M);

    fprintf('\nLOCAL CUBESAT CM IN C AXES [mm]\n');
    disp(R.TablaCM);

    fprintf('Euclidean CM error - ideal   = %.6f mm\n', ...
        R.error_CM_norma_ideal_mm);
    fprintf('Euclidean CM error - noisy     = %.6f mm\n', ...
        R.error_CM_norma_medida_mm);

    fprintf('\nCUBESAT INERTIA tensor ABOUT ITS CM, C AXES\n');
    disp(R.TablaJ);

    fprintf('error Frobenius J - ideal  = %.6f %%\n',R.error_J_F_ideal_pct);
    fprintf('error Frobenius J - noisy    = %.6f %%\n',R.error_J_F_medida_pct);

    if ~C.productos_inercia_reportados
        fprintf([ ...
            '\nNOTE: source products of inertia are unavailable for %s. ', ...
            'Zero products are used only to close the numerical plant; ', ...
            'unavailable source-comparison entries remain NaN.\n' ...
        ],C.nombre);
    end


    fprintf('\nEstimated friction coefficient [N*m*s/rad]\n');
    fprintf('  BASE 10P           = %.9g\n',R.friccion.BASE_10P);
    fprintf('  BASE+CubeSat 10P   = %.9g\n',R.friccion.CON_CUBESAT_10P);
    fprintf('  BASE 7P            = %.9g\n',R.friccion.BASE_7P);
    fprintf('  BASE+CubeSat 7P    = %.9g\n',R.friccion.CON_CUBESAT_7P);
    fprintf('  b is not subtracted; it belongs to the test system.\n');

    fprintf('=============================================================\n');

 % =============================================================
 % PLOTS of comparison
 % =============================================================

    figure( ...
        'Name',['CubeSat ' C.nombre ' - CM comparison'], ...
        'NumberTitle','off' ...
    );

    bar([cm_original_mm,cm_ideal_mm,cm_medida_mm]);
    grid on;
    xticks(1:3);
    xticklabels(ejes);
    ylabel('Position [mm]');
    legend('Reference','Ideal','Noisy','Location','best');
    title(['Local center of mass - ' C.nombre]);

    figure( ...
        'Name',['CubeSat ' C.nombre ' - inertia comparison'], ...
        'NumberTitle','off' ...
    );

    bar([J_original_6,J_ideal_6,J_medida_6]);
    grid on;
    xticks(1:6);
    xticklabels(nombresJ);
    ylabel('Inertia [kg*m^2]');
    legend('Reference','Ideal','Noisy','Location','best');
    title(['Tensor about the CM - ' C.nombre]);

end


function validar_campana(S,nombre)

    if ~isstruct(S) || ~isfield(S,'canonica_final')
        error('Campaign %s does not contain canonica_final.',nombre);
    end

    campos = {'ideal','medida','real'};

    for k = 1:numel(campos)
        if ~isfield(S.canonica_final,campos{k})
            error('Campaign %s does not contain mode %s.', ...
                nombre,campos{k});
        end
    end

end


function v = matriz_a_6(J)

    v = [ ...
        J(1,1); ...
        J(2,2); ...
        J(3,3); ...
        J(1,2); ...
        J(1,3); ...
        J(2,3) ...
    ];

end


function v = matriz_a_6_con_nan(J)
% Same as matriz_a_6; the explicit name reminds that J may contain
% NaN values when source products of inertia are unavailable.

    v = matriz_a_6(J);

end


function e = porcentaje_error(estimado,original)

    e = nan(size(original));

    valido = isfinite(original) & abs(original) > eps;
    e(valido) = 100*abs(estimado(valido)-original(valido))./abs(original(valido));

end
