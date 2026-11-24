function R = estimate_inertia_parameters(datos,P,A,Mr_conocido,d_ref)
%ESTIMATE_INERTIA_PARAMETERS
% Second estimation of:
%
% theta_J = [Jxx Jyy Jzz Jxy Jxz Jyz b]'
%
% DIFFERENCE WITH RESPECT TO THE PREVIOUS version
% ------------------------------------------
% During this second stage, the X/and masses are not LONGER stationary.
% They move rapidly to excite omega_x and omega_y.
%
% Therefore, the estimator explicitly considers:
%
% J(t) = J_ref + DeltaJ(t)
% Mr(t)
% h_rel(t)
% h_RW(t)
%
% The equation used by the plant is:
%
% J*w_dot = tau_g - b*w ...
% - J_dot*w - hdot_rel - hdot_RW ...
% - w x (J*w + h_rel + h_RW)
%
% with:
%
% tau_g = Mr x g = -[g]x Mr
%
% Defining:
%
% J(t) = J_ref + DeltaJ(t)
% Q(t) = DeltaJ(t)*w(t) + h_rel(t) + h_RW(t)
%
% we obtain:
%
% J_ref*w_dot + w x (J_ref*w) + b*w ...
% + [g]x*Mr(t) + dQ/dt + w x Q = 0
%
% and, after integration, there is not need to numerically differentiate omega or the
% mass velocities.
%
% The reference configuration is the one that exists immediately before starting
% the second excitation. d_ref = [dX_ref dY_ref dZ_ref]'.
%
% The selector is retained:
%
% to.modo_datos_est2 = 'noisy'
% to.modo_datos_est2 = 'filtered'

    if nargin < 2 || isempty(P)
        P = platform_parameters();
    end

    if nargin < 3 || isempty(A)
        A = actuator_torque_parameters();
    end

    modo_est2 = lower(strtrim(A.modo_datos_est2));

    switch modo_est2

        case 'ruidosa'
            etiqueta_medida = 'Ruidosa';

        case 'filtrada'
            etiqueta_medida = 'Ruidosa + SG';

        otherwise
            error('Unknown second-estimation data mode: %s', ...
                A.modo_datos_est2);

    end


 %% ============================================================
 % DATA
 % ============================================================

    t = datos.t-datos.t(1);

    w_ideal = [ ...
        datos.wx_true, ...
        datos.wy_true, ...
        datos.wz_true ...
    ];

    w_medida = [ ...
        datos.wx, ...
        datos.wy, ...
        datos.wz ...
    ];

    gB_ideal = datos.gB_true;
    gB_medida = datos.gB;

 % Relative wheel speeds. They are treated as known
 % from the encoder/model and are not numerically differentiated.
    Omega_RW = datos.Omega_RW;

    dX = datos.dX;
    dY = datos.dY;
    dZ = datos.dZ;

    vX = datos.vX;
    vY = datos.vY;
    vZ = datos.vZ;


 %% ============================================================
 % 7P ESTIMATES WITH MOVING MASSES
 % ============================================================

    ideal = estimar7_masas_moviles( ...
        t,w_ideal,gB_ideal,Omega_RW, ...
        dX,dY,dZ,vX,vY,vZ, ...
        Mr_conocido.ideal,d_ref,A.stride_est2,P ...
    );

    medida = estimar7_masas_moviles( ...
        t,w_medida,gB_medida,Omega_RW, ...
        dX,dY,dZ,vX,vY,vZ, ...
        Mr_conocido.ruido,d_ref,A.stride_est2,P ...
    );


 %% ============================================================
 % REFERENCE CONFIGURATION AND RECONSTRUCTION OF original J
 % ============================================================

    I3 = eye(3);

    rX_ref = P.rX0 + P.eX*d_ref(1);
    rY_ref = P.rY0 + P.eY*d_ref(2);
    rZ_ref = P.rZ0 + P.eZ*d_ref(3);

    DeltaJ_ref = ...
        cambio_inercia(P.mX,rX_ref,P.rX0,I3) ...
        + cambio_inercia(P.mY,rY_ref,P.rY0,I3) ...
        + cambio_inercia(P.mZ,rZ_ref,P.rZ0,I3);

    DeltaS_ref = ...
        P.mX*(rX_ref-P.rX0) ...
        + P.mY*(rY_ref-P.rY0) ...
        + P.mZ*(rZ_ref-P.rZ0);

 % =============================================================
 % CANONICAL CONFIGURATION
 % =============================================================
 %
 % The common canonical configuration for all campaigns is:
 %
 % dX = dY = dZ = 0
 %
 % The 7P estimator identifies J_ref in the configuration that exists
 % immediately before the second excitation. To enable future comparison of to
 % BASE campaign with to BASE + CubeSat campaign, J is always reconstructed
 % in the same canonical configuration.
 %
 % The 7P estimator does NOT estimate S again. The first moment used
 % in this stage comes from the 10P estimate and is transported
 % algebraically from the reference configuration to d=0.

    ideal.J_hat_canonica = ...
        ideal.J_hat_ref-DeltaJ_ref;

    medida.J_hat_canonica = ...
        medida.J_hat_ref-DeltaJ_ref;

    ideal.S_hat_ref = Mr_conocido.ideal;
    medida.S_hat_ref = Mr_conocido.ruido;

    ideal.S_hat_canonica = ...
        ideal.S_hat_ref-DeltaS_ref;

    medida.S_hat_canonica = ...
        medida.S_hat_ref-DeltaS_ref;

 % Alias using the historical nomenclature.
    ideal.Mr_hat_ref = ideal.S_hat_ref;
    medida.Mr_hat_ref = medida.S_hat_ref;

    ideal.Mr_hat_canonica = ideal.S_hat_canonica;
    medida.Mr_hat_canonica = medida.S_hat_canonica;

    ideal.rCM_hat_canonica = ...
        ideal.S_hat_canonica/P.M_total;

    medida.rCM_hat_canonica = ...
        medida.S_hat_canonica/P.M_total;

 % Alias retained for compatibility with previous scripts.
    ideal.J_hat_original = ideal.J_hat_canonica;
    medida.J_hat_original = medida.J_hat_canonica;

    ideal.Mr_hat_original = ideal.S_hat_canonica;
    medida.Mr_hat_original = medida.S_hat_canonica;

    ideal.J_hat_compensada = ideal.J_hat_ref;
    medida.J_hat_compensada = medida.J_hat_ref;


 %% ============================================================
 % TRUE VALUES FOR VALIDATION
 % ============================================================

    J_real_ref = P.J0 + DeltaJ_ref;
    J_real_original = P.J0;

    S_real_ref = P.Mr0 + DeltaS_ref;
    S_real_canonica = P.Mr0;
    rCM_real_canonica = P.rCM0;


 %% ============================================================
 % TABLES
 % ============================================================

    nombres = { ...
        'Jxx';'Jyy';'Jzz'; ...
        'Jxy';'Jxz';'Jyz'; ...
        'b' ...
    };

    theta_real_ref = [ ...
        J_real_ref(1,1); ...
        J_real_ref(2,2); ...
        J_real_ref(3,3); ...
        J_real_ref(1,2); ...
        J_real_ref(1,3); ...
        J_real_ref(2,3); ...
        P.b ...
    ];

    TablaReferencia = table( ...
        nombres, ...
        theta_real_ref, ...
        ideal.theta_hat, ...
        medida.theta_hat, ...
        'VariableNames',{ ...
            'Parametro','Real','Ideal','Medida' ...
        } ...
    );

    theta_real_original = matriz_a_6(J_real_original);
    theta_J_orig_ideal = matriz_a_6(ideal.J_hat_original);
    theta_J_orig_medida = matriz_a_6(medida.J_hat_original);

    nombresJ = { ...
        'Jxx';'Jyy';'Jzz'; ...
        'Jxy';'Jxz';'Jyz' ...
    };

    TablaOriginal = table( ...
        nombresJ, ...
        theta_real_original, ...
        theta_J_orig_ideal, ...
        theta_J_orig_medida, ...
        'VariableNames',{ ...
            'Parametro','Real','Ideal','Medida' ...
        } ...
    );

    ejes = {'X';'Y';'Z'};

    TablaSCanonica = table( ...
        ejes, ...
        S_real_canonica, ...
        ideal.S_hat_canonica, ...
        medida.S_hat_canonica, ...
        'VariableNames',{ ...
            'Eje','Real','Ideal','Medida' ...
        } ...
    );

    TablaCMCanonico = table( ...
        ejes, ...
        1000*rCM_real_canonica, ...
        1000*ideal.rCM_hat_canonica, ...
        1000*medida.rCM_hat_canonica, ...
        'VariableNames',{ ...
            'Eje','Real_mm','Ideal_mm','Medida_mm' ...
        } ...
    );


 %% ============================================================
 % DISPLAY RESULTS
 % ============================================================

    fprintf('\n');
    fprintf('=============================================================\n');
    fprintf(' second estimation - MOVING MASSES + moment angular RW\n');
    fprintf('=============================================================\n');

    fprintf('\nRango ideal   = %d / 7\n',ideal.rank_Phi);
    fprintf('Cond. ideal   = %.6e\n',ideal.cond_Phi);
    fprintf('RMSE ideal    = %.6e\n',ideal.RMSE_batch);

    fprintf('\nValores singulares ideal:\n');
    disp(ideal.singular_values.');

    if ~isempty(ideal.singular_values) && ideal.singular_values(1) > 0
        fprintf('values singulares ideal normalizados:\n');
        disp((ideal.singular_values./ideal.singular_values(1)).');
    end

    fprintf('\nData used in measured estimation: %s\n',etiqueta_medida);
    fprintf('Rango measured = %d / 7\n',medida.rank_Phi);
    fprintf('Cond. measured = %.6e\n',medida.cond_Phi);
    fprintf('RMSE measured  = %.6e\n',medida.RMSE_batch);

    fprintf('\nValores singulares measured:\n');
    disp(medida.singular_values.');

    if ~isempty(medida.singular_values) && medida.singular_values(1) > 0
        fprintf('values singulares measured normalizados:\n');
        disp((medida.singular_values./medida.singular_values(1)).');
    end

    fprintf('\nSECOND-STAGE REFERENCE CONFIGURATION\n');
    disp(TablaReferencia);

    fprintf('\nCANONICAL CONFIGURATION: dX = dY = dZ = 0\n');
    disp(TablaOriginal);

    fprintf('\nTRUE CANONICAL J [kg*m^2]\n');
    disp(J_real_original);

    fprintf('ESTIMATED CANONICAL J - IDEAL [kg*m^2]\n');
    disp(ideal.J_hat_canonica);

    fprintf('ESTIMATED CANONICAL J - %s [kg*m^2]\n',upper(etiqueta_medida));
    disp(medida.J_hat_canonica);

    fprintf('\nS CANONICO [kg*m]\n');
    disp(TablaSCanonica);

    fprintf('\nCANONICAL CENTER OF MASS [mm]\n');
    disp(TablaCMCanonico);


 %% ============================================================
 % PLOT
 % ============================================================

    figure( ...
        'Name','Segunda estimacion - matriz de inercia', ...
        'NumberTitle','off' ...
    );

    bar([ ...
        theta_real_original, ...
        theta_J_orig_ideal, ...
        theta_J_orig_medida ...
    ]);

    grid on;

    xticks(1:6);
    xticklabels(nombresJ);

    xlabel('Parametro');
    ylabel('Inercia [kg*m^2]');

    legend('Real','Ideal',etiqueta_medida,'Location','best');

    title('Segunda estimacion con masas moviles y momento angular de ruedas');


 %% ============================================================
 % OUTPUT
 % ============================================================

    R = struct;

    R.ideal = ideal;
    R.ruido = medida;       % nombre conservado por compatibilidad
    R.medida = medida;

    R.J_real_referencia = J_real_ref;
    R.J_real_compensada = J_real_ref;   % alias
    R.J_real_original = J_real_original;

    R.S_real_referencia = S_real_ref;
    R.S_real_canonica = S_real_canonica;
    R.Mr_real_canonica = S_real_canonica;

    R.DeltaJ_ref = DeltaJ_ref;
    R.DeltaJ_final = DeltaJ_ref;         % alias
    R.DeltaS_ref = DeltaS_ref;
    R.DeltaMr_ref = DeltaS_ref;          % alias

    R.TablaReferencia = TablaReferencia;
    R.TablaCompensada = TablaReferencia; % alias
    R.TablaOriginal = TablaOriginal;
    R.TablaSCanonica = TablaSCanonica;
    R.TablaCMCanonico = TablaCMCanonico;

    R.Mr_conocido = Mr_conocido;
    R.d_ref = d_ref;
    R.d_final = d_ref;                   % alias

 % ------------------------------------------------------------
 % STANDARDIZED OUTPUT IN THE CANONICAL CONFIGURATION
 % ------------------------------------------------------------

    R.configuracion_canonica = struct;
    R.configuracion_canonica.d = [0;0;0];
    R.configuracion_canonica.descripcion = ...
        'Masas X/Y/Z en sus posiciones nominales: dX=dY=dZ=0';

    R.referencia = struct;
    R.referencia.d = d_ref(:);
    R.referencia.DeltaJ_desde_canonica = DeltaJ_ref;
    R.referencia.DeltaS_desde_canonica = DeltaS_ref;

    R.canonica = struct;

    R.canonica.ideal = struct;
    R.canonica.ideal.J = ideal.J_hat_canonica;
    R.canonica.ideal.S = ideal.S_hat_canonica;
    R.canonica.ideal.Mr = ideal.S_hat_canonica;
    R.canonica.ideal.rCM = ideal.rCM_hat_canonica;
    R.canonica.ideal.b = ideal.b_hat;

    R.canonica.ruido = struct;
    R.canonica.ruido.J = medida.J_hat_canonica;
    R.canonica.ruido.S = medida.S_hat_canonica;
    R.canonica.ruido.Mr = medida.S_hat_canonica;
    R.canonica.ruido.rCM = medida.rCM_hat_canonica;
    R.canonica.ruido.b = medida.b_hat;

 % Same content with to generic name, useful when the data have been
 % filtered and it is not longer appropriate to call them strictly 'noise'.
    R.canonica.medida = R.canonica.ruido;

    R.canonica.real = struct;
    R.canonica.real.J = J_real_original;
    R.canonica.real.S = S_real_canonica;
    R.canonica.real.Mr = S_real_canonica;
    R.canonica.real.rCM = rCM_real_canonica;
    R.canonica.real.b = P.b;

    R.canonica.fuente_S = ...
        'S proviene de la primera estimacion 10P y se transporta a d=0';

    R.modo_datos_est2 = modo_est2;
    R.etiqueta_medida = etiqueta_medida;

    if isfield(datos,'filtrado_sg')
        R.filtrado_sg = datos.filtrado_sg;
    else
        R.filtrado_sg = struct('aplicado',false);
    end

    save(A.archivo_resultado2,'R');

end


function out = estimar7_masas_moviles( ...
    t,w,gB,Omega_RW, ...
    dX,dY,dZ,vX,vY,vZ, ...
    Mr_ref,d_ref,stride,P)
%ESTIMAR7_MASAS_MOVILES
%
% Estimates J_ref and b assuming the following are known:
%
% - mass displacements
% - mass velocities
% - variation of Mr(t)
% - variation DeltaJ(t)
% - relative angular momentum h_rel(t)
% - angular momentum stored in the wheels h_RW(t)
%
% Integrated equation:
%
% Phi_J*theta_J + Phi_b*b =
% - I_gMr
% - (Q(t)-Q(0))
% - I_wQ
%
% where:
%
% Q = DeltaJ*w + h_rel + h_RW
% I_gMr = integral( [g]x Mr(t) dt )
% I_wQ = integral( [w]x Q dt )

    N = length(t);
    I3 = eye(3);

    if N < 2
        error('The second stage requires at least two samples.');
    end

    rX_ref = P.rX0 + P.eX*d_ref(1);
    rY_ref = P.rY0 + P.eY*d_ref(2);
    rZ_ref = P.rZ0 + P.eZ*d_ref(3);

    DeltaJ = zeros(3,3,N);
    h_rel = zeros(N,3);
    h_RW = Omega_RW.*P.J_RW(:)';
    Mr_t = zeros(N,3);
    Q = zeros(N,3);

    for k = 1:N

        rX = P.rX0 + P.eX*dX(k);
        rY = P.rY0 + P.eY*dY(k);
        rZ = P.rZ0 + P.eZ*dZ(k);

        rX_dot = P.eX*vX(k);
        rY_dot = P.eY*vY(k);
        rZ_dot = P.eZ*vZ(k);

        DeltaJ(:,:,k) = ...
            cambio_inercia(P.mX,rX,rX_ref,I3) ...
            + cambio_inercia(P.mY,rY,rY_ref,I3) ...
            + cambio_inercia(P.mZ,rZ,rZ_ref,I3);

        h_rel(k,:) = ( ...
            P.mX*cross(rX,rX_dot) ...
            + P.mY*cross(rY,rY_dot) ...
            + P.mZ*cross(rZ,rZ_dot) ...
        )';

        DeltaMr = ...
            P.mX*(rX-rX_ref) ...
            + P.mY*(rY-rY_ref) ...
            + P.mZ*(rZ-rZ_ref);

        Mr_t(k,:) = (Mr_ref + DeltaMr)';

        Q(k,:) = ( ...
            DeltaJ(:,:,k)*w(k,:)' ...
            + h_rel(k,:)' ...
            + h_RW(k,:)' ...
        )';

    end


 %% ============================================================
 % INTEGRATED REGRESSOR
 % ============================================================

    idx_est = 2:stride:N;
    numBlocks = numel(idx_est);

    Phi_B = zeros(3*numBlocks,7);
    y_B = zeros(3*numBlocks,1);

    skew = @(a) [ ...
         0,    -a(3),  a(2); ...
         a(3),  0,    -a(1); ...
        -a(2),  a(1),  0 ...
    ];

    Omega = @(ww) [ ...
        ww(1),0,0,ww(2),ww(3),0; ...
        0,ww(2),0,ww(1),0,ww(3); ...
        0,0,ww(3),0,ww(1),ww(2) ...
    ];

    w_prev = w(1,:)';
    Q_prev = Q(1,:)';

    Omega0 = Omega(w_prev);
    Omega_prev = Omega0;

    A_prev = skew(w_prev)*Omega_prev;
    gMr_prev = skew(gB(1,:)')*Mr_t(1,:)';
    C_prev = skew(w_prev)*Q_prev;

    I_wOmega = zeros(3,6);
    I_w = zeros(3,1);
    I_gMr = zeros(3,1);
    I_wQ = zeros(3,1);

    block = 0;
    next_store = idx_est(1);

    for k = 2:N

        dt = t(k)-t(k-1);

        w_k = w(k,:)';
        Q_k = Q(k,:)';

        Omega_k = Omega(w_k);

        A_k = skew(w_k)*Omega_k;
        gMr_k = skew(gB(k,:)')*Mr_t(k,:)';
        C_k = skew(w_k)*Q_k;

        I_wOmega = I_wOmega + 0.5*dt*(A_prev+A_k);
        I_w = I_w + 0.5*dt*(w_prev+w_k);
        I_gMr = I_gMr + 0.5*dt*(gMr_prev+gMr_k);
        I_wQ = I_wQ + 0.5*dt*(C_prev+C_k);

        if k == next_store

            block = block+1;

            Phi_J = (Omega_k-Omega0) + I_wOmega;
            Phi_b = I_w;

            Phi_k = [Phi_J Phi_b];

            y_k = ...
                - I_gMr ...
                - (Q_k-Q(1,:)') ...
                - I_wQ;

            rows = (3*block-2):(3*block);

            Phi_B(rows,:) = Phi_k;
            y_B(rows) = y_k;

            if block < numBlocks
                next_store = idx_est(block+1);
            end

        end

        w_prev = w_k;
        Omega_prev = Omega_k;

        A_prev = A_k;
        Q_prev = Q_k;

        gMr_prev = gMr_k;
        C_prev = C_k;

    end


 %% ============================================================
 % SOLUCION
 % ============================================================

    theta_hat = Phi_B\y_B;

    J_hat_ref = [ ...
        theta_hat(1),theta_hat(4),theta_hat(5); ...
        theta_hat(4),theta_hat(2),theta_hat(6); ...
        theta_hat(5),theta_hat(6),theta_hat(3) ...
    ];

    [~,S,~] = svd(Phi_B,'econ');
    svals = diag(S);

    out = struct;

    out.theta_hat = theta_hat;

    out.J_hat_ref = J_hat_ref;
    out.b_hat = theta_hat(7);

    out.rank_Phi = rank(Phi_B);

    if isempty(svals) || svals(end) <= eps(max(svals))
        out.cond_Phi = Inf;
    else
        out.cond_Phi = svals(1)/svals(end);
    end

    out.singular_values = svals;

    residual = y_B-Phi_B*theta_hat;
    out.RMSE_batch = sqrt(mean(residual.^2));

 % diagnostic information.
    out.Mr_inicio = Mr_t(1,:)';
    out.Mr_final = Mr_t(end,:)';

    out.max_DeltaJ_F = 0;
    for k = 1:N
        out.max_DeltaJ_F = max( ...
            out.max_DeltaJ_F,norm(DeltaJ(:,:,k),'fro') ...
        );
    end

    out.max_h_rel = max(sqrt(sum(h_rel.^2,2)));
    out.max_h_RW = max(sqrt(sum(h_RW.^2,2)));
    out.Omega_RW_inicial = Omega_RW(1,:)';
    out.Omega_RW_final = Omega_RW(end,:)';

end


function DJ = cambio_inercia(m,r_new,r_ref,I3)

    DJ = m*( ...
        ((r_new'*r_new)-(r_ref'*r_ref))*I3 ...
        - (r_new*r_new'-r_ref*r_ref') ...
    );

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
