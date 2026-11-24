function R = estimate_platform_parameters(datos,P,A)
%ESTIMATE_PLATFORM_PARAMETERS
% First estimation: 10 parameters.
%
% theta =
% [Jxx Jyy Jzz Jxy Jxz Jyz Mrx Mry Mrz b]'
%
% In the complete sequence, the main quantities of interest from this estimate are:
%
% rCM_x
% rCM_y
% rCM_z
%
% It is computed using both ideal and noisy data.
%
% REACTION-WHEEL MODEL:
% The wheels are treated as internal actuators. Their stored angular momentum
% h_RW = J_RW .* Omega_RW is incorporated into Q.
% The reaction torque of the wheels is not longer added as an external torque
% in the integrated equation.

    if nargin < 2 || isempty(P)
        P = platform_parameters();
    end

    if nargin < 3 || isempty(A)
        A = actuator_torque_parameters();
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

    w_ruido = [ ...
        datos.wx, ...
        datos.wy, ...
        datos.wz ...
    ];

    gB_ideal = datos.gB_true;
    gB_ruido = datos.gB;

 % Known torque that remains on the right-hand side of the 10P model:
 % only the incremental gravitational torque due to the known motion
 % of the masses relative to the reference configuration.
    if isfield(datos,'tau_est1')
        tau_known = datos.tau_est1;
    else
        tau_known = datos.tau_masaX + datos.tau_masaY;
    end

 % Relative wheel speeds, known from the encoder/model.
    Omega_RW = datos.Omega_RW;

    dX = datos.dX;
    dY = datos.dY;
    dZ = datos.dZ;

    vX = datos.vX;
    vY = datos.vY;
    vZ = datos.vZ;

    dX_ref = datos.cfg.dX_bal;
    dY_ref = datos.cfg.dY_bal;
    dZ_ref = datos.cfg.dZ_bal;

    rX_ref = P.rX0 + P.eX*dX_ref;
    rY_ref = P.rY0 + P.eY*dY_ref;
    rZ_ref = P.rZ0 + P.eZ*dZ_ref;


 %% ============================================================
 % CHANGE PRODUCED BY THE FIRST BALANCING MANEUVER
 % ============================================================

    I3 = eye(3);

    DeltaJ_balanceo = ...
        cambio_inercia(P.mX,rX_ref,P.rX0,I3) ...
        + cambio_inercia(P.mY,rY_ref,P.rY0,I3) ...
        + cambio_inercia(P.mZ,rZ_ref,P.rZ0,I3);

    DeltaMr_balanceo = ...
        P.mX*(rX_ref-P.rX0) ...
        + P.mY*(rY_ref-P.rY0) ...
        + P.mZ*(rZ_ref-P.rZ0);


 %% ============================================================
 % ESTIMATES
 % ============================================================

    ideal = estimar10( ...
        t,w_ideal,gB_ideal,tau_known,Omega_RW, ...
        dX,dY,dZ,vX,vY,vZ, ...
        A.stride_est1,P, ...
        rX_ref,rY_ref,rZ_ref, ...
        DeltaJ_balanceo,DeltaMr_balanceo ...
    );

    ruido = estimar10( ...
        t,w_ruido,gB_ruido,tau_known,Omega_RW, ...
        dX,dY,dZ,vX,vY,vZ, ...
        A.stride_est1,P, ...
        rX_ref,rY_ref,rZ_ref, ...
        DeltaJ_balanceo,DeltaMr_balanceo ...
    );

 % =============================================================
 % CANONICAL CONFIGURATION
 % =============================================================
 %
 % to partir of this version, toda output destinada to comparisons
 % between campaigns expresa also in to configuration canonical
 % comun:
 %
 % dX = dY = dZ = 0
 %
 % The 10P estimation is performed around the balanced configuration,
 % balanceada, but to the final removes exactamente the efecto of the
 % is removed exactly at the end.
 %
 % the campos *_ini retain by compatibilidad; fisicamente
 % they represent this canonical configuration.

    ideal.J_hat_canonica = ideal.J_hat_ini;
    ruido.J_hat_canonica = ruido.J_hat_ini;

    ideal.S_hat_canonica = ideal.Mr_hat_ini;
    ruido.S_hat_canonica = ruido.Mr_hat_ini;

 % Alias with the nomenclature antigua.
    ideal.Mr_hat_canonica = ideal.S_hat_canonica;
    ruido.Mr_hat_canonica = ruido.S_hat_canonica;

    ideal.rCM_hat_canonica = ideal.S_hat_canonica/P.M_total;
    ruido.rCM_hat_canonica = ruido.S_hat_canonica/P.M_total;

    ideal.rCM_hat_ini = ideal.rCM_hat_canonica;
    ruido.rCM_hat_ini = ruido.rCM_hat_canonica;


 %% ============================================================
 % TRUE VALUES FOR VALIDATION
 % ============================================================

    theta_real = [ ...
        P.J0(1,1); ...
        P.J0(2,2); ...
        P.J0(3,3); ...
        P.J0(1,2); ...
        P.J0(1,3); ...
        P.J0(2,3); ...
        P.Mr0; ...
        P.b ...
    ];

    nombres = { ...
        'Jxx';'Jyy';'Jzz'; ...
        'Jxy';'Jxz';'Jyz'; ...
        'Mrx';'Mry';'Mrz';'b' ...
    };

    Tabla10P = table( ...
        nombres, ...
        theta_real, ...
        ideal.theta_hat_ini, ...
        ruido.theta_hat_ini, ...
        'VariableNames',{ ...
            'Parametro','Real','Ideal','Ruidosa' ...
        } ...
    );

    ejes = {'X';'Y';'Z'};

    TablaCM = table( ...
        ejes, ...
        1000*P.rCM0, ...
        1000*ideal.rCM_hat_ini, ...
        1000*ruido.rCM_hat_ini, ...
        'VariableNames',{ ...
            'Eje','Real_mm','Ideal_mm','Ruidosa_mm' ...
        } ...
    );


 %% ============================================================
 % display
 % ============================================================

    fprintf('\n');
    fprintf('=============================================================\n');
    fprintf(' first estimation - 10 parameters\n');
    fprintf('=============================================================\n');

    fprintf('\nRango ideal   = %d / 10\n',ideal.rank_Phi);
    fprintf('Cond. ideal   = %.6e\n',ideal.cond_Phi);
    fprintf('RMSE ideal    = %.6e\n',ideal.RMSE_batch);

    fprintf('\nRango noise   = %d / 10\n',ruido.rank_Phi);
    fprintf('Cond. noise   = %.6e\n',ruido.cond_Phi);
    fprintf('RMSE noise    = %.6e\n\n',ruido.RMSE_batch);

    disp(Tabla10P);

    fprintf('\nCANONICAL CONFIGURATION: dX = dY = dZ = 0\n');
    fprintf('ESTIMATED CANONICAL J - IDEAL [kg*m^2]\n');
    disp(ideal.J_hat_canonica);

    fprintf('ESTIMATED CANONICAL J - NOISY [kg*m^2]\n');
    disp(ruido.J_hat_canonica);

    fprintf('S CANONICO ESTIMADO - ideal [kg*m]\n');
    disp(ideal.S_hat_canonica);

    fprintf('S CANONICO ESTIMADO - noisy [kg*m]\n');
    disp(ruido.S_hat_canonica);

    fprintf('\nCANONICAL CENTER OF MASS\n');
    disp(TablaCM);


 %% ============================================================
 % OUTPUT
 % ============================================================

    R = struct;

    R.ideal = ideal;
    R.ruido = ruido;

    R.theta_real = theta_real;
    R.Tabla10P = Tabla10P;
    R.TablaCM = TablaCM;

    R.rCM_real = P.rCM0;

    R.DeltaJ_balanceo = DeltaJ_balanceo;
    R.DeltaMr_balanceo = DeltaMr_balanceo;

 % ------------------------------------------------------------
 % STANDARDIZED OUTPUT IN THE CANONICAL CONFIGURATION
 % ------------------------------------------------------------
 %
 % this estructura sera the that utilice posteriormente for
 % compare the BASE campaign with the BASE + CubeSat campaign.

    R.configuracion_canonica = struct;
    R.configuracion_canonica.d = [0;0;0];
    R.configuracion_canonica.descripcion = ...
        'Masas X/Y/Z en sus posiciones nominales: dX=dY=dZ=0';

    R.referencia = struct;
    R.referencia.d = [dX_ref;dY_ref;dZ_ref];
    R.referencia.DeltaJ_desde_canonica = DeltaJ_balanceo;
    R.referencia.DeltaS_desde_canonica = DeltaMr_balanceo;

    R.canonica = struct;

    R.canonica.ideal = struct;
    R.canonica.ideal.J = ideal.J_hat_canonica;
    R.canonica.ideal.S = ideal.S_hat_canonica;
    R.canonica.ideal.Mr = ideal.S_hat_canonica;
    R.canonica.ideal.rCM = ideal.rCM_hat_canonica;
    R.canonica.ideal.b = ideal.b_hat;

    R.canonica.ruido = struct;
    R.canonica.ruido.J = ruido.J_hat_canonica;
    R.canonica.ruido.S = ruido.S_hat_canonica;
    R.canonica.ruido.Mr = ruido.S_hat_canonica;
    R.canonica.ruido.rCM = ruido.rCM_hat_canonica;
    R.canonica.ruido.b = ruido.b_hat;

    R.canonica.real = struct;
    R.canonica.real.J = P.J0;
    R.canonica.real.S = P.Mr0;
    R.canonica.real.Mr = P.Mr0;
    R.canonica.real.rCM = P.rCM0;
    R.canonica.real.b = P.b;


    save(A.archivo_resultado1,'R');

end


function out = estimar10( ...
    t,w,gB,tau_known,Omega_RW, ...
    dX,dY,dZ,vX,vY,vZ, ...
    stride,P, ...
    rX_ref,rY_ref,rZ_ref, ...
    DeltaJ_balanceo,DeltaMr_balanceo ...
)

    N = length(t);
    I3 = eye(3);

    DeltaJ = zeros(3,3,N);
    h_rel = zeros(N,3);
    h_RW = Omega_RW.*P.J_RW(:)';

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

    end

 % Q_total agrupa all the moments angulares known that not
 % pertenecen to J_ref*w:
 %
 % Q_total = DeltaJ*w + h_rel + h_RW
    Q = zeros(N,3);

    for k = 1:N
        Q(k,:) = ( ...
            DeltaJ(:,:,k)*w(k,:)' ...
            + h_rel(k,:)' ...
            + h_RW(k,:)' ...
        )';
    end

    idx_est = 2:stride:N;
    numBlocks = numel(idx_est);

    Phi_B = zeros(3*numBlocks,10);
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
    g_prev = gB(1,:)';
    tau_prev = tau_known(1,:)';
    Q_prev = Q(1,:)';

    Omega0 = Omega(w_prev);
    Omega_prev = Omega0;

    A_prev = skew(w_prev)*Omega_prev;
    G_prev = skew(g_prev);
    C_prev = skew(w_prev)*Q_prev;

    I_wOmega = zeros(3,6);
    I_g = zeros(3,3);
    I_w = zeros(3,1);
    I_tau = zeros(3,1);
    I_wQ = zeros(3,1);

    block = 0;
    next_store = idx_est(1);

    for k = 2:N

        dt = t(k)-t(k-1);

        w_k = w(k,:)';
        g_k = gB(k,:)';
        tau_k = tau_known(k,:)';
        Q_k = Q(k,:)';

        Omega_k = Omega(w_k);

        A_k = skew(w_k)*Omega_k;
        G_k = skew(g_k);
        C_k = skew(w_k)*Q_k;

        I_wOmega = I_wOmega + 0.5*dt*(A_prev+A_k);
        I_g = I_g + 0.5*dt*(G_prev+G_k);
        I_w = I_w + 0.5*dt*(w_prev+w_k);
        I_tau = I_tau + 0.5*dt*(tau_prev+tau_k);
        I_wQ = I_wQ + 0.5*dt*(C_prev+C_k);

        if k == next_store

            block = block+1;

            Phi_J = (Omega_k-Omega0) + I_wOmega;
            Phi_mr = I_g;
            Phi_b = I_w;

            Phi_k = [Phi_J Phi_mr Phi_b];

            y_k = ...
                I_tau ...
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
        G_prev = G_k;

        tau_prev = tau_k;
        Q_prev = Q_k;
        C_prev = C_k;

    end

    theta_hat_bal = Phi_B\y_B;

    J_hat_bal = [ ...
        theta_hat_bal(1),theta_hat_bal(4),theta_hat_bal(5); ...
        theta_hat_bal(4),theta_hat_bal(2),theta_hat_bal(6); ...
        theta_hat_bal(5),theta_hat_bal(6),theta_hat_bal(3) ...
    ];

    Mr_hat_bal = theta_hat_bal(7:9);
    b_hat = theta_hat_bal(10);

    J_hat_ini = J_hat_bal-DeltaJ_balanceo;
    Mr_hat_ini = Mr_hat_bal-DeltaMr_balanceo;

    theta_hat_ini = [ ...
        J_hat_ini(1,1); ...
        J_hat_ini(2,2); ...
        J_hat_ini(3,3); ...
        J_hat_ini(1,2); ...
        J_hat_ini(1,3); ...
        J_hat_ini(2,3); ...
        Mr_hat_ini; ...
        b_hat ...
    ];

    [~,S,~] = svd(Phi_B,'econ');

    svals = diag(S);

    out = struct;

    out.theta_hat_bal = theta_hat_bal;
    out.theta_hat_ini = theta_hat_ini;

    out.J_hat_bal = J_hat_bal;
    out.J_hat_ini = J_hat_ini;

    out.Mr_hat_bal = Mr_hat_bal;
    out.Mr_hat_ini = Mr_hat_ini;

 % nomenclature explicit of reference/canonical.
    out.J_hat_ref = J_hat_bal;
    out.S_hat_ref = Mr_hat_bal;
    out.Mr_hat_ref = Mr_hat_bal;

    out.J_hat_canonica = J_hat_ini;
    out.S_hat_canonica = Mr_hat_ini;
    out.Mr_hat_canonica = Mr_hat_ini;
    out.rCM_hat_canonica = Mr_hat_ini/P.M_total;
    out.theta_hat_canonica = theta_hat_ini;

    out.b_hat = b_hat;

    out.rank_Phi = rank(Phi_B);

    if svals(end) > 0
        out.cond_Phi = svals(1)/svals(end);
    else
        out.cond_Phi = Inf;
    end

    out.singular_values = svals;

    residual = y_B-Phi_B*theta_hat_bal;
    out.RMSE_batch = sqrt(mean(residual.^2));

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
