function datos = simulate_platform(tspan,x0,cfg,P,A)
%SIMULATE_PLATFORM
% Runs ONE stage of the experimental sequence.
%
% Coordination between stages is performed in:
% secuencia_experimental.m
%
% Inputs:
% tspan : vector of absolute times
% x0 : initial condition
% cfg : stage configuration
% P : platform_parameters()
% to : actuator_torque_parameters()
%
% Output:
% data : states, attitude, gravity, CM, Omega_RW, h_RW, and torques.
    if nargin < 4 || isempty(P)
        P = platform_parameters();
    end

    if nargin < 5 || isempty(A)
        A = actuator_torque_parameters();
    end


 %% ============================================================
 % INTEGRATION
 % ============================================================

    opts = odeset( ...
        'RelTol',1e-8, ...
        'AbsTol',1e-10, ...
        'MaxStep',A.dt ...
    );

    [t,x] = ode45( ...
        @(tt,xx) platform_dynamics_ode(tt,xx,cfg,P,A), ...
        tspan, ...
        x0, ...
        opts ...
    );

    N = length(t);


 %% ============================================================
 % STATES
 % ============================================================

    q = x(:,1:4);

    qnorm = sqrt(sum(q.^2,2));
    q = q./qnorm;

    qx = q(:,1);
    qy = q(:,2);
    qz = q(:,3);
    qw = q(:,4);

    wx_true = x(:,5);
    wy_true = x(:,6);
    wz_true = x(:,7);

    dX = x(:,8);
    dY = x(:,9);
    dZ = x(:,10);

    vX = x(:,11);
    vY = x(:,12);
    vZ = x(:,13);

 % integral states of the X/and position PI controller.
    xiX = x(:,14);
    xiY = x(:,15);

 % Relative reaction-wheel speeds.
    Omega_RW = x(:,16:18);

 % Relative angular momentum stored in the wheels.
    h_RW = Omega_RW.*P.J_RW(:)';


 %% ============================================================
 % EULER ANGLES
 % ============================================================

    sin_theta = 2.*(qw.*qy-qz.*qx);
    sin_theta = max(-1,min(1,sin_theta));

    phi = atan2( ...
        2.*(qw.*qx+qy.*qz), ...
        1-2.*(qx.^2+qy.^2) ...
    );

    theta = asin(sin_theta);

    psi = atan2( ...
        2.*(qw.*qz+qx.*qy), ...
        1-2.*(qy.^2+qz.^2) ...
    );


 %% ============================================================
 % ideal GRAVITY
 % ============================================================

    gx_true = -P.g.*2.*(qx.*qz-qy.*qw);
    gy_true = -P.g.*2.*(qy.*qz+qx.*qw);
    gz_true = -P.g.*(1-2.*(qx.^2+qy.^2));

    gB_true = [gx_true gy_true gz_true];


 %% ============================================================
 % NOISY MEASUREMENTS
 % ============================================================

    wx = wx_true + A.sigma_gyro.*randn(N,1);
    wy = wy_true + A.sigma_gyro.*randn(N,1);
    wz = wz_true + A.sigma_gyro.*randn(N,1);

    gx = gx_true + A.sigma_g.*randn(N,1);
    gy = gy_true + A.sigma_g.*randn(N,1);
    gz = gz_true + A.sigma_g.*randn(N,1);

    gB = [gx gy gz];


 %% ============================================================
 % TRUE CENTER OF MASS
 % ============================================================

    rCM = zeros(N,3);

    for k = 1:N

        rX = P.rX0 + P.eX*dX(k);
        rY = P.rY0 + P.eY*dY(k);
        rZ = P.rZ0 + P.eZ*dZ(k);

        Mr = ...
            P.M_p*P.r_P ...
            + P.mX*rX ...
            + P.mY*rY ...
            + P.mZ*rZ;

        rCM(k,:) = (Mr/P.M_total)';

    end


 %% ============================================================
 % KNOWN torques
 % ============================================================

    tau_RW = zeros(N,3);
    tau_masaX = zeros(N,3);
    tau_masaY = zeros(N,3);

    for k = 1:N

        if strcmpi(cfg.modo,'excitacion_1')

            tau_RW(k,:) = A.torque_RW1(t(k))';

            DeltaMr_X = ...
                P.mX*(dX(k)-cfg.dX_bal)*P.eX;

            DeltaMr_Y = ...
                P.mY*(dY(k)-cfg.dY_bal)*P.eY;

            tau_masaX(k,:) = ...
                cross(DeltaMr_X,gB_true(k,:)')';

            tau_masaY(k,:) = ...
                cross(DeltaMr_Y,gB_true(k,:)')';


        elseif strcmpi(cfg.modo,'excitacion_2')

            tau_RW(k,:) = ...
                A.torque_RW2(t(k)-cfg.t0)';

 % incremental gravitational torques produced by the
 % mass displacements relative to the
 % second-stage reference configuration.
 %
 % They are stored for diagnostics and plots. The second
 % estimator does NOT add them to tau_RW; it incorporates them through
 % Mr(t), DeltaJ(t), and h_rel(t).

            DeltaMr_X = ...
                P.mX*(dX(k)-cfg.dX_ref2)*P.eX;

            DeltaMr_Y = ...
                P.mY*(dY(k)-cfg.dY_ref2)*P.eY;

            tau_masaX(k,:) = ...
                cross(DeltaMr_X,gB_true(k,:)')';

            tau_masaY(k,:) = ...
                cross(DeltaMr_Y,gB_true(k,:)')';

        end

    end

 % Known incremental gravitational torque due to the motion
 % of the X/and masses. This is the only known torque appearing on the
 % right-hand side of the first 10P estimation.
    tau_masas = tau_masaX + tau_masaY;

 % tau_exc is retained as to diagnostic signal equivalent to the
 % previous version: wheel reaction torque + incremental torque
 % from the masses. It is not longer used as an input to the 10P estimator.
    tau_exc = tau_RW + tau_masas;

 % Explicit name of the known input to the 10P estimator.
    tau_est1 = tau_masas;


 %% ============================================================
 % PLOTS FOR THIS SIMULATION
 %
 % They are generated automatically after EACH call to
 % simulate_platform.m.
 %
 % figure 1:
 % ideal and noisy angular velocity
 %
 % figure 2:
 % ideal and noisy gravity vector
 % ============================================================

    nombre_etapa = strrep(cfg.modo,'_',' ');

 % ------------------------------------------------------------
 % angular VELOCITY
 % ------------------------------------------------------------

    figure( ...
        'Name',['Velocidad angular - ' nombre_etapa], ...
        'NumberTitle','off' ...
    );

    subplot(3,1,1);

    plot(t,wx_true,'LineWidth',1.3);
    hold on;
    plot(t,wx,'LineWidth',0.6);

    grid on;

    ylabel('\omega_x [rad/s]');

    title([ ...
        'Velocidad angular - ', ...
        nombre_etapa ...
    ]);

    legend('Ideal','Ruidosa','Location','best');


    subplot(3,1,2);

    plot(t,wy_true,'LineWidth',1.3);
    hold on;
    plot(t,wy,'LineWidth',0.6);

    grid on;

    ylabel('\omega_y [rad/s]');


    subplot(3,1,3);

    plot(t,wz_true,'LineWidth',1.3);
    hold on;
    plot(t,wz,'LineWidth',0.6);

    grid on;

    ylabel('\omega_z [rad/s]');
    xlabel('Tiempo [s]');


 % ------------------------------------------------------------
 % GRAVITY VECTOR
 % ------------------------------------------------------------

    figure( ...
        'Name',['Vector de gravedad - ' nombre_etapa], ...
        'NumberTitle','off' ...
    );

    subplot(3,1,1);

    plot(t,gx_true,'LineWidth',1.3);
    hold on;
    plot(t,gx,'LineWidth',0.6);

    grid on;

    ylabel('g_x [m/s^2]');

    title([ ...
        'Vector de gravedad en BODY - ', ...
        nombre_etapa ...
    ]);

    legend('Ideal','Ruidosa','Location','best');


    subplot(3,1,2);

    plot(t,gy_true,'LineWidth',1.3);
    hold on;
    plot(t,gy,'LineWidth',0.6);

    grid on;

    ylabel('g_y [m/s^2]');


    subplot(3,1,3);

    plot(t,gz_true,'LineWidth',1.3);
    hold on;
    plot(t,gz,'LineWidth',0.6);

    grid on;

    ylabel('g_z [m/s^2]');
    xlabel('Tiempo [s]');


 %% ============================================================
 % MOVING-MASS POSITIONS
 %
 % The relative displacements dX, dY, and dZ are shown
 % relative to their initial positions defined in
 % platform_parameters.m.
 % ============================================================

    figure( ...
        'Name',['Posicion de masas - ' nombre_etapa], ...
        'NumberTitle','off' ...
    );

    plot(t,1000*dX,'LineWidth',1.3);
    hold on;

    plot(t,1000*dY,'LineWidth',1.3);
    plot(t,1000*dZ,'LineWidth',1.3);

    grid on;

    xlabel('Tiempo [s]');
    ylabel('Desplazamiento [mm]');

    legend( ...
        'd_X', ...
        'd_Y', ...
        'd_Z', ...
        'Location','best' ...
    );

    title([ ...
        'Posicion de las masas moviles - ', ...
        nombre_etapa ...
    ]);


 %% ============================================================
 % OUTPUT
 % ============================================================

    datos = struct;

    datos.t = t;
    datos.x = x;

    datos.q = q;

    datos.phi = phi;
    datos.theta = theta;
    datos.psi = psi;

    datos.phi_deg = rad2deg(phi);
    datos.theta_deg = rad2deg(theta);
    datos.psi_deg = rad2deg(psi);

    datos.wx_true = wx_true;
    datos.wy_true = wy_true;
    datos.wz_true = wz_true;

    datos.wx = wx;
    datos.wy = wy;
    datos.wz = wz;

    datos.gx_true = gx_true;
    datos.gy_true = gy_true;
    datos.gz_true = gz_true;
    datos.gB_true = gB_true;

    datos.gx = gx;
    datos.gy = gy;
    datos.gz = gz;
    datos.gB = gB;

    datos.dX = dX;
    datos.dY = dY;
    datos.dZ = dZ;

    datos.vX = vX;
    datos.vY = vY;
    datos.vZ = vZ;

    datos.xiX = xiX;
    datos.xiY = xiY;

    datos.Omega_RW = Omega_RW;
    datos.Omega_RW_rpm = Omega_RW*(60/(2*pi));
    datos.h_RW = h_RW;

    datos.rCM_real = rCM;

    datos.tau_RW = tau_RW;
    datos.tau_masaX = tau_masaX;
    datos.tau_masaY = tau_masaY;
    datos.tau_masas = tau_masas;
    datos.tau_est1 = tau_est1;
    datos.tau_exc = tau_exc;

    datos.cfg = cfg;

end
