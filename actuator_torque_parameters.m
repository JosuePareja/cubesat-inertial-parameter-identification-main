function A = actuator_torque_parameters()
%ACTUATOR_TORQUE_PARAMETERS
% source unica of parameters of control, actuators, excitaciones,
%
% configuration of time REDUCIDO
% --------------------------------
% this version removes the espera libre initial: t = 0 corresponds to the
% instant in that starts the balancing.
%
% Values validated with the standalone balancing test:
% Ktheta = Kphi = 0.08 m/(s*rad)
% vmax_X = vmax_Y = 1 mm/s
% amax_X = amax_Y = 0.5 mm/s^2
%
% the timeline of the excitaciones compacta for reducir the time
% simulated time without changing the estimator structure.
%
% timing, noise, and file settings.
%
% the sequence complete es:
%
% 1) balancing initial X/and
% 2) first excitation + estimation 10P
% 3) rebalancing X/and
% 4) compensation of the CM in Z
% 5) Espera until masses stationary
% 6) second excitation: masses X/and fast + wheel RW-Z
% 7) second estimation of J and b with J(t) and masses moving

 %% ============================================================
 % control of balancing X/and
 % ============================================================

 % control PI angular. the output es to reference of velocity
 % for the masses X/and.
    A.Kp_bal = 0.08;
    A.Ki_bal = 0.002;

 % Alias conservados by compatibilidad with scripts antiguos.
    A.Ktheta = A.Kp_bal;
    A.Kphi   = A.Kp_bal;

    A.error_lim_xy = deg2rad(0.01);     % [rad]

 % discharges suave of the states integrales within of the band dead.
    A.tau_reset_int_bal = 2.0;          % [s]

 % start of the first balancing
    A.T_balanceo1 = 0;               % [s]


 %% ============================================================
 % parameters of the actuator of masses
 % ============================================================

    A.vmax_X = 0.001;                  % [m/s]
    A.vmax_Y = 0.001;                  % [m/s]
    A.vmax_Z = 0.002;                  % [m/s]

    A.amax_X = 5.0e-4;                  % [m/s^2]
    A.amax_Y = 5.0e-4;                  % [m/s^2]
    A.amax_Z = 2.0e-3;                  % [m/s^2]

    A.tau_mX = 1.0;                     % [s]
    A.tau_mY = 1.0;                     % [s]
    A.tau_mZ = 0.5;                     % [s]

 % Local position control
    A.Kpos_X = 0.05;                    % [1/s]
    A.Kpos_Y = 0.05;                    % [1/s]
    A.Kpos_Z = 0.10;                    % [1/s]

 % limits relativos with respect to to the positions initial
    A.dX_min = -0.160;                   % [m]
    A.dX_max = +0.160;                   % [m]

    A.dY_min = -0.160;                   % [m]
    A.dY_max = +0.160;                   % [m]

 % travel allowed of the mass Z in this experiment.
 % the target rCM,z = -46.9 mm requires nominalmente ~297 mm.
    A.dZ_min = -0.300;                   % [m]
    A.dZ_max = +1.000;                   % [m]


 %% ============================================================
 % INITIAL CONDITION of the wheels of reaction
 % ============================================================
 %
 % Omega_RW son velocities relative of the volantes with respect to to the
 % platform. to wheel can conservar Omega_RW != 0 although its torque
 % instantaneo sea cero.

    A.Omega_RW0 = [0;0;0];              % [rad/s], ruedas X/Y/Z


 %% ============================================================
 % first stage of excitation
 %
 % in this version the first estimation excita only through
 % the motion of the masses X/and. the wheel Z remains desactivada
 % during this stage (to.tauZ1 = 0).
 % ============================================================

    A.T_exc1 = 100;                    % [s]

 % wheel of reaction Z desactivada in the first estimation.
 % T_Z1 retains as reference temporal of the stage.
    A.T_Z1 = A.T_exc1;                  % [s]
    A.tauZ1 = 0.0;                      % [N*m]
    A.DT_Z1 = 0.3;                      % [s]

 % masses X/and for generate excitaciones gravitational
    A.T_X1 = A.T_exc1 + 20;             % [s]
    A.T_Y1 = A.T_X1 + 120;             % [s]

    A.deltaX1 = +0.015;                % [m]
    A.deltaY1 = +0.015;                % [m]

    A.T_fin_exc1 = A.T_Y1 + 120;               % [s]


 %% ============================================================
 % rebalancing after of the first estimation
 % ============================================================

    A.T_rebalanceo_max = 150;          % [s]

 % Condicion numerica of "mass stationary"
    A.v_tol = 1e-6;                    % [m/s]

 % time that must mantenerse the condicion of balance
    A.t_hold_xy = 20;                  % [s]


 %% ============================================================
 % compensation of the center of mass in Z
 % ============================================================

 % Z-axis CM targets for the two experiment campaigns.
 %
 % BASE:
 % retains the value validado in the platform without CubeSat.
 %
 % BASE + CubeSat:
 % after of the first estimation 10P the mass Z desplaza
 % automatically until the global CM is left 5 mm below CR.
    A.rCMz_obj_base = -27.9e-3;         % [m]
    A.rCMz_obj_cubesat = -5.0e-3;       % [m]

 % value by defecto. run_identification_campaign receives to copia of to and the
 % script maestro replaces this campo according to the campaign.
    A.rCMz_obj = A.rCMz_obj_base;

 % Safety/acceptance criterion for starting the second test.
    A.rCMz_lim = 51.0e-3;               % [m]

    A.pos_tol_Z = 2.0e-4;              % [m]
    A.t_hold_Z = 20;                    % [s]

    A.T_compZ_max = 500;              % [s]


 %% ============================================================
 % ESPERA after of COMPENSAR Z
 % ============================================================

    A.T_espera = 30;                  % [s]


 %% ============================================================
 % second stage of excitation
 %
 % target:
 % - generate omega_x and omega_y of the order of 0.1 rad/s usando
 % movimientos fast of the masses and and X, respectively;
 % - retain the Z reaction wheel, which already works well for
 % excitar omega_z.
 %
 % the mass and excita principalmente the axis X.
 % the mass X excita principalmente the axis and.
 %
 % Times are RELATIVE to the start of the second test.
 % ============================================================

    A.omega_obj_exc2 = 0.10;            % [rad/s], objetivo aproximado

 % ------------------------------------------------------------
 % wheel of reaction Z
 % ------------------------------------------------------------

    A.tZ2 = 20;                        % [s]
    A.tauZ2 = 4e-3;                    % [N*m] = 5 mNm
    A.DT_Z2 = 5.0;                     % [s]

 % second pulse Z, same magnitud and signo contrario.
    A.tZ2_opuesto = 100;       % [s]

 % ------------------------------------------------------------
 % excitation of omega_x through the mass and
 %
 % profile:
 % dY_ref -> dY_ref + deltaY2 -> dY_ref
 % dY_ref -> dY_ref - deltaY2 -> dY_ref
 % ------------------------------------------------------------

    A.deltaY2 = 0.05;                 % [m] = 7 mm

    A.tY2_pos_ini = 180;               % [s]
    A.tY2_pos_fin = 205;   % [s]

    A.tY2_neg_ini = 280;               % [s]
    A.tY2_neg_fin = 305;   % [s]

 % ------------------------------------------------------------
 % excitation of omega_y through the mass X
 % ------------------------------------------------------------

    A.deltaX2 = 0.05;                 % [m] = 7 mm

    A.tX2_pos_ini = 380;              % [s]
    A.tX2_pos_fin = 405;   % [s]

    A.tX2_neg_ini = 480;              % [s]
    A.tX2_neg_fin = 505;   % [s]

 % ------------------------------------------------------------
 % dynamics fast of the masses X/and during the second stage.
 %
 % these values not modifican the balancing ni the first estimation.
 % only use when cfg.modo = 'excitacion_2'.
 %
 % with delta = 15 mm and rCM,z ~ -47 mm espera to respuesta of the
 % order of 0.1 rad/s. The true peak is printed at the end of the test.
 % ------------------------------------------------------------

    A.vmax_X2 = 0.05;                  % [m/s] = 50 mm/s
    A.vmax_Y2 = 0.05;                  % [m/s] = 50 mm/s

    A.amax_X2 = 0.50;                   % [m/s^2]
    A.amax_Y2 = 0.50;                   % [m/s^2]

    A.tau_mX2 = 0.03;                   % [s]
    A.tau_mY2 = 0.03;                   % [s]

 % ------------------------------------------------------------
 % control PI of position - second excitation
 % ------------------------------------------------------------
 % the output of the PI es the reference of velocity of each mass:
 %
 % v_cmd = Kp*e + Ki*integral(e)
 %
 % only applies to X/and during cfg.modo = 'excitacion_2'.

    A.Kp_X2 = 10.0;                     % [1/s]
    A.Kp_Y2 = 10.0;                     % [1/s]

    A.Ki_X2 = 5.0;                      % [1/s^2]
    A.Ki_Y2 = 5.0;                      % [1/s^2]

 % tolerance of position realista.
 % if |d_ref-d| <= 0.1 mm, the command of the motor hace cero.
    A.tol_pos_X2 = 0.1e-3;              % [m] = 0.1 mm
    A.tol_pos_Y2 = 0.1e-3;              % [m] = 0.1 mm

 % within of the band dead discharges the state integral of way
 % suave for evitar memory of the pulse previous.
    A.tau_reset_int_X2 = 0.20;           % [s]
    A.tau_reset_int_Y2 = 0.20;           % [s]

    A.T_exc2_total = 600;              % [s]


 %% ============================================================
 % INTEGRATION and noise
 % ============================================================

    A.dt = 0.02;                       % [s]
    A.fs = 1/A.dt;                     % [Hz]

    A.sigma_gyro = 0.0003;              % [rad/s]
    A.sigma_gyro_deg_s = rad2deg(A.sigma_gyro);          % [deg/s]

    A.sigma_g = 0.0003;                 % [m/s^2]

    A.rng_seed = 1;                     % reproducible noise


 %% ============================================================
 % FILTRO SAVITZKY-GOLAY - second estimation
 %
 % applies after of the second simulation and before of
 % estimate_inertia_parameters.m.
 %
 % only filtran the signals measurements:
 % wx, wy, wz
 % gx, gy, gz
 %
 % not filtran:
 % - the signals ideales *_true
 % - the torques known of the wheels
 %
 % with dt = 0.02 s, to window of 25 samples corresponds
 % approximately to 0.5 s.
 % ============================================================

    A.SG_orden = 3;
    A.SG_ventana = 25;                  % [muestras] = 0.50 s con dt=0.02
    A.SG_filtrar_gravedad = true;


 %% ============================================================
 % estimators
 % ============================================================

    A.stride_est1 = 1;
    A.stride_est2 = 1;

 % ------------------------------------------------------------
 % DATA used in the second estimation
 %
 % Opciones:
 %
 % 'noisy' -> uses directamente wx, wy, wz, gx, gy, gz
 % with noise, without Savitzky-Golay.
 %
 % 'filtered' -> uses esas same measurements after of apply
 % the Savitzky-Golay filter.
 %
 % the signals ideales siguen utilizando only as reference
 % for compare the desempeno of the estimator.
 % ------------------------------------------------------------

    A.modo_datos_est2 = 'ruidosa';     % 'ruidosa' o 'filtrada'


 %% ============================================================
 % files
 % ============================================================

    A.archivo_etapa1 = 'stage1_data.mat';
    A.archivo_resultado1 = '10P_estimation_result.mat';

    A.archivo_etapa2 = 'stage2_data.mat';
    A.archivo_etapa2_filtrado = 'stage2_data_SG_filtered.mat';
    A.archivo_resultado2 = 'inertia_estimation_result.mat';

    A.archivo_secuencia = 'complete_sequence_result.mat';


 %% ============================================================
 % REACTION-TORQUE FUNCTIONS ABOUT THE PLATFORM
 % ============================================================
 %
 % Convencion of signo:
 % tau_RW es the torque that the wheel ejerce about the platform.
 % by accion-reaction, the torque about the volante es -tau_RW and:
 %
 % J_RW .* Omega_RW_dot = -tau_RW
 %
 % in the dynamics rotacional tau_RW already not suma as torque externo;
 % aparece to traves of -d(h_RW)/dt and of the termino giroscopico h_RW.
 % ============================================================

    TZ1 = A.T_Z1;
    DTZ1 = A.DT_Z1;
    tauZ1 = A.tauZ1;

    A.torque_RW1 = @(t) torque_rw1(t,TZ1,DTZ1,tauZ1);

    tZ2 = A.tZ2;
    tZ2_opuesto = A.tZ2_opuesto;
    DTZ2 = A.DT_Z2;
    tauZ2 = A.tauZ2;

    A.torque_RW2 = @(trel) ...
        torque_rw2(trel,tZ2,tZ2_opuesto,DTZ2,tauZ2);


 %% ============================================================
 % CHECKS
 % ============================================================

    if A.T_balanceo1 >= A.T_exc1
        error('T_balanceo1 must be smaller than T_exc1.');
    end

    if A.T_fin_exc1 <= A.T_Y1
        error('T_fin_exc1 must be later than T_Y1.');
    end

    if abs(A.rCMz_obj) >= A.rCMz_lim
        error('rCMz_obj must satisfy |rCMz_obj| < rCMz_lim.');
    end

    if mod(A.SG_ventana,2) == 0
        error('A.SG_ventana must be an odd number.');
    end

    if A.SG_ventana <= A.SG_orden
        error('A.SG_ventana must be greater than A.SG_orden.');
    end

    if A.tY2_pos_fin <= A.tY2_pos_ini || A.tY2_neg_fin <= A.tY2_neg_ini
        error('the intervalos of excitation of the mass and not son validos.');
    end

    if A.tX2_pos_fin <= A.tX2_pos_ini || A.tX2_neg_fin <= A.tX2_neg_ini
        error('the intervalos of excitation of the mass X not son validos.');
    end

    if max([A.tY2_neg_fin,A.tX2_neg_fin,A.tZ2_opuesto+A.DT_Z2]) ...
            >= A.T_exc2_total
        error('any excitation of the second stage remains outside of T_exc2_total.');
    end

    modo_est2 = lower(strtrim(A.modo_datos_est2));

    if ~ismember(modo_est2,{'ruidosa','filtrada'})
        error([ ...
            'A.modo_datos_est2 debe ser ''ruidosa'' o ''filtrada''.' ...
        ]);
    end

    if A.Kp_X2 < 0 || A.Kp_Y2 < 0 || A.Ki_X2 < 0 || A.Ki_Y2 < 0
        error('The PI Kp/Ki gains must be nonnegative.');
    end

    if A.tol_pos_X2 <= 0 || A.tol_pos_Y2 <= 0
        error('The X/Y position tolerances must be positive.');
    end

    if A.tau_reset_int_X2 <= 0 || A.tau_reset_int_Y2 <= 0
        error('The integrator discharge time constants must be positive.');
    end

end


function tau = torque_rw1(t,TZ,DTZ,tauZ)

    tau = zeros(3,1);

    if t >= TZ && t < TZ+DTZ
        tau(3) = tauZ;
    end

end


function tau = torque_rw2(t,tZ,tZ_opuesto,DTZ,tauZ)

    tau = zeros(3,1);

 % during the second stage only actua the wheel Z.
    if t >= tZ && t < tZ+DTZ
        tau(3) = tauZ;

    elseif t >= tZ_opuesto && t < tZ_opuesto+DTZ
        tau(3) = -tauZ;
    end

end
