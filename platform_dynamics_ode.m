function dx = platform_dynamics_ode(t,x,cfg,P,A)
%PLATFORM_DYNAMICS_ODE
% dynamics complete of the platform with masses moving.
%
% states:
%
% x = [qx qy qz qw ...
% wx wy wz ...
% dX dY dZ ...
% vX vY vZ ...
% xiX xiY ...
% OmegaRW_X OmegaRW_Y OmegaRW_Z]'
%
% xiX, xiY son states integrales reutilizados by stage:
% - balancing/rebalancing/compensation: PI angular X/and;
% - second excitation: PI of position X/and.
% run_identification_campaign resets explicitly these states to the change of
% control function to avoid memory between stages.
%
% OmegaRW_X/and/Z son the velocities angulares relative of the volantes
% with respect to to the platform. the wheels can keep OmegaRW != 0 still
% when their instantaneous torque is zero.
%
% Modos available:
%
% balanceo_inicial
% excitacion_1
% rebalanceo_xy
% compensacion_z
% espera
% excitacion_2

    if nargin < 4 || isempty(P)
        P = platform_parameters();
    end

    if nargin < 5 || isempty(A)
        A = actuator_torque_parameters();
    end


 %% ============================================================
 % STATES
 % ============================================================

    q = x(1:4);
    w = x(5:7);

    dX = x(8);
    dY = x(9);
    dZ = x(10);

    vX = x(11);
    vY = x(12);
    vZ = x(13);

 % states integrales of the controllers PI X/and.
    xiX = x(14);
    xiY = x(15);

 % Relative speeds of the three reaction wheels.
    Omega_RW = x(16:18);


 %% ============================================================
 % QUATERNION AND ANGLES
 % ============================================================

    q = q/norm(q);

    qx = q(1);
    qy = q(2);
    qz = q(3);
    qw = q(4);

    sin_theta = 2*(qw*qy-qz*qx);
    sin_theta = max(-1,min(1,sin_theta));

    phi = atan2( ...
        2*(qw*qx+qy*qz), ...
        1-2*(qx^2+qy^2) ...
    );

    theta = asin(sin_theta);


 %% ============================================================
 % MASS COMMANDS
 % ============================================================

    vX_cmd = 0;
    vY_cmd = 0;
    vZ_cmd = 0;

 % by defecto the integrators remain frozen.
    xiX_dot = 0;
    xiY_dot = 0;


    switch lower(cfg.modo)

        case 'balanceo_inicial'

            if t >= A.T_balanceo1

                [vX_cmd,vY_cmd,xiX_dot,xiY_dot] = ...
                    control_balanceo_xy_PI(phi,theta,xiX,xiY,A);

            end


        case 'excitacion_1'

            dX_target = cfg.dX_bal;
            dY_target = cfg.dY_bal;
            dZ_target = cfg.dZ_bal;

            if t >= A.T_X1
                dX_target = cfg.dX_bal + A.deltaX1;
            end

            if t >= A.T_Y1
                dY_target = cfg.dY_bal + A.deltaY1;
            end

            vX_cmd = A.Kpos_X*(dX_target-dX);
            vY_cmd = A.Kpos_Y*(dY_target-dY);
            vZ_cmd = A.Kpos_Z*(dZ_target-dZ);


        case 'rebalanceo_xy'

            [vX_cmd,vY_cmd,xiX_dot,xiY_dot] = ...
                control_balanceo_xy_PI(phi,theta,xiX,xiY,A);

            vZ_cmd = A.Kpos_Z*(cfg.dZ_hold-dZ);


        case 'compensacion_z'

 % keeps the balancing X/and while Z desplaza.
            [vX_cmd,vY_cmd,xiX_dot,xiY_dot] = ...
                control_balanceo_xy_PI(phi,theta,xiX,xiY,A);

            vZ_cmd = A.Kpos_Z*(cfg.dZ_target-dZ);


        case 'espera'

            vX_cmd = A.Kpos_X*(cfg.dX_hold-dX);
            vY_cmd = A.Kpos_Y*(cfg.dY_hold-dY);
            vZ_cmd = A.Kpos_Z*(cfg.dZ_hold-dZ);


        case 'excitacion_2'

 % second identification:
 % - mass and excita principalmente omega_x;
 % - mass X excita principalmente omega_y;
 % - Z mass remains at its compensation position;
 % - the Z reaction wheel generates omega_z excitation.
 %
 % the movimientos X/and son fast and bipolar with respect to to
 % the reference configuration at the start of this stage.

            trel = t-cfg.t0;

            dX_target = objetivo_masa_exc2( ...
                trel,cfg.dX_ref2,A.deltaX2, ...
                A.tX2_pos_ini,A.tX2_pos_fin, ...
                A.tX2_neg_ini,A.tX2_neg_fin ...
            );

            dY_target = objetivo_masa_exc2( ...
                trel,cfg.dY_ref2,A.deltaY2, ...
                A.tY2_pos_ini,A.tY2_pos_fin, ...
                A.tY2_neg_ini,A.tY2_neg_fin ...
            );

            dX_target = max(A.dX_min,min(A.dX_max,dX_target));
            dY_target = max(A.dY_min,min(A.dY_max,dY_target));

 % ----------------------------------------------------
 % control PI of position X/and with band dead
 % ----------------------------------------------------
 %
 % the PI entrega to reference of velocity for each mass.
 % When |position error| <= 0.1 mm (configurable value),
 % the command of the motor hace cero. the integrator deja of
 % acumular error and discharges suavemente for evitar to golpe
 % when cambie the next reference.
 %
 % Anti-windup is also applied through conditional integration
 % when the reference of velocity reaches vmax_X2/Y2.

            [vX_cmd,xiX_dot] = control_posicion_PI_deadband( ...
                dX_target,dX,xiX, ...
                A.Kp_X2,A.Ki_X2,A.vmax_X2, ...
                A.tol_pos_X2,A.tau_reset_int_X2 ...
            );

            [vY_cmd,xiY_dot] = control_posicion_PI_deadband( ...
                dY_target,dY,xiY, ...
                A.Kp_Y2,A.Ki_Y2,A.vmax_Y2, ...
                A.tol_pos_Y2,A.tau_reset_int_Y2 ...
            );

 % Z remains estacionaria in the position alcanzada during
 % center-of-mass compensation.
            vZ_cmd = A.Kpos_Z*(cfg.dZ_ref2-dZ);


        otherwise

            error('Modo desconocido: %s',cfg.modo);

    end


 %% ============================================================
 % parameters of the actuator ACTIVOS in this stage
 % ============================================================

    if strcmpi(cfg.modo,'excitacion_2')

        vmax_X_use = A.vmax_X2;
        vmax_Y_use = A.vmax_Y2;

        amax_X_use = A.amax_X2;
        amax_Y_use = A.amax_Y2;

        tau_mX_use = A.tau_mX2;
        tau_mY_use = A.tau_mY2;

    else

        vmax_X_use = A.vmax_X;
        vmax_Y_use = A.vmax_Y;

        amax_X_use = A.amax_X;
        amax_Y_use = A.amax_Y;

        tau_mX_use = A.tau_mX;
        tau_mY_use = A.tau_mY;

    end

 % The Z mass retains its nominal dynamics.
    vmax_Z_use = A.vmax_Z;
    amax_Z_use = A.amax_Z;
    tau_mZ_use = A.tau_mZ;


 %% ============================================================
 % SATURACIONES of velocity
 % ============================================================

    vX_cmd = max(-vmax_X_use,min(vmax_X_use,vX_cmd));
    vY_cmd = max(-vmax_Y_use,min(vmax_Y_use,vY_cmd));
    vZ_cmd = max(-vmax_Z_use,min(vmax_Z_use,vZ_cmd));


 %% ============================================================
 % limits of travel
 % ============================================================

    if dX <= A.dX_min && vX_cmd < 0
        vX_cmd = 0;
    end

    if dX >= A.dX_max && vX_cmd > 0
        vX_cmd = 0;
    end

    if dY <= A.dY_min && vY_cmd < 0
        vY_cmd = 0;
    end

    if dY >= A.dY_max && vY_cmd > 0
        vY_cmd = 0;
    end

    if dZ <= A.dZ_min && vZ_cmd < 0
        vZ_cmd = 0;
    end

    if dZ >= A.dZ_max && vZ_cmd > 0
        vZ_cmd = 0;
    end


 %% ============================================================
 % dynamics of the actuators
 % ============================================================

    aX = (vX_cmd-vX)/tau_mX_use;
    aY = (vY_cmd-vY)/tau_mY_use;
    aZ = (vZ_cmd-vZ)/tau_mZ_use;

    aX = max(-amax_X_use,min(amax_X_use,aX));
    aY = max(-amax_Y_use,min(amax_Y_use,aY));
    aZ = max(-amax_Z_use,min(amax_Z_use,aZ));


 %% ============================================================
 % CINEMATICA RELATIVA
 % ============================================================

    dX_dot = vX;
    dY_dot = vY;
    dZ_dot = vZ;

    r_X = P.rX0 + P.eX*dX;
    r_Y = P.rY0 + P.eY*dY;
    r_Z = P.rZ0 + P.eZ*dZ;

    rX_dot = P.eX*vX;
    rY_dot = P.eY*vY;
    rZ_dot = P.eZ*vZ;

    rX_ddot = P.eX*aX;
    rY_ddot = P.eY*aY;
    rZ_ddot = P.eZ*aZ;


 %% ============================================================
 % first moment TOTAL
 % ============================================================

    Mr_total = ...
        P.M_p*P.r_P ...
        + P.mX*r_X ...
        + P.mY*r_Y ...
        + P.mZ*r_Z;


 %% ============================================================
 % MATRIZ of inertia TOTAL
 % ============================================================

    I3 = eye(3);

    J_X = P.J_X_CM + ...
        P.mX*((r_X'*r_X)*I3-r_X*r_X');

    J_Y = P.J_Y_CM + ...
        P.mY*((r_Y'*r_Y)*I3-r_Y*r_Y');

    J_Z = P.J_Z_CM + ...
        P.mZ*((r_Z'*r_Z)*I3-r_Z*r_Z');

    J = P.J_P + J_X + J_Y + J_Z;


 %% ============================================================
 % J_dot
 % ============================================================

    Jdot_X = P.mX*( ...
        2*(r_X'*rX_dot)*I3 ...
        - rX_dot*r_X' ...
        - r_X*rX_dot' ...
    );

    Jdot_Y = P.mY*( ...
        2*(r_Y'*rY_dot)*I3 ...
        - rY_dot*r_Y' ...
        - r_Y*rY_dot' ...
    );

    Jdot_Z = P.mZ*( ...
        2*(r_Z'*rZ_dot)*I3 ...
        - rZ_dot*r_Z' ...
        - r_Z*rZ_dot' ...
    );

    J_dot = Jdot_X + Jdot_Y + Jdot_Z;


 %% ============================================================
 % moment angular relative
 % ============================================================

    h_rel = ...
        P.mX*cross(r_X,rX_dot) ...
        + P.mY*cross(r_Y,rY_dot) ...
        + P.mZ*cross(r_Z,rZ_dot);

    hdot_rel = ...
        P.mX*cross(r_X,rX_ddot) ...
        + P.mY*cross(r_Y,rY_ddot) ...
        + P.mZ*cross(r_Z,rZ_ddot);


 %% ============================================================
 % gravity in BODY
 % ============================================================

    gx = -P.g*2*(qx*qz-qy*qw);
    gy = -P.g*2*(qy*qz+qx*qw);
    gz = -P.g*(1-2*(qx^2+qy^2));

    g_B = [gx;gy;gz];


 %% ============================================================
 % torque GRAVITACIONAL
 % ============================================================

    tau_g = cross(Mr_total,g_B);


 %% ============================================================
 % wheels of reaction: torque, velocity and moment angular
 % ============================================================
 %
 % tau_RW es the torque of reaction ejercido about the platform.
 % the torque about each volante tiene signo contrario. Therefore:
 %
 % J_RW .* Omega_RW_dot = -tau_RW
 %
 % and the moment angular relative almacenado in the wheels es:
 %
 % h_RW = J_RW .* Omega_RW
 %
 % the wheel not fuerza to detenerse when tau_RW = 0. in ese caso
 % Omega_RW remains constant and h_RW remains part of the total
 % angular momentum of the system.

    tau_RW = zeros(3,1);

    if strcmpi(cfg.modo,'excitacion_1')

        tau_RW = A.torque_RW1(t);

    elseif strcmpi(cfg.modo,'excitacion_2')

        trel = t-cfg.t0;
        tau_RW = A.torque_RW2(trel);

    end

    Omega_RW_dot = -tau_RW./P.J_RW;

    h_RW = P.J_RW.*Omega_RW;
    hdot_RW = P.J_RW.*Omega_RW_dot;


 %% ============================================================
 % dynamics ROTACIONAL
 % ============================================================
 %
 % moment angular total:
 %
 % H = J*w + h_rel + h_RW
 %
 % the wheels son actuators INTERNOS. by that tau_RW not aparece as
 % torque externo independiente: its efecto this incluido in -hdot_RW.

    w_dot = J \ ( ...
        tau_g ...
        - P.b*w ...
        - J_dot*w ...
        - hdot_rel ...
        - hdot_RW ...
        - cross(w,J*w+h_rel+h_RW) ...
    );


 %% ============================================================
 % CINEMATICA of CUATERNION
 % ============================================================

    qv = [qx;qy;qz];

    S_qv = [ ...
         0,   -qz,  qy; ...
         qz,   0,  -qx; ...
        -qy,   qx,  0 ...
    ];

    qv_dot = 0.5*(qw*eye(3)+S_qv)*w;

    qw_dot = -0.5*qv'*w;


 %% ============================================================
 % OUTPUT
 % ============================================================

    dx = [ ...
        qv_dot; ...
        qw_dot; ...
        w_dot; ...
        dX_dot; ...
        dY_dot; ...
        dZ_dot; ...
        aX; ...
        aY; ...
        aZ; ...
        xiX_dot; ...
        xiY_dot; ...
        Omega_RW_dot ...
    ];

end


function [v_cmd,xi_dot] = control_posicion_PI_deadband( ...
    d_target,d,xi,Kp,Ki,vmax,tol_pos,tau_reset)
%CONTROL_POSICION_PI_DEADBAND
% PI of position cuya output es to reference of velocity.
%
% e = d_target-d
% v_unsat = Kp*e + Ki*xi
%
% Caracteristicas:
% 1) position deadband: if |e| <= tol_pos, v_cmd = 0;
% 2) within the deadband, the integral state is smoothly discharged;
% 3) anti-windup by integration conditional when v_cmd satura.
%
% the mass not detiene of way instantanea to the entrar in the band.
% Simplemente ordena v_cmd = 0 and the dynamics of the actuator determina
% the deceleration through tau_m and amax.

    e = d_target-d;

 % ------------------------------------------------------------
 % motor / control off within of the tolerance
 % ------------------------------------------------------------
    if abs(e) <= tol_pos

        v_cmd = 0;

 % Smoothly discharges the integrator. This avoids retaining a large memory
 % of the previous pulse when d_target changes later.
        xi_dot = -xi/tau_reset;
        return;

    end

 % ------------------------------------------------------------
 % PI outside of the band dead
 % ------------------------------------------------------------
    v_unsat = Kp*e + Ki*xi;
    v_cmd = max(-vmax,min(vmax,v_unsat));

 % ------------------------------------------------------------
 % ANTI-WINDUP by INTEGRATION conditional
 % ------------------------------------------------------------
 % if not hay saturation, integramos normalmente.
 % if hay saturation and the error empuja still more toward the saturation,
 % congelamos the integral. if the error ayuda to salir of saturation,
 % permitimos integrate for descargarlo.

    if abs(v_unsat) <= vmax

        xi_dot = e;

    elseif (v_unsat > vmax && e < 0) || ...
           (v_unsat < -vmax && e > 0)

        xi_dot = e;

    else

        xi_dot = 0;

    end

end


function d_target = objetivo_masa_exc2( ...
    trel,d_ref,delta,t_pos_ini,t_pos_fin,t_neg_ini,t_neg_fin)
%OBJETIVO_MASA_EXC2
% profile bipolar with respect to to to position of reference:
%
% reference -> +delta -> reference -> -delta -> reference

    d_target = d_ref;

    if trel >= t_pos_ini && trel < t_pos_fin

        d_target = d_ref + delta;

    elseif trel >= t_neg_ini && trel < t_neg_fin

        d_target = d_ref - delta;

    end

end


function [vX_cmd,vY_cmd,xiX_dot,xiY_dot] = ...
    control_balanceo_xy_PI(phi,theta,xiX,xiY,A)
%CONTROL_BALANCEO_XY_PI
% Angular PI controller with velocity-reference output.
%
% the mass X corrects pitch and the mass and corrects roll:
% eX = -theta
% eY = phi
%
% v_cmd = Kp_bal*e + Ki_bal*xi
%
% includes band dead angular, discharges suave of the integrator and
% anti-windup by integration conditional.

    eX = -theta;
    eY =  phi;

    [vX_cmd,xiX_dot] = control_PI_deadband_generico( ...
        eX,xiX,A.Kp_bal,A.Ki_bal,A.vmax_X, ...
        A.error_lim_xy,A.tau_reset_int_bal);

    [vY_cmd,xiY_dot] = control_PI_deadband_generico( ...
        eY,xiY,A.Kp_bal,A.Ki_bal,A.vmax_Y, ...
        A.error_lim_xy,A.tau_reset_int_bal);

end


function [u,xi_dot] = control_PI_deadband_generico( ...
    e,xi,Kp,Ki,umax,deadband,tau_reset)
%CONTROL_PI_DEADBAND_GENERICO PI with saturation and anti-windup.

    if abs(e) <= deadband
        u = 0;
        xi_dot = -xi/tau_reset;
        return;
    end

    u_unsat = Kp*e + Ki*xi;
    u = max(-umax,min(umax,u_unsat));

    if abs(u_unsat) <= umax
        xi_dot = e;
    elseif (u_unsat > umax && e < 0) || ...
           (u_unsat < -umax && e > 0)
        xi_dot = e;
    else
        xi_dot = 0;
    end

end
