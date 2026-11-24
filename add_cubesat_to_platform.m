function P = add_cubesat_to_platform(P_base,C)
%ADD_CUBESAT_TO_PLATFORM
% Returns to new P structure representing the original rigid platform
% plus the rigidly mounted CubeSat, while leaving the three
% moving masses unchanged.
%
% The CubeSat is incorporated into the rigid set P.M_p/P.r_P/P.J_P.
% This way, platform_dynamics_ode.m, simulate_platform.m, and the estimators do not
% require to different dynamic formulation for the BASE+CubeSat campaign.

    if nargin < 2
        error('P_base and C are required.');
    end

    P = P_base;

    M_rigido_base = P_base.M_p;
    S_rigido_base = M_rigido_base*P_base.r_P;
    J_rigido_base = P_base.J_P;

    S_cubesat = C.M*C.r_CR_CM_B;
    J_cubesat_CR = C.J_CR_B;

    P.M_p = M_rigido_base + C.M;
    P.r_P = (S_rigido_base + S_cubesat)/P.M_p;
    P.J_P = J_rigido_base + J_cubesat_CR;

    P.tiene_cubesat = true;
    P.nombre_cubesat = C.nombre;
    P.cubesat = C;

 % =============================================================
 % RECONSTRUCTION OF DERIVED PARAMETERS AT dX=dY=dZ=0
 % =============================================================

    P.M_total = P.M_p + P.mX + P.mY + P.mZ;

    P.Mr0 = ...
        P.M_p*P.r_P ...
        + P.mX*P.rX0 ...
        + P.mY*P.rY0 ...
        + P.mZ*P.rZ0;

    P.rCM0 = P.Mr0/P.M_total;

    I3 = eye(3);

    JX0 = P.J_X_CM + ...
        P.mX*((P.rX0'*P.rX0)*I3-P.rX0*P.rX0');

    JY0 = P.J_Y_CM + ...
        P.mY*((P.rY0'*P.rY0)*I3-P.rY0*P.rY0');

    JZ0 = P.J_Z_CM + ...
        P.mZ*((P.rZ0'*P.rZ0)*I3-P.rZ0*P.rZ0');

    P.J0 = P.J_P + JX0 + JY0 + JZ0;

 % Check values used to verify the addition algebraically.
    P.cubesat.S_CR_B = S_cubesat;
    P.cubesat.J_CR_B = J_cubesat_CR;

end
