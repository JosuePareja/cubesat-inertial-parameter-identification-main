function P = platform_parameters()
%PLATFORM_PARAMETERS
% Single source of physical plant parameters.
%
% Convention:
% - positions with respect to the center of rotation (CR)
% - BODY coordinates
% - scalar-last quaternion
%
% IMPORTANT:
% P.J_P corresponds ONLY to the platform, without moving masses.
% P.J0 corresponds to the platform + moving masses at their
% initial positions.

 %% ============================================================
 % PHYSICAL CONSTANTS
 % ============================================================

    P.g = 9.81;                         % [m/s^2]
    P.b = 0.00242;                      % [N*m*s/rad]

 % Metadata. The add_cubesat_to_platform.m function changes these
 % fields when to rigid CubeSat is mounted on the platform.
    P.tiene_cubesat = false;
    P.nombre_cubesat = 'NINGUNO';


 %% ============================================================
 % REACTION WHEELS - ROTOR AXIAL INERTIA
 % ============================================================
 %
 % Each wheel is aligned with one BODY axis of the platform.
 % P.J_RW(i) is the rotor axial moment of inertia of the wheel
 % corresponding to X, and, and Z.
 %
 % IMPORTANT:
 % - this parameter is used for the RELATIVE angular momentum
 % h_RW = J_RW .* Omega_RW
 % - P.J_P must continue to represent the complete rigid inertia of
 % the platform about CR. If the wheel hardware was already
 % included in P.J_P, P.J_RW must NOT be added to J again.
 %
 % Current value: same rotor on all three axes. Change here if the
 % current wheels have different inertias.

% P.J_RW = [ ...
% 8.14e-6; ...
% 8.14e-6; ...
% 8.14e-6 ...
% ]; % [kg*m^2]

    P.J_RW = [ ...
        9.5493e-5; ...
        9.5493e-5; ...
        9.5493e-5 ...
    ];                                  % [kg*m^2]
 %% ============================================================
 % RIGID PLATFORM WITHOUT MOVING MASSES OR CUBESAT
 % ============================================================
 % P.M_p, P.r_P, and P.J_P ALREADY physically include the three complete
 % 375 g reaction wheels at their current locations.
 % P.J_RW is NOT added again to P.J_P; it is only used for the
 % relative angular momentum h_RW = J_RW .* Omega_RW.

    P.M_p = 3.667375;                   % [kg]

 % CR -> CM vector of the rigid platform
    P.r_P = [ ...
         0.001320; ...
         0.002952; ...
        -0.040999 ...
    ];                                  % [m]

 % Inertia matrix of the rigid platform about CR
    P.J_P = [ ...
         0.062580,  0.007738, -0.002098; ...
         0.007738,  0.066110,  0.001621; ...
        -0.002098,  0.001621,  0.085380 ...
    ];                                  % [kg*m^2]


 %% ============================================================
 % MOVING MASSES
 % ============================================================

    P.mX = 0.500;                       % [kg]
    P.mY = 0.500;                       % [kg]
    P.mZ = 0.500;                       % [kg]

    P.R_X = 0.025;                      % [m]
    P.R_Y = 0.025;                      % [m]
    P.R_Z = 0.025;                      % [m]

    P.h_X = 0.032439;                   % [m]
    P.h_Y = 0.032439;                   % [m]
    P.h_Z = 0.032439;                   % [m]

 % Allowed travels used in this experiment.
 % X/and: +/-160 mm relative to the initial position.
 % Z : up to approximately 300 mm to raise the CM.
    P.dX_max_abs = 0.160;               % [m]
    P.dY_max_abs = 0.160;               % [m]
    P.dZ_max_abs = 1.000;               % [m]


 %% ============================================================
 % CENTROIDAL INERTIAS OF THE MOVING MASSES
 % ============================================================

    Ipar_X  = 0.5*P.mX*P.R_X^2;
    Iperp_X = (1/12)*P.mX*(3*P.R_X^2 + P.h_X^2);

    Ipar_Y  = 0.5*P.mY*P.R_Y^2;
    Iperp_Y = (1/12)*P.mY*(3*P.R_Y^2 + P.h_Y^2);

    Ipar_Z  = 0.5*P.mZ*P.R_Z^2;
    Iperp_Z = (1/12)*P.mZ*(3*P.R_Z^2 + P.h_Z^2);

 % X mass: cylinder axis parallel to X
    P.J_X_CM = diag([Ipar_X Iperp_X Iperp_X]);

 % and mass: cylinder axis parallel to and
    P.J_Y_CM = diag([Iperp_Y Ipar_Y Iperp_Y]);

 % Z mass: cylinder axis parallel to Z
    P.J_Z_CM = diag([Iperp_Z Iperp_Z Ipar_Z]);


 %% ============================================================
 % INITIAL POSITIONS OF THE MOVING MASSES
 % ============================================================

    P.rX0 = [ ...
         0.077754; ...
        -0.125000; ...
        -0.170000 ...
    ];                                  % [m]

    P.rY0 = [ ...
         0.125000; ...
        -0.077754; ...
        -0.170000 ...
    ];                                  % [m]

    P.rZ0 = [ ...
        -0.169421; ...
         0.170177; ...
        -0.141246 ...
    ];                                  % [m]


 %% ============================================================
 % MOTION AXES
 % ============================================================

    P.eX = [1;0;0];
    P.eY = [0;1;0];
    P.eZ = [0;0;1];


 %% ============================================================
 % PARAMETERS DERIVED FROM THE original CONFIGURATION
 % ============================================================

    P.M_total = P.M_p + P.mX + P.mY + P.mZ;

    P.Mr0 = ...
        P.M_p*P.r_P ...
        + P.mX*P.rX0 ...
        + P.mY*P.rY0 ...
        + P.mZ*P.rZ0;

    P.rCM0 = P.Mr0/P.M_total;

    I3 = eye(3);

    JX0 = P.J_X_CM + ...
        P.mX*((P.rX0'*P.rX0)*I3 - P.rX0*P.rX0');

    JY0 = P.J_Y_CM + ...
        P.mY*((P.rY0'*P.rY0)*I3 - P.rY0*P.rY0');

    JZ0 = P.J_Z_CM + ...
        P.mZ*((P.rZ0'*P.rZ0)*I3 - P.rZ0*P.rZ0');

    P.J0 = P.J_P + JX0 + JY0 + JZ0;


 %% ============================================================
 % CHECKS
 % ============================================================

    if numel(P.J_RW) ~= 3 || any(P.J_RW <= 0)
        error('P.J_RW must contain three positive axial inertias.');
    end

    P.J_RW = P.J_RW(:);

end
