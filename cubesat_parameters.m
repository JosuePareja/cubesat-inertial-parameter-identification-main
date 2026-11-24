function C = cubesat_parameters(name)
%CUBESAT_PARAMETERS Reference mass properties for the validation spacecraft.
%
%   C = cubesat_parameters(NAME)
%
% Supported cases:
%   INNOCUBE, TPA1, SOURCE, M3SAT, INFLATESAIL
%
% Common spacecraft-to-platform mounting convention
% -------------------------------------------------
% The normalized spacecraft frame C is mounted on the platform BODY frame B as
%
%   +Y_C = +Y_B
%   +Z_C = +X_B
%   +X_C = -Z_B
%
% hence
%
%   [x_B y_B z_B]' = [z_C y_C -x_C]'
%   J_B = R_BC*J_C*R_BC'
%
% The platform rotation center CR is the BODY-frame origin.  The common
% geometric spacecraft mounting point is
%
%   r_CR->G_B = [0 0 62]' mm.
%
% The reference tensor used by the identification comparison is centroidal,
% J_CM, for every validation case.  SOURCE and M3Sat do not provide source
% products of inertia suitable for the final comparison; their reported
% comparison masks therefore retain only the three diagonal moments.  For the
% numerical plant, missing products are completed as described in the relevant
% case below.

    if nargin < 1 || isempty(name)
        error('The CubeSat name must be specified.');
    end

    R_BC = [ ...
         0, 0, 1; ...
         0, 1, 0; ...
        -1, 0, 0 ...
    ];

    key = upper(regexprep(char(string(name)),'[^A-Za-z0-9]',''));

    C = struct;
    C.R_BC = R_BC;
    C.r_CR_G_B = [0;0;62e-3];
    C.tensor_referencia_modelo = 'CM';
    C.tensor_referencia_fuente_confirmada = true;
    C.cm_is_assumed = false;

    switch key

        case {'INNOCUBE','INNO'}
            C.id = 'INNOCUBE';
            C.nombre = 'InnoCube';
            C.M = 4.20; % [kg]

            % Local CM offset after the source-frame-to-C-frame orientation
            % adopted for the common mounting configuration.
            C.r_G_CM_C = 1e-3*[ ...
                -6.83; ...
                -9.92; ...
                 5.66 ...
            ];

            % Centroidal inertia tensor in C axes [kg m^2].
            C.J_reportado_C = 1e-3*[ ...
                189,  38, -12; ...
                 38,  34, -39; ...
                -12, -39, 188 ...
            ];
            C.J_modelo_C = C.J_reportado_C;
            C.productos_inercia_reportados = true;
            C.inertia_metric_mode = 'full';
            C.data_note = 'Team-provided mass-property data; tensor referred to CM.';

        case {'TPA1','TPA'}
            C.id = 'TPA1';
            C.nombre = 'TPA-1';
            C.M = 4.50; % [kg]

            C.r_G_CM_C = 1e-3*[ ...
                 1.60; ...
                -2.60; ...
                 0.70 ...
            ];

            % The source spreadsheet reports the moments about the payload CG
            % and defines products with the positive-integral convention
            % Pxy = +int(x*y dm).  Tensor off-diagonal elements therefore carry
            % the opposite sign before the source-to-C frame rotation.  The
            % resulting centroidal tensor in C axes is below.
            C.J_reportado_C = [ ...
                0.042432754880694146,  0.00029284164859002165, -0.000029284164859002168; ...
                0.00029284164859002165, 0.005612798264642082,   0.0003709327548806942; ...
               -0.000029284164859002168,0.0003709327548806942,  0.04217895878524946 ...
            ];
            C.J_modelo_C = C.J_reportado_C;
            C.productos_inercia_reportados = true;
            C.inertia_metric_mode = 'full';
            C.data_note = 'Flight mass properties; inertia reported about payload CG.';

        case 'SOURCE'
            C.id = 'SOURCE';
            C.nombre = 'SOURCE';
            C.M = 5.02; % [kg]

            C.r_G_CM_C = 1e-3*[ ...
                 0.30; ...
                -8.47; ...
                 0.43 ...
            ];

            % The SOURCE team provided moments relative to the CoM. Products of
            % inertia were not supplied. Exact available diagonal values are
            % retained here instead of the rounded values used previously.
            C.J_reportado_C = [ ...
                0.05255722, NaN,        NaN; ...
                NaN,        0.01008568, NaN; ...
                NaN,        NaN,        0.05245392 ...
            ];
            C.J_modelo_C = diag([0.05255722,0.01008568,0.05245392]);
            C.productos_inercia_reportados = false;
            C.inertia_metric_mode = 'diagonal';
            C.data_note = 'CAD-derived mass properties; reported moments relative to CoM.';

        case {'M3SAT','M3'}
            C.id = 'M3SAT';
            C.nombre = 'M3Sat';
            C.M = 4.5624; % [kg], component-mass sum in the supplied worksheet

            % Supplied worksheet:
            %   COV = [51, 54, 170.25] mm
            %   CM  = [49.04000399, 52.20677776, 157.7160245] mm
            % After the adopted source-to-C orientation:
            C.r_G_CM_C = 1e-3*[ ...
                -1.959996009; ...
               -12.53397554; ...
                 1.793222243 ...
            ];

            % The supplied worksheet calculates the three moments about COV.
            % Source-frame moments about COV, reordered into C axes:
            J_G_C = diag([ ...
                0.01672892768, ...
                0.0005224709483, ...
                0.01696419624 ...
            ]);

            % Products about COV are unavailable; zero values close the
            % numerical ground truth. Translate that model from G/COV to CM.
            r_G_CM = C.r_G_CM_C;
            I3 = eye(3);
            shift_G_to_CM = C.M*((r_G_CM'*r_G_CM)*I3-r_G_CM*r_G_CM');
            C.J_modelo_C = J_G_C-shift_G_to_CM;

            % Only the centroidal diagonal moments are retained as available
            % reference components for the reported error metric.
            C.J_reportado_C = [ ...
                C.J_modelo_C(1,1), NaN, NaN; ...
                NaN, C.J_modelo_C(2,2), NaN; ...
                NaN, NaN, C.J_modelo_C(3,3) ...
            ];
            C.productos_inercia_reportados = false;
            C.inertia_metric_mode = 'diagonal';
            C.data_note = 'Component worksheet moments translated from COV to CM.';

        case {'INFLATESAIL','INFLATE','GB06'}
            C.id = 'INFLATESAIL';
            C.nombre = 'InflateSail';
            C.M = 3.4989; % [kg], measured flight-model mass

            % Design-report geometric center and CoG in QB50 axes:
            %   G  = [50.00, 50.00, 163.75] mm
            %   CM = [50.32, 50.20, 177.00] mm
            % The adopted physical orientation is
            %   x_C = x_QB50, y_C = -z_QB50, z_C = y_QB50,
            % which gives x_B > 0 and y_B < 0 under the common platform mount.
            C.r_G_CM_C = 1e-3*[ ...
                 0.32; ...
               -13.25; ...
                 0.20 ...
            ];

            % QB50 principal moments reordered into the normalized C axes and
            % used directly as the centroidal reference tensor.
            C.J_reportado_C = diag([0.0290,0.0059,0.0290]);
            C.J_modelo_C = C.J_reportado_C;
            C.productos_inercia_reportados = true;
            C.inertia_metric_mode = 'full';
            C.data_note = 'Mission mass-property data expressed in the common C frame.';

        otherwise
            error(['Unknown CubeSat: %s. Supported options: ', ...
                'INNOCUBE, TPA1, SOURCE, M3SAT, INFLATESAIL.'],name);
    end

    % =============================================================
    % TRANSFORMATION TO PLATFORM BODY AXES
    % =============================================================

    C.r_G_CM_B = C.R_BC*C.r_G_CM_C;
    C.r_CR_CM_B = C.r_CR_G_B+C.r_G_CM_B;

    C.J_CM_B = C.R_BC*C.J_modelo_C*C.R_BC';

    I3 = eye(3);
    r = C.r_CR_CM_B;
    C.J_CR_B = C.J_CM_B+C.M*((r'*r)*I3-r*r');
    C.S_CR_B = C.M*C.r_CR_CM_B;

    % Source-supported comparison tensor in B axes. For partial tensors,
    % preserve NaN in components not available from the source-derived data.
    if C.productos_inercia_reportados
        C.J_reportado_B = C.R_BC*C.J_reportado_C*C.R_BC';
    else
        C.J_reportado_B = nan(3);
        C.J_reportado_B(1,1) = C.J_reportado_C(3,3);
        C.J_reportado_B(2,2) = C.J_reportado_C(2,2);
        C.J_reportado_B(3,3) = C.J_reportado_C(1,1);
    end

    C.mascara_J_C = isfinite(C.J_reportado_C);
    C.mascara_J_B = isfinite(C.J_reportado_B);
    C.nota_tensor = C.data_note;

    % Sanity check for the selected common mounting orientation.
    if ~(C.r_G_CM_B(1) > 0 && C.r_G_CM_B(2) < 0)
        error(['The selected orientation for %s violates the common mounting ', ...
            'criterion x_B > 0 and y_B < 0.'],C.nombre);
    end
end
