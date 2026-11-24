function T = validate_spacecraft_database()
%VALIDATE_SPACECRAFT_DATABASE Validate the five spacecraft reference cases.
%
% Checks:
%   1) the five expected validation spacecraft can be loaded;
%   2) the selected mounting orientation satisfies x_B > 0 and y_B < 0;
%   3) centroidal model tensors are symmetric and positive definite;
%   4) CR-referenced tensors reconstructed by the parallel-axis theorem are
%      symmetric and positive definite;
%   5) the source-supported comparison mask is nonempty.
%
% The returned table is also written to
%   spacecraft_reference_database.csv
% in the current folder.

    ids = {'INNOCUBE','TPA1','SOURCE','M3SAT','INFLATESAIL'};
    n = numel(ids);

    Spacecraft = strings(n,1);
    Mass_kg = zeros(n,1);
    xB_mm = zeros(n,1);
    yB_mm = zeros(n,1);
    zB_mm = zeros(n,1);
    MinEig_JCM = zeros(n,1);
    MinEig_JCR = zeros(n,1);
    ReferenceComponents = zeros(n,1);
    Status = strings(n,1);

    for k = 1:n
        C = cubesat_parameters(ids{k});

        assert(norm(C.J_modelo_C-C.J_modelo_C','fro') < 1e-12, ...
            '%s J_CM is not symmetric.',C.nombre);
        assert(norm(C.J_CR_B-C.J_CR_B','fro') < 1e-12, ...
            '%s J_CR is not symmetric.',C.nombre);

        eig_cm = eig((C.J_modelo_C+C.J_modelo_C')/2);
        eig_cr = eig((C.J_CR_B+C.J_CR_B')/2);
        assert(all(eig_cm > 0),'%s J_CM is not positive definite.',C.nombre);
        assert(all(eig_cr > 0),'%s J_CR is not positive definite.',C.nombre);
        assert(C.r_G_CM_B(1) > 0 && C.r_G_CM_B(2) < 0, ...
            '%s violates x_B > 0, y_B < 0.',C.nombre);
        assert(any(C.mascara_J_C(:)), ...
            '%s has no source-supported inertia components.',C.nombre);

        Spacecraft(k) = string(C.nombre);
        Mass_kg(k) = C.M;
        xB_mm(k) = 1e3*C.r_G_CM_B(1);
        yB_mm(k) = 1e3*C.r_G_CM_B(2);
        zB_mm(k) = 1e3*C.r_G_CM_B(3);
        MinEig_JCM(k) = min(eig_cm);
        MinEig_JCR(k) = min(eig_cr);
        ReferenceComponents(k) = nnz(C.mascara_J_C);
        Status(k) = "OK";
    end

    T = table(Spacecraft,Mass_kg,xB_mm,yB_mm,zB_mm, ...
        MinEig_JCM,MinEig_JCR,ReferenceComponents,Status);

    writetable(T,'spacecraft_reference_database.csv');
    disp(T);
end
