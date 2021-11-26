function [perm] = upscalingK(cg,ug,rck,tipo_up)
    N = cg.cells.num;
    d = cg.griddim;
    map  = @(b) find(cg.partition == b);
    perm = zeros(N,d);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Boundary conditiions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % name tags to locate the correct pair of boundary faces for each flow problem
    bcsides = {'XMin', 'XMax'; 'YMin', 'YMax'; 'ZMin', 'ZMax'};
    for j = 1:d
        bcl{j} = pside([], ug, bcsides{j,1}, 0); %store a template of all boundary conditions
        bcr{j} = pside([], ug, bcsides{j,2}, 0);
    end
    Dp = {1.0, 0.0}; % pressure drop
    L  = max(ug.faces.centroids) - min(ug.faces.centroids); %characteristic length in the axial direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fluid = initSingleFluid('mu',1.0,'rho', 1.0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if tipo_up == 'sealing_'
        for j=1:N
            ck = rck.perm(map(j),:);
            crock = makeRock(ug, ck, 1.0);
            perm(j,:) = upscalePermeabilityFixed_mrb(ug,Dp{1},fluid,crock,bcl,bcr,L);
        end
    end
    if tipo_up == 'periodic'
        [Gp, bcp] =  makePeriodicGridMulti3d(ug,bcl,bcr,Dp);
        for el =1:N
            ck = rck.perm(map(el),:);
            crock = makeRock(ug, ck, 1.0);
            kk = upscalePermeabilityPeriodic_mrb(Gp, ug, bcp, Dp{1}, fluid, crock, L)
            perm(el,:) = diag(kk).';
        end
    end
end

