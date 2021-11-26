function [out] = rockupscaling(Y,physical_dim,fine_mesh,coarse_mesh,tipo)
    if tipo == 'perm'
        tipo_upsc = 'sealing_';
        %% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        G = computeGeometry(cartGrid([fine_mesh],[physical_dim]*meter^3));
        uG = computeGeometry(cartGrid([ceil(fine_mesh(1)/coarse_mesh(1)) ...
            ceil(fine_mesh(2)/coarse_mesh(2)) ...
            ceil(fine_mesh(3)/coarse_mesh(3))], [1.0 1.0 1.0]));
        p  = partitionUI(G, [coarse_mesh]);
        coarseG = computeGeometry(cartGrid([coarse_mesh],[physical_dim]*meter^3));
        p  = compressPartition(p);
        CG = generateCoarseGrid(G, p);
        map_blocks = @(b) find(CG.partition == b);
        clear p
        rho  = 1.0;
        beta = 1.0;
        K    = [beta*exp(rho*Y) beta*exp(rho*Y) beta*exp(rho*Y)];
        rock = makeRock(G, K, 1.0);
        Kup  = upscalingK(CG,uG,rock,tipo_upsc);
        out  = (log(Kup(:,1)) - log(beta))/rho;
    end
    if tipo == 'poro'
        %% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        G    = computeGeometry(cartGrid([fine_mesh],[physical_dim]*meter^3));
        rho  = 1.0;
        beta = 1.0;
        poro = [beta*exp(rho*Y)];
        rock = makeRock(G, 1, poro);
        p    = partitionUI(G, [coarse_mesh]);
        p    = compressPartition(p);
        CG   = generateCoarseGrid(G, p);
        crock.poro = accumarray(p, rock.poro)./accumarray(p,1);
        out  = (log(crock.poro) - log(beta))/rho;
    end
end

