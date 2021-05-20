function [K] = load_perm(G,namex,namey,namez,depth,nini,nD)
    dp = depth;
    TOL= 1.0e-07;
    L  = max(G.nodes.coords);
    Lx = L(1); Ly = L(2); Lz = L(3);
    L  = min(G.nodes.coords);
    Lx0 = L(1); Ly0 = L(2); Lz0 = L(3);
    n  = G.cartDims;
    nx = int16(n(1)); ny = int16(n(2)); nz = int16(n(3));
    dx = Lx/double(nx);
    dy = Ly/double(ny);
    dz = Lz/double(nz);
    [yx, Llx, Lly, Llz, nnx, nny, nnz] = perm_reader(namex,nini,nD);
    if(nD == '2D') Llz = Lz; dp = 0.0; end;
    if(abs(Llx-Lx)>TOL || abs(Lly-Ly)>TOL || abs(Llz-Lz+dp)>TOL )
        error('Wrong dimention'); 
    end
    [yy, Llx, Lly, Llz, nnx, nny, nnz] = perm_reader(namey,nini,nD);
    if(nD == '2D') Llz = Lz; dp = 0.0; end;
    if(abs(Llx-Lx)>TOL || abs(Lly-Ly)>TOL || abs(Llz-Lz+dp)>TOL )
        error('Wrong dimention');
    end
    [yz, Llx, Lly, Llz, nnx, nny, nnz] = perm_reader(namez,nini,nD);
    if(nD == '2D') Llz = Lz; dp = 0.0; end;
    if(abs(Llx-Lx)>TOL || abs(Lly-Ly)>TOL || abs(Llz-Lz+dp)>TOL )
        error('Wrong dimention');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Log-permeability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    K = [yx yy yz];
    clear yx yy yz
    if(nx == nnx && ny == nny && nz==nnz)
        fprintf('\nFile number %d OK\n',nini);
    else
        if(nnx > nx || nny > ny || nnz > nz )
            error('Upscaling required')
        else
            Geo  = cartGrid([nnx nny nnz],[Lx Ly (Lz-Lz0)]*meter^3);
            Geo.nodes.coords(:, 3) = depth + Geo.nodes.coords(:, 3)*meter;
            Geo.nodes.coords(:, 2) = Geo.nodes.coords(:, 2)*meter;
            Geo.nodes.coords(:, 1) = Geo.nodes.coords(:, 1)*meter;
            Geo  = computeGeometry(Geo);
            K = mappingGeo(Geo,G,K,(Lx-Lx0)/double(nnx),...
                (Ly-Ly0)/double(nny),(Lz-Lz0)/double(nnz));
        end
        clear Geo
    end
end

