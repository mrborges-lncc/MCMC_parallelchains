function [G, rock] = load_geology(Lx,Ly,Lz,nx,ny,nz,namex,namey,namez,...
    phi,rho,beta,fat,nini)
    TOL= 1.0e-07;
    dx = Lx/double(nx);
    dy = Ly/double(ny);
    dz = Lz/double(nz);
    G  = computeGeometry(cartGrid([nx ny nz],[Lx Ly Lz]*meter^3));
    [yx, Llx, Lly, Llz, nnx, nny, nnz] = perm_reader(namex,nini,'3D');
    if(abs(Llx-Lx)>TOL || abs(Lly-Ly)>TOL || abs(Llz-Lz)>TOL ) error('\nWrong dimention\n'); end
    [yy, Llx, Lly, Llz, nnx, nny, nnz] = perm_reader(namey,nini,'3D');
    if(abs(Llx-Lx)>TOL || abs(Lly-Ly)>TOL || abs(Llz-Lz)>TOL ) error('\nWrong dimention\n'); end
    [yz, Llx, Lly, Llz, nnx, nny, nnz] = perm_reader(namez,nini,'3D');
    if(abs(Llx-Lx)>TOL || abs(Lly-Ly)>TOL || abs(Llz-Lz)>TOL ) error('\nWrong dimention\n'); end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Log-permeability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    K = [beta*exp(rho*yx) beta*exp(rho*yy) beta*exp(rho*yz)]*fat*milli*darcy;
    clear yx yy yz
    if(nx == nnx && ny == nny && nz==nnz)
        fprintf('\nOK\n');
    else
        if(nnx > nx || nny > ny || nnz > nz )
            error('Upscaling required')
        else
            Geo = computeGeometry(cartGrid([nnx nny nnz],[Lx Ly Lz]*meter^3));
            K = mappingGeo(Geo,G,K,Lx/double(nnx),Ly/double(nny),Lz/double(nnz));
        end
        clear Geo
    end
    rock = makeRock(G, K, phi);
    clear K
end

