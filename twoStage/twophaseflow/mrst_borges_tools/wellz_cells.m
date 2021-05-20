function [wcells] = wellz_cells(px,py,G,Lx,Ly,nx,ny)
    TOL = 0.0e-6;
    wcells = [];
    p_inj = [px py];
    dx = 0.5*Lx/nx+TOL;
    dy = 0.5*Ly/ny+TOL;
    sz = G.cells.num;
    for i=1:sz
        x = G.cells.centroids(i,1);
        y = G.cells.centroids(i,2);
        if(p_inj(1)<=(x+dx) && p_inj(1)>(x-dx))
            if(p_inj(2)<=(y+dy) && p_inj(2)>(y-dy))
            wcells = [wcells i];
        end
    end
end

