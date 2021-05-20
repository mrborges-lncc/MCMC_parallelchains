function [inj,prod1,prod2,prod3,prod4] = fivespot_wells(G,Lx,Ly,nx,ny)
    df  = 0.0;
    dx  = 0.25*Lx/nx;
    dy  = 0.25*Ly/ny;
    Lx  = Lx-dx;
    Ly  = Ly-dy;
    px  = Lx/2; py = Ly/2;
    inj = wellz_cells(px,py,G,Lx,Ly,nx,ny);
    px  = Lx - df; py = Ly - df;
    prod1 = wellz_cells(px,py,G,Lx,Ly,nx,ny);
    px  = dx + df; py = Ly - df;
    prod2 = wellz_cells(px,py,G,Lx,Ly,nx,ny);
    px  = Lx - df; py = dy + df;
    prod3 = wellz_cells(px,py,G,Lx,Ly,nx,ny);
    px  = dx + df; py = dy + df;
    prod4 = wellz_cells(px,py,G,Lx,Ly,nx,ny);
end
    