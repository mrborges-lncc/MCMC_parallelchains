function [inj,prod1,prod2,prod3,prod4] = mrb_wells(G,Lx,Ly,nx,ny,...
    p0,p1,p2,p3,p4)
    px    = p0(1); py = p0(2);
    inj   = wellz_cells(px,py,G,Lx,Ly,nx,ny);
    px    = p1(1); py = p1(2);
    prod1 = wellz_cells(px,py,G,Lx,Ly,nx,ny);
    px    = p2(1); py = p2(2);
    prod2 = wellz_cells(px,py,G,Lx,Ly,nx,ny);
    px    = p3(1); py = p3(2);
    prod3 = wellz_cells(px,py,G,Lx,Ly,nx,ny);
    px    = p4(1); py = p4(2);
    prod4 = wellz_cells(px,py,G,Lx,Ly,nx,ny);
end
    