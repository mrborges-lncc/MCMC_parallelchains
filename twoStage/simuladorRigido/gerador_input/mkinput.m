clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nome = '../exp_coarse/input.in';
dim = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% phisical dimensions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx = 1.00;%304.8;  % [L]
Ly = 1.00;%0.6096;   % [L]
Lz = 1;%304.8;  % [L]
%%%% mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = 50;
ny = 50;
nz = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
porosity = 0.20;
if(dim==2)
    VolTotalPoroso = Lx*Ly*porosity
end
if(dim==3)
    VolTotalPoroso = Lx*Ly*Lz*porosity
end
T = 365.40; %TEMPO
pV = 0.5; % porcentagem do Volume poroso total
Q=VolTotalPoroso*pV/T % L^3/T
Q=1.0
%Q=50.0; %500BBL/DAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(dim==3)
%    slabx(Lx,Ly,Lz,nx,ny,nz,Q,nome)
%    slaby(Lx,Ly,Lz,nx,ny,nz,Q,nome)
    slabz(Lx,Ly,Lz,nx,ny,nz,Q,nome)
%    fivespot(Lx,Ly,Lz,nx,ny,nz,Q,nome)
%    fivespotz(Lx,Ly,Lz,nx,ny,nz,Q,nome)
end
if(dim==2)
%    slabx2(Lx,Ly,Lz,nx,ny,nz,Q,nome)
%    slaby2(Lx,Ly,Lz,nx,ny,nz,Q,nome)
    fivespot2(Lx,Ly,Lz,nx,ny,nz,Q,nome)
%    fivespot12(Lx,Ly,Lz,nx,ny,nz,Q,nome)
%    fivespot21(Lx,Ly,Lz,nx,ny,nz,Q,nome)
end
%clear
