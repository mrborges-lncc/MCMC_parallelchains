% Single phase flow Simulator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 25/03/2020
% Based on MRST 2019
% Author: Marcio Borges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

tstart =  tic;

addpath ~/Dropbox/mrst_borges_tools/
addpath ~/Dropbox/mrst-2021a/
startup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx  = 510.0;
Ly  = 510.0;
Lz  = 20.0;
nx  = 51;
ny  = 51;
nz  = 5;
depth = 1e3;
ini = 0;
fim = 0;
prt = 1; % print for simulation
ckc = 2.5e-11;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nome= 'permKC';
permrho = 0.42;
permbeta= 5.7e-14;      %% Factor to permeability
permvar = '\kappa';
phirho  = 0.275;
phibeta = 0.125;
phivar  = '\phi';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx  = Lx/double(nx);
dy  = Ly/double(ny);
dz  = Lz/double(nz);
G   = cartGrid([nx ny nz],[Lx Ly Lz]*meter^3);
G.nodes.coords(:, 3) = depth + G.nodes.coords(:, 3)*meter;
G.nodes.coords(:, 2) = G.nodes.coords(:, 2)*meter;
G.nodes.coords(:, 1) = G.nodes.coords(:, 1)*meter;
G   = computeGeometry(G);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
home = '~/Dropbox/PROJETO_MCMC_RIGID/MCMC_parallelchains/MonteCarlo/twophaseflow/fields';
permname = 'perm';
phiname  = 'phi';
[FILENAME, PATHNAME] =uigetfile({'~/fields/*.dat'}, 'LOAD DATA PORO FIELD');
%[FILENAME, PATHNAME] =uigetfile({'~/MCMC_parallelchains/twophaseflow/exp/fields/*.dat'}, 'LOAD DATA PORO FIELD');
%[FILENAME, PATHNAME] =uigetfile({'~/Dropbox/mrstBorges/out/*.dat'}, 'LOAD DATA');
filen=sprintf('%s%s', PATHNAME,FILENAME);
lf = length(filen);
filem = filen(1:end-4);
k = length(filem);
while filem(k) ~= '_'
    k = k-1;
end
phifilen = filen(1:k);
nini = int32(str2num(filem(k+1:end)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
   require incomp ad-mechanics
catch %#ok<CTCH>
   mrstModule add incomp ad-mechanics
end
verbose = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nD = '3D';
color = 'none';
%color = 'k';
vw  = [-35 20];
et = 3;
param  = [];
medias = [];
desvios= [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MPHI =[];
MPER =[];
for n = ini:fim
    snum    = num2str(n,'%d');
    Y   = load_perm(G,phifilen,phifilen,phifilen,depth,n,nD);
    phi = phibeta * exp(phirho * Y(:,1));
    if (max(phi) > 0.95 || min(phi) < 0.0)
        error('Problema nos valores max e min de phi')
    end
    perm = ckc * (phi.^3) ./ ((1.0 - phi).^2);
    sv   = 1e-03;
    perm = perm + lhsnorm(0,sv,G.cells.num).*perm;
    mk   = mean((perm));
    vk   = var((perm));
    mlk   = mean(log(perm));
    rho   = var(log(perm));
    beta  = exp(mlk);
    param = [param; beta rho];
    medias= [medias; mean(phi) mean(perm)];
    desvios= [desvios; std(phi) std(perm)];
    fprintf('\n==============================================================\n')
    fprintf('Mean phi.....: %4.3f    \t | \t std phi....: %4.2e   \n',mean(phi),std(phi));
    fprintf('Mean perm....: %4.3f mD \t | \t std K......: %4.3f mD\n',mk/(milli * darcy()),sqrt(vk)/(milli * darcy()));
    fprintf('==============================================================\n')
    savefields(Lx,Ly,Lz,nx,ny,nz,1,reverseKlog(perm,beta,rho),...
        reverseKlog(phi,phibeta,phirho),n,home,permname,phiname,prt);
    MPHI = [MPHI; phi];
    MPER = [MPER; mk];
end
pbeta = 1.0; % pbeta = permbeta;
plotKC(phi,perm/pbeta,0)
base=['../figuras/perm_phi_KC_' nome];
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r600', base);
beta = mean(param(:,1));
rho  = mean(param(:,2));
fprintf('\n==============================================================\n')
fprintf('==============================================================\n')
fprintf('PHI =====> Mean = %5.4f\t Std = %5.4f\n',mean(medias(:,1)),mean(desvios(:,1)))
fprintf('PERM ====> Mean = %5.4e\t Std = %5.4e\n',mean(medias(:,2)),mean(desvios(:,2)))
fprintf('==============================================================\n')
fprintf('==============================================================\n')
fprintf('PHI =====> beta = %5.4f\t rho = %f\n',phibeta,phirho)
fprintf('PERM ====> beta = %5.4e\t rho = %f\n',beta,rho)
fprintf('==============================================================\n')
fprintf('==============================================================\n')
fprintf('<PHI>  = %5.4f   \t sigPHI = %5.4f\n',mean(MPHI),std(MPHI))
fprintf('<PERM> = %5.4f mD\t sigK   = %5.4f mD\n',mean(MPER)/(milli * darcy()),std(MPER)/(milli * darcy()))
fprintf('==============================================================\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(22)
hist(param(:,1),30);
figure(23)
hist(param(:,2),30);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
