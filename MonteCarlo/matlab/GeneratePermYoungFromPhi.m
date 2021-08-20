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
addpath ../../MATLAB_POSPROC/
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
fim = 1999;
prt = 1; % print for simulation
printa = 10;
Ss  = 68798.072;
c   = 10;
ckc = 1/(c * Ss^2);
a   = 7.12;
E0  = 5e10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nome= 'permKC';
permrho = 0.81;
permbeta= 8.950e-14;      %% Factor to permeability
permvar = '\kappa';
phirho  = 0.23;
phibeta = 0.146;
phivar  = '\phi';
Erho    = 0.23;
Ebeta   = 1.722e10;
Evar    = '\mathsf{E}';
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
home = '~/fields/campos/';
permname = 'perm';
phiname  = 'phi';
Ename    = 'E';
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
paramE = [];
medias = [];
desvios= [];
fat = 1/(milli * darcy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup material parameters for Biot and mechanics %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = ini:fim
    snum    = num2str(n,'%d');
    Y   = load_perm(G,phifilen,phifilen,phifilen,depth,n,nD);
%     Y   = Y(:,1);
    Y   = Y(:,1) - mean(Y(:,1));
    Y   = Y/std(Y);
    if printa == 1, NORMAL(Y,mean(Y),std(Y),'$Y$'); end 
    meanY = mean(Y);
    stdY  = std(Y);
    phi = phibeta * exp(phirho * Y);
    if printa == 1, LOGNORMAL(phi,mean(phi),std(phi),'$\phi$'); end 
    if (max(phi) > 0.95 || min(phi) < 0.0)
        error('Problema nos valores max e min de phi')
    end
    %
    perm  = ckc * (phi.^3) ./ ((1.0 - phi).^2);
    sv    = 1e-03;
    perm  = perm + lhsnorm(0,sv,G.cells.num).*perm;
    if printa == 1, LOGNORMAL(perm,mean(perm),std(perm),'$\kappa$'); end 
    mk    = mean((perm));
    vk    = var((perm));
    mlk   = mean(log(perm));
    rho   = var(log(perm));
    beta  = exp(mlk);
    param = [param; beta rho];
    fprintf('\n==============================================================\n')
    fprintf('Mean Y......: %4.3f    \t | \t std Y......: %4.2f   \n',meanY,stdY);
    fprintf('Mean phi....: %4.2f    \t\t | \t std phi....: %4.2f   \n',mean(phi),std(phi));
    fprintf('Mean perm...: %4.3f mD \t | \t std perm...: %4.3f mD\n',mk*fat,sqrt(vk)*fat);
    %% Sprigg’s representation
    E  = E0 * exp(-a * phi);
    sv = 2.5e-05;
    E  = E + lhsnorm(0,sv,G.cells.num).*E;
    mu_E = mean(E);
    std_E= std(E);
    muY  = log(mu_E) - (1/2) * log((std_E^2)/(mu_E^2) + 1);
    varY = log((std_E^2)/(mu_E^2) + 1);
    stdY = sqrt(varY);
    betaE= exp(muY);
    rhoE = sqrt(varY);
    if printa == 1, NORMAL(E,mean(E),std(E),'$\mathsf{E}$'); end 
    paramE = [paramE; betaE rhoE];
    fprintf('Mean E.......: %4.3fe9 N/m^2 \t | \t std E......: %4.3fe9 N/m^2\n',mean(E)/1e9,std(E)/1e9);
    fprintf('==============================================================\n')
    %
    if printa == 1, NORMAL(reverseKlog(E,Ebeta,Erho),...
            mean(reverseKlog(E,Ebeta,Erho)),std(reverseKlog(E,Ebeta,Erho)),'$Y_{\mathsf{E}}$'); end
    savefields3(Lx,Ly,Lz,nx,ny,nz,1,reverseKlog(perm,beta,rho),...
        reverseKlog(phi,phibeta,phirho),reverseKlog(E,Ebeta,Erho),n,...
        home,permname,phiname,Ename,prt)
    medias= [medias; mean(phi) mean(perm) mean(E)];
    desvios= [desvios; std(phi) std(perm) std(E)];
end
pbeta = 1.0; % pbeta = permbeta;
plotKC(phi,perm/pbeta,0)
base=['../../figuras/perm_phi_KC_' nome];
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r600', base);

plotEphi(phi,E,E0);
base=['../../figuras/perm_phi_E_' nome];
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r600', base);
beta = mean(param(:,1));
rho  = mean(param(:,2));
betaE= mean(paramE(:,1));
rhoE = mean(paramE(:,2));
fprintf('\n==============================================================\n')
fprintf('==============================================================\n')
fprintf('PHI =====> Mean = %5.4f\t Std = %5.4f\n',mean(medias(:,1)),mean(desvios(:,1)))
fprintf('PERM ====> Mean = %5.4e\t Std = %5.4e\n',mean(medias(:,2)),mean(desvios(:,2)))
fprintf('E =======> Mean = %5.4e\t Std = %5.4e\n',mean(medias(:,3)),mean(desvios(:,3)))
fprintf('==============================================================\n')
fprintf('==============================================================\n')
fprintf('PHI =====> beta = %5.4f\t rho = %5.3f\n',phibeta,phirho)
fprintf('PERM ====> beta = %5.4e\t rho = %5.3f\n',beta,rho)
fprintf('E =======> beta = %5.4e\t rho = %5.3f\n',betaE,rhoE)
fprintf('==============================================================\n')
fprintf('==============================================================\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(22)
hist(param(:,1),30);
figure(23)
hist(param(:,2),30);
figure(24)
hist(paramE(:,1),30);
figure(25)
hist(paramE(:,2),30);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
