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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nome= 'permKC';
permrho = 0.597;
permbeta= 9.1098e-14;      %% Factor to permeability
permvar = '\kappa';
phirho  = 0.23;
phibeta = 0.146;
phivar  = '\phi';
Erho    = 0.201;
Ebeta   = 9.9313e09;
E0      = 2.5e10;
Evar    = '\mathsf{E}';
home = '~/Dropbox/PROJETO_MCMC_RIGID/MCMC_parallelchains/';
homef= '~/Dropbox/PROJETO_MCMC_RIGID/paper/figuras/';
home = '~/MCMC_parallelchains/';
homef= '~/MCMC_parallelchains/figuras/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[FILENAME, PATHNAME] =uigetfile({[home 'twophaseflow/exp/fields/*.dat']}, 'LOAD DATA PORO FIELD');
[FILENAME, PATHNAME] =uigetfile({'~/fields/campos/phi*.dat'}, 'LOAD DATA PORO FIELD');
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
%[FILENAME, PATHNAME] =uigetfile({[home 'twophaseflow/exp/fields/*.dat']}, 'LOAD DATA PERM FIELD');
[FILENAME, PATHNAME] =uigetfile({'~/fields/campos/pe*.dat'}, 'LOAD DATA PERM FIELD');
filen=sprintf('%s%s', PATHNAME,FILENAME);
lf = length(filen);
filem = filen(1:end-4);
k = length(filem);
while filem(k) ~= '_'
    k = k-1;
end
permfilen = filen(1:k);
nini2 = int32(str2num(filem(k+1:end)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[FILENAME, PATHNAME] =uigetfile({[home 'twophaseflow/exp/fields/*.dat']}, 'LOAD DATA PERM FIELD');
[FILENAME, PATHNAME] =uigetfile({'~/fields/campos/E*.dat'}, 'LOAD DATA PERM FIELD');
filen=sprintf('%s%s', PATHNAME,FILENAME);
lf = length(filen);
filem = filen(1:end-4);
k = length(filem);
while filem(k) ~= '_'
    k = k-1;
end
Efilen = filen(1:k);
nini3 = int32(str2num(filem(k+1:end)));
if nini ~= nini2 || nini ~= nini3
    error('Campos não compatíveis')
end
try
   require incomp ad-mechanics
catch %#ok<CTCH>
   mrstModule add incomp ad-mechanics
end
verbose = true;

nD = '3D';
color = 'none';
%color = 'k';
vw  = [-35 20];
et = 3;
%%
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
%% GEOLOGIC MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHI
Y    = load_perm(G,phifilen,phifilen,phifilen,depth,nini,nD);
phiY = Y(:,1);
phi  = phibeta .* exp(phirho * phiY);
%% PERM
Y    = load_perm(G,permfilen,permfilen,permfilen,depth,nini,nD);
permY= Y(:,1);
perm = permbeta .* exp(permrho * permY);
%% YOUNG
Y    = load_perm(G,Efilen,Efilen,Efilen,depth,nini,nD);
EY= Y(:,1);
E = Ebeta .* exp(Erho * EY);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rock model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rock = makeRock(G, perm, phi);
mK   = mean(rock.perm/(milli*darcy));
sK   = std(rock.perm/(milli*darcy));
fprintf('\n==============================================================\n')
fprintf('Mean Yphi....: %4.3f    \t | \t std phi....: %4.3f   \n',mean(phi),std(phi));
fprintf('Mean Yperm...: %4.3f mD \t | \t std K......: %4.3f mD\n',mK,sK);
fprintf('==============================================================\n')
clear K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotKC(phi,perm,0)
base=[homef 'perm_phi_KC_' nome]
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r600', base);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotEphi(phi,E,E0);
base=[homef 'perm_phi_E_' nome]
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r600', base);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotKC_Y(phiY,permY)
base=[homef 'Yperm_Yphi_KC_' nome]
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r600', base);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotE_Y(phiY,EY)
base=[homef 'Yperm_Yphi_YE_' nome]
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r600', base);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all