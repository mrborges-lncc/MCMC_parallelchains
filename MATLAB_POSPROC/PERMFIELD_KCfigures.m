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
permrho = 0.4432;
permbeta= 5.8675e-14;      %% Factor to permeability
permvar = '\kappa';
phirho  = 0.2;
phibeta = 0.12;
phivar  = '\phi';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[FILENAME, PATHNAME] =uigetfile({'~/MCMC_parallelchains/twophaseflow/exp/fields/*.dat'}, 'LOAD DATA PORO FIELD');
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
[FILENAME, PATHNAME] =uigetfile({'~/MCMC_parallelchains/twophaseflow/exp/fields/*.dat'}, 'LOAD DATA PERM FIELD');
%[FILENAME, PATHNAME] =uigetfile({'~/Dropbox/mrstBorges/out/*.dat'}, 'LOAD DATA');
filen=sprintf('%s%s', PATHNAME,FILENAME);
lf = length(filen);
filem = filen(1:end-4);
k = length(filem);
while filem(k) ~= '_'
    k = k-1;
end
permfilen = filen(1:k);
nini = int32(str2num(filem(k+1:end)));
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rock model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rock = makeRock(G, perm, phi);
mK   = mean(rock.perm/(milli*darcy));
sK   = std(rock.perm/(milli*darcy));
fprintf('\n==============================================================\n')
fprintf('Mean Yphi....: %4.1f    \t | \t std phi....: %4.1f   \n',mean(phi),std(phi));
fprintf('Mean Yperm...: %4.1f mD \t | \t std K......: %4.1f mD\n',mK,sK);
fprintf('==============================================================\n')
clear K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotKC(phi,perm,0)
base=['../figuras/perm_phi_KC_' nome]
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r600', base);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
