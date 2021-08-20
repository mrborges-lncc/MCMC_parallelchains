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
ini = 0;
fim = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho = 0.81;
beta= 8.950e-14;      %% Factor to permeability
nome= 'permREF';
variav = '\kappa';
% rho = 0.20;
% beta= 0.12;
% nome= 'phiREF';
% variav = '\phi';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[FILENAME, PATHNAME] =uigetfile({'~/Dropbox/PROJETO_MCMC_RIGID/MCMC_parallelchains/twoStage/select_fields/*.dat'}, 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile({'~/MCMC_parallelchains/twophaseflow/exp/fields/*.dat'}, 'LOAD DATA');
[FILENAME, PATHNAME] =uigetfile({'~/fields/campos/*.dat'}, 'LOAD DATA');
filen=sprintf('%s%s', PATHNAME,FILENAME);
lf = length(filen);
filem = filen(1:end-4);
k = length(filem);
while filem(k) ~= '_'
    k = k-1;
end
filen = filen(1:k);
nini = int32(str2num(filem(k+1:end)))
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
for n = ini:fim
    snum = num2str(n,'%d');
    Y = load_perm(G,filen,filen,filen,depth,n,nD);
    Y = Y(:,1);
    K = beta .* exp(rho * Y);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Rock model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rock = makeRock(G, K, 0.20);
    mK   = mean(rock.perm/(milli*darcy));
    sK   = std(rock.perm/(milli*darcy));
    fprintf('\n==============================================================\n')
    fprintf('Mean Y....: %4.1f    \t | \t std Y....: %4.1f   \n',mean(Y),std(Y));
    fprintf('Mean K....: %4.1f mD \t | \t std K....: %4.1f mD\n',mK,sK);
    fprintf('==============================================================\n')
    clear K
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lim = [0 0];
    %plot_rock(Y,G,'Yn','$Y$',color,lim,vw,1);
    plot_rock_poro(rock.perm(:,1),G,'Y',beta,rho,['$Y_{' variav '}$'],...
        color,lim,vw,2);
    base=['../../figuras/Y_' nome '_' snum];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    pause(1); close all
% %     lim = [(mean(rock.perm(:,1)) - 0.0125*std(rock.perm(:,1)))...
% %         (mean(rock.perm(:,1)) + 0.25*std(rock.perm(:,1)))];
% %     plot_rock(rock.perm(:,1),G,'Yn',variav,color,lim,vw,2);
%     plot_rock_poro(rock.perm(:,1),G,'Yn',beta,rho,['$' variav '$'],...
%         color,lim,vw,12);
%     base=['../figuras/' nome '_' snum];
%     set(gcf,'PaperPositionMode','auto');
%     print('-depsc','-r600', base);
%     pause(1); close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
