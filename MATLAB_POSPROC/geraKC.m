clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath ~/Dropbox/mrst_borges_tools/
addpath ~/Dropbox/mrst-2021a/
startup

cd ~/Dropbox/mrstBorges/

try
   require ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui
catch %#ok<CTCH>
   mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui
end
verbose = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx  = 510.0;
Ly  = 510.0;
Lz  = 20.0;
nx  = 51;
ny  = 51;
nz  = 5;
depth = 1.0e03*meter;   %% depth until the top of reservoir
dx  = Lx/double(nx);
dy  = Ly/double(ny);
dz  = Lz/double(nz);
G   = cartGrid([nx ny nz],[Lx Ly Lz]*meter^3);
G.nodes.coords(:, 3) = depth + G.nodes.coords(:, 3)*meter;
G.nodes.coords(:, 2) = G.nodes.coords(:, 2)*meter;
G.nodes.coords(:, 1) = G.nodes.coords(:, 1)*meter;
G   = computeGeometry(G);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filenx = '~/MCMC_parallelchains/twophaseflow/exp/fields/perm_ref_';
fileny = '~/MCMC_parallelchains/twophaseflow/exp/fields/perm_ref_';
filenz = '~/MCMC_parallelchains/twophaseflow/exp/fields/perm_ref_';
phi  = 0.2;                       %% Porosity
KC   = true;                      %% Kozeny-Carman relation
K = [100 100 100]*milli*darcy;
vK= ([100 100 100]*milli*darcy.^2);
beta = K; %% Factor to permeability
rho  = 1.0;
phibeta = phi;
phirho  = 0.23;
spriggs = 3.50; % E = E0 * exp(-spriggs * phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup material parameters for Biot and mechanics %%%%%%%%%%%%%%%%%%%%%%%
E      = 1 * giga * Pascal;         %% Young's module
nu     = 0.3;                       %% Poisson's ratio
alpha  = 1.0;                       %% Biot's coefficient
Bulk   = mean(E)/(3.0*(1.0 - 2*mean(nu))); %% Bulk modulus
CompR  = 1.0/Bulk;                  %% Rock compressibility
CompW  = 2.0e-06/psia;              %% Water compressibility
CompO  = 1.0e-05/psia;              %% Oil compressibility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GEOLOGIC MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nini   = 0;
nD     = '3D';
heter  = 1;      %% if 0 => homog.; 1 => heter.; 2 => read *.mat
printa = 1;
nome   = 'ref';
[K, phi, E, nu, alpha, Bulkm, CompR] = geologiaMECH(G,nini,filenx,...
        fileny,filenz,nD,K,vK,beta,rho,phibeta,phirho,E,spriggs,nu,...
        alpha,depth,KC,heter,printa,nome);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rock model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rock = makeRock(G, K, phi);
pv   = poreVolume(G, rock);
mK   = mean(rock.perm/(milli*darcy));
sK   = std(rock.perm/(milli*darcy));
fprintf('\n==============================================================\n\n')
fprintf('Mean k_x....: %4.1f mD \t | \t std k_x....: %4.1f mD\n',mK(1),sK(1));
fprintf('Mean k_y....: %4.1f mD \t | \t std k_y....: %4.1f mD\n',mK(2),sK(2));
fprintf('Mean k_z....: %4.1f mD \t | \t std k_z....: %4.1f mD\n',mK(3),sK(3));
fprintf('\n==============================================================\n')
clear K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color = 'none';
%color = 'k';
vw  = [-35 20];
lim = [1e-15 1e-12];
et  = 1;
if printa == 1 && heter == 1 
    plot_rock(rock.perm(:,1),G,'Yn','$\kappa_x$',color,lim,vw,1);
    base=['./figuras/permKx_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    plot_rock(rock.perm(:,2),G,'Yn','$\kappa_y$',color,lim,vw,2);
    base=['./figuras/permKy_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    plot_rock(rock.perm(:,3),G,'Yn','$\kappa_z$',color,lim,vw,3);
    base=['./figuras/permKz_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    %plot_rock(rock.poro,G,'Yn','$\phi$',color,lim,vw,4);
    plot_rock_poro(rock.poro,G,'Yn',1,1,'$\phi$',color,[0 0],vw,4)
    base=['./figuras/phi_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    pause(et); clf; close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
home = '~/MCMC_parallelchains/twophaseflow/exp/fields';
prt  = 1;
n    = 0;
name = 'perm_ref';
save_field(Lx,Ly,Lz,nx,ny,nz,1,...
    reverseKlog(rock.perm(:,1),beta,rho),...
    reverseKlog(rock.perm(:,2),beta,rho),...
    reverseKlog(rock.perm(:,3),beta,rho),n,home,name,prt);
name = 'poro_ref';
save_field(Lx,Ly,Lz,nx,ny,nz,1,...
        reverseKlog(rock.poro(:,1),phibeta,phirho),...
        reverseKlog(rock.poro(:,1),phibeta,phirho),...
        reverseKlog(rock.poro(:,1),phibeta,phirho),n,home,name,prt);
clear;