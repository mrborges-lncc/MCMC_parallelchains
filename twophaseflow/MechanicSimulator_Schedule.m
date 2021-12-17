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

%cd ~/Dropbox/mrstBorges/

try
   require ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui
catch %#ok<CTCH>
   mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui
end

nD = '3D';
color = 'none';
% color = 'k';
vw  = [-35 20];
%vw  = [0 90];
verbose = true;
grav  = 1;
lim   = [0 0];
heter = 1;      %% if 0 => homog.; 1 => heter.; 2 => read *.mat
printa= 1;
salva = 1;
monitorpres = 1;  %% if == 1 pressure monitors at some points
monitorsat  = 1;  %% if == 1 saturation monitors at some points
monitordisp = 1;  %% if == 1 displacement monitors at some points
%nome  = 'mech';
%nome  = 'MechRigid';
%nome  = 'mech';
nome  = 'teste150';
%nome   = 'mech_fluid_incomp';
et    = 0;
isCompr = true;
if exist('isCompr','var')~=1
    isCompr = false;
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx  = 510.0;
Ly  = 510.0;
Lz  = 20.0;
nx  = 51;
ny  = 51;
nz  = 5;
prod_dirichlet = false; %% if true, production wells defined as dirichlet bc
TT        = 600.00;     %% days
Tinjstart = 100;        %% days
Tinjup    = 200;        %% days
Tinjstop  = 1700;        %% days
nstep     = 150;        %% number of time steps for pressure-velocity system
nprint    = 20;         %% Number of impressions
ndata     = 150;
ndt       = 20;
[nprint nprjump] = ajusteImpress(nprint,nstep);
[ndata  njump]   = ajusteImpress(ndata,nstep);
well_r = 0.125;         %% well radius
p_at_surface = 1.0*atm;             %% Pressure at 0m cote
depth  = 7.00e03*meter;              %% depth until the top of reservoir
rhoR   = 2.70e03*kilogram/meter^3;  %% mean density of overload rocks
PRbhp  = 0.0e03*barsa;              %% production well bh pressure
vinj   = 5.0e2/day;                 %% Injection rate
vprod  = 0.0*vinj/4;                %% Production rate
overburden= 0.0e3*atm;              %% Load (overburden)
%% Setup hidrogeological parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read Gaussian fileds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filenx = '~/Dropbox/Upscaling/KLgenerator/campos/e51x51x5_20_';
fileny = '~/Dropbox/Upscaling/KLgenerator/campos/e51x51x5_20_';
filenz = '~/Dropbox/Upscaling/KLgenerator/campos/e51x51x5_20_';
filenx = '~/MCMC_par/trunk/simuladorRigido/exp/fields/e510x510x50_51x51x5_l40x40x20_';
fileny = '~/MCMC_par/trunk/simuladorRigido/exp/fields/e510x510x50_51x51x5_l40x40x20_';
filenz = '~/MCMC_par/trunk/simuladorRigido/exp/fields/e510x510x50_51x51x5_l40x40x20_';
filenx = '~/fields/e510x510x20_51x51x5_l50x50x10_';
fileny = '~/fields/e510x510x20_51x51x5_l50x50x10_';
filenz = '~/fields/e510x510x20_51x51x5_l50x50x10_';
% filenx = '~/MCMC_par/trunk/simuladorRigido/exp/fields/amostra_';
% fileny = '~/MCMC_par/trunk/simuladorRigido/exp/fields/amostra_';
% filenz = '~/MCMC_par/trunk/simuladorRigido/exp/fields/amostra_';
nini = 0;
phi  = 0.15;                      %% Porosity
KC   = true;                      %% Kozeny-Carman relation
K    = [125 125 125]*milli*darcy;
vK   = ([118 118 118]*milli*darcy).^2;
beta = K; %% Factor to permeability
rho  = 1.0;
phibeta = phi;
phirho  = 0.15;
spriggs = 2.11; % E = E0 * exp(-spriggs * phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup material parameters for Biot and mechanics %%%%%%%%%%%%%%%%%%%%%%%
E      = 1.0 * giga * Pascal;         %% Young's module
nu     = 0.2;                       %% Poisson's ratio
alpha  = 1.0;                       %% Biot's coefficient
Bulk   = mean(E)/(3.0*(1.0 - 2*mean(nu))); %% Bulk modulus
CompR  = 1.0/Bulk;                  %% Rock compressibility
CompW  = 2.0e-16/psia;              %% Water compressibility
CompO  = 1.0e-15/psia;              %% Oil compressibility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = Lx/double(nx);
dy = Ly/double(ny);
dz = Lz/double(nz);
G  = cartGrid([nx ny nz],[Lx Ly Lz]*meter^3);
G.nodes.coords(:, 3) = depth + G.nodes.coords(:, 3)*meter;
G.nodes.coords(:, 2) = G.nodes.coords(:, 2)*meter;
G.nodes.coords(:, 1) = G.nodes.coords(:, 1)*meter;
G  = computeGeometry(G);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_monitores = monitors(G, monitorpres, 'exp010/in/input_pres.in');
s_monitores = monitors(G, monitorsat , 'exp010/in/input_sat.in');
d_monitores = monitors(G, monitordisp, 'exp010/in/input_disp.in');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GEOLOGIC MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[K, phi, E, nu, alpha, Bulkm, CompR] = geologiaMECH(G,nini,filenx,...
        fileny,filenz,nD,K,vK,beta,rho,phibeta,phirho,E,spriggs,nu,...
        alpha,depth,KC,heter,printa,nome);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rock model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rock = makeRock(G, K, phi);
rock.alpha = alpha; clear alpha;
mK         = mean(rock.perm/(milli*darcy));
sK         = std(rock.perm/(milli*darcy));
fprintf('\n==============================================================\n\n')
fprintf('Mean k_x....: %4.1f mD \t | \t std k_x....: %4.1f mD\n',mK(1),sK(1));
fprintf('Mean k_y....: %4.1f mD \t | \t std k_y....: %4.1f mD\n',mK(2),sK(2));
fprintf('Mean k_z....: %4.1f mD \t | \t std k_z....: %4.1f mD\n',mK(3),sK(3));
fprintf('Mean phi....: %4.2f    \t | \t std phi....: %4.1e   \n',mean(rock.poro),std(rock.poro));
fprintf('Mean E......: %4.1e  \t | \t std E......: %4.1e   \n',mean(E),std(E));
fprintf('\n==============================================================\n')
clear K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if printa == 1 && (heter == 1 || heter == 2)
    plot_rock(rock.perm(:,1),G,'Yn','$\kappa_x$',color,lim,vw,1);
    base=['figuras/permKx_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    plot_rock(rock.perm(:,2),G,'Yn','$\kappa_y$',color,lim,vw,2);
    base=['figuras/permKy_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    plot_rock(rock.perm(:,3),G,'Yn','$\kappa_z$',color,lim,vw,3);
    base=['figuras/permKz_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    plot_rock(phi,G,'Yn','$\phi$',color,[0 0],vw,10)
    base=['figuras/poro_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    plot_rock(E,G,'Yn','$\mathsf{E}$',color,[0 0],vw,11)
    base=['figuras/young_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    pause(et); clf; close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RESERVOIR STATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if grav == 1
    gravity reset on;
else
    gravity off;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write intro text for each case
opt.verbose    = verbose;
writeIntroText = @(opt)(fprintf(...
    '\n*** Start new simulation\n* fluid model : %s\n* method : %s\n\n',...
    opt.fluid_model, opt.method));

%  Fluid models ===========================================================
opt.fluid_model = 'oil water';      %  Two phase oil water phases case
%opt.fluid_model = 'water';         %  water case
%opt.fluid_model = 'blackoil';      %  Three phases Black-Oil phases case

% Type of mechanical coupling =============================================
opt.method  = 'fully coupled';
opt.method  = 'fixed stress splitting';

% Boundary conditions =====================================================
opt.bc_case = 'bottom fixed';
%opt.bc_case = 'no displacement';
opt.bc_case = 'load + roller';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
writeIntroText(opt);
optvals = cellfun(@(x) opt.(x), fieldnames(opt), 'uniformoutput', false);
optlist = reshape(vertcat(fieldnames(opt)', optvals'), [], 1);
[model, states, initState, wellSols, W, dt] = runMechanicSchedule(G,rock,E,nu,...
    overburden,p_at_surface,depth,rhoR,grav,vinj,vprod,PRbhp,...
    well_r*meter,TT,Tinjstart,Tinjup,Tinjstop,nstep,CompW,CompO,CompR,...
    ndt,printa,prod_dirichlet,optlist{:});
plotMechanic(model, states, opt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
telapsed = toc(tstart);
Totaltime= seconds(telapsed);
Totaltime.Format = 'hh:mm:ss';
fprintf('\n =================================================\n')
fprintf(' Total time elapsed.......: %s\n',char(Totaltime))
fprintf(' =================================================\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if salva == 1
    save_data(initState,G,W,wellSols{1},salva,nome,0.0,0,0,0,0,...
        s_monitores,p_monitores,d_monitores);
    t = 0.0;
    for n=1:numel(states)
        t   = t + dt(1,n);
        sol = states{n};
        save_data(sol,G,W,wellSols{n},salva,nome,(t/day),n,njump,ndt,nini,...
            s_monitores,p_monitores,d_monitores);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if printa == 1
    maxp = -1.0e30;
    minp = -maxp;
    for  n = 1:numel(states)
        sol = states{n};
        pres= sol.pressure;
        maxp= max(maxp, max(pres));
        minp= min(minp, min(pres));
    end
    lim = [minp maxp]/barsa;
    lim = [80 200];
    lim = [0 0];
    npk = PandSfigures(initState,G,W,printa,vw,nome,et,0,-1,0,0,-1,lim);
    t   = 0.0; npk = 1;
    for n = 1:numel(states)
        t   = t + dt(1,n);
        sol = states{n};
        npk = PandSfigures(sol,G,W,printa,vw,nome,et,n,nprjump, ...
            (t/day),npk,ndt,lim);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if printa == 1
    uu0   = initState.uu;
    maxuu = max(uu0(:,3));
    minuu = min(uu0(:,3));
    for n = 1:numel(states)
        sol   = states{n};
        uu    = sol.uu - uu0;
        maxuu = max(maxuu, max(uu(:,3)));
        minuu = min(minuu, min(uu(:,3)));
    end
    lim = [minuu, maxuu];
    npk = VerticalDispfigures(initState.uu,G,W,printa,vw,nome,et,0,-1,0,0,-1,lim);
    t   = 0.0; npk = 1;
    for n = 1:numel(states)
        t   = t + dt(1,n);
        sol = states{n};
        npk = VerticalDispfigures(sol.uu,G,W,printa,vw,nome,et,n, ...
            nprjump,(t/day),npk,ndt,lim);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting with GUI from ad-core
if printa == 10
    mrstModule add ad-core
    plotWellSols(wellSols,cumsum(dt))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
