% Single phase flow Simulator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 25/08/2021
% Based on MRST 2021
% Author: Marcio Borges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

exper  = 'exp';
tstart =  tic;

addpath ../../twophaseflow/mrst_borges_tools/
addpath ../../twophaseflow/mrst-2021a/
startup

try
   require ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui
catch %#ok<CTCH>
   mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui
end

nD    = '3D';
color = 'none';
% color = 'k';
vw  = [-35 20];
%vw  = [0 90];
verbose = true;
grav  = 1;
lim   = [0 0];
permheter  = 1;   %% if 0 => homog.; 1 => heter.
poroheter  = 1;   %% if 0 => homog.; 1 => heter.
youngheter = 1;   %% if 0 => homog.; 1 => heter.
printa= 10;
salva = 1;
monitorpres = 1;  %% if == 1 pressure monitors at some points
monitorsat  = 1;  %% if == 1 saturation monitors at some points
monitordisp = 1;  %% if == 1 displacement monitors at some points

nome  = 'ref';
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
TT        = 800.00;     %% days
Tinjstart = 100;        %% days
Tinjup    = 200;        %% days
Tinjstop  = 1700;       %% days
nstep     = 200;        %% number of time steps for pressure-velocity system
nprint    = 20;         %% Number of impressions
ndata     = 200;
ndt       = 20;
[nprint nprjump] = ajusteImpress(nprint,nstep);
[ndata  njump]   = ajusteImpress(ndata,nstep);
well_r = 0.125;         %% well radius
p_at_surface = 1.0*atm;             %% Pressure at 0m cote
depth  = 7.00e03*meter;             %% depth until the top of reservoir
rhoR   = 2.70e03*kilogram/meter^3;  %% mean density of overload rocks
PRbhp  = 0.0e03*barsa;              %% production well bh pressure
vinj   = 5.0e2/day;                 %% Injection rate
vprod  = 0.0*vinj/4;                %% Production rate
overburden= 0.0e3*atm;              %% Load (overburden)
%% Setup hidrogeological parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read Gaussian fileds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filenx  = [exper '/fields/perm_' nome '_'];
fileny  = [exper '/fields/perm_' nome '_'];
filenz  = [exper '/fields/perm_' nome '_'];
filephi = [exper '/fields/poro_' nome '_'];
fileE   = [exper '/fields/young_' nome '_'];
nini    = 0;
phibeta = 0.146;
phirho  = 0.23;
permbeta= 9.1098e-14;
permrho = 0.597;
Ebeta   = 1.0225e10;
Erho    = 0.457;
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
p_monitores = monitors(G, monitorpres, [exper '/in/input_pres.in']);
s_monitores = monitors(G, monitorsat , [exper '/in/input_sat.in']);
d_monitores = monitors(G, monitordisp, [exper '/in/input_disp.in']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GEOLOGIC MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nD  = '3D';
%% permeability
if permheter == 1
    K = load_perm(G,filenx,fileny,filenz,depth,nini,nD);
else
    K = zeros(G.cells.num,3);
end
K   = permbeta * exp(permrho * K);
%% porosity
if poroheter == 1
    phi = load_poro(G,filephi,depth,nini,nD);
else
    phi = zeros(G.cells.num,1);
end
phi = phibeta * exp(phirho * phi);
%% Young's modulus
if youngheter == 1
    E   = load_poro(G,fileE,depth,nini,nD);
else
    E   = zeros(G.cells.num,1);
end
E   = Ebeta * exp(Erho * E);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rock model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu     = 0.2;                         %% Poisson's ratio
nu_u   = 0.3;
E0     = 1.0e11;                      %% E(zero)
Ks     = E0/(3.0*(1.0 - 2*mean(nu_u))); %% Bulk modulus (grain)
Bulk   = mean(E)/(3.0*(1.0 - 2*mean(nu))); %% Bulk modulus
alpha  = 1.0 - Bulk/Ks;             %% Biot's coefficient
rock       = makeRock(G, K, phi);
rock.E     = E;
rock.nu    = nu*ones(G.cells.num,1);
rock.alpha = alpha*ones(G.cells.num,1); clear alpha;
mK         = mean(rock.perm/(milli*darcy));
sK         = std(rock.perm/(milli*darcy));
fprintf('\n==============================================================\n\n')
fprintf('Mean k_x....: %4.1f mD \t\t | std k_x....: %4.1f mD\n',mK(1),sK(1));
fprintf('Mean k_y....: %4.1f mD \t\t | std k_y....: %4.1f mD\n',mK(2),sK(2));
fprintf('Mean k_z....: %4.1f mD \t\t | std k_z....: %4.1f mD\n',mK(3),sK(3));
fprintf('Mean phi....: %4.2f    \t\t | std phi....: %4.1e   \n',mean(rock.poro),std(rock.poro));
fprintf('Mean E......: %4.1f GPa\t\t | std E......: %4.1f GPa\n',mean(E)/1e9,std(E)/1e9);
fprintf('\n==============================================================\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup material parameters for Biot and mechanics %%%%%%%%%%%%%%%%%%%%%%%
CompR  = 1.0/Bulk;                  %% Rock compressibility
CompW  = 2.0e-06/psia;              %% Water compressibility
CompO  = 7.50e-06/psia;             %% Oil compressibility
params = poroParams(mean(rock.poro), true, 'E', mean(rock.E),...
    'nu', mean(rock.nu), 'alpha', mean(rock.alpha), 'K_f', 1/CompO);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if printa == 1
    plot_rock_poro(rock.perm(:,1),G,'Yn',permbeta,permrho,...
        '$\kappa_x$',color,lim,vw,1);
    base=['../../figuras/permKx_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    plot_rock_poro(rock.perm(:,2),G,'Yn',permbeta,permrho,...
        '$\kappa_y$',color,lim,vw,2);
    base=['../../figuras/permKy_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    plot_rock_poro(rock.perm(:,3),G,'Yn',permbeta,permrho,...
        '$\kappa_z$',color,lim,vw,3);
    base=['../../figuras/permKz_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    plot_rock_poro(rock.poro,G,'Yn',phibeta,phirho,...
        '$\phi$',color,[0 0],vw,10);
    base=['../../figuras/poro_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    plot_rock_poro(rock.E,G,'Yn',Ebeta,Erho,...
        '$\mathsf{E}$',color,[0 0],vw,11);
    base=['../../figuras/young_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    pause(et); clf; close all
end
clear K phi E nu
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
%opt.method  = 'fixed stress splitting';

% Boundary conditions =====================================================
opt.bc_case = 'bottom fixed';
%opt.bc_case = 'no displacement';
opt.bc_case = 'load + roller';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
writeIntroText(opt);
optvals = cellfun(@(x) opt.(x), fieldnames(opt), 'uniformoutput', false);
optlist = reshape(vertcat(fieldnames(opt)', optvals'), [], 1);
[model, states, initState, wellSols, W, dt] = runMechanicSchedule(G,rock,...
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
    save_data(initState,G,W,wellSols{1},salva,nome,exper,0.0,0,0,0,0,...
        s_monitores,p_monitores,d_monitores);
    t = 0.0;
    for n=1:numel(states)
        t   = t + dt(1,n);
        sol = states{n};
        save_data(sol,G,W,wellSols{n},salva,nome,exper,(t/day),n,njump,ndt,nini,...
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
if printa == 1
    mrstModule add ad-core
    plotWellSols(wellSols,cumsum(dt))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
