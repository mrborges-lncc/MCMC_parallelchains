function [pres oilprod] = Simulator(Y, physicaldim, meshg)
% Single phase flow Simulator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 25/03/2020
% Based on MRST 2019
% Author: Marcio Borges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tstart =  tic;
TOL = 1e-07;
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
%vw  = [0 90];
grav= 1;
lim = [0 0];
%lim = 1.0e+02 * [0.85   2.0];
heter    = 1;   %% if 0 => homog.; 1 => heter.; 2 => read *.mat
phiheter = 1;
printa   = 10;
salva    = 1;        %% if == 1 save well informations
monitorpres = 1;  %% if == 1 pressure monitors at some points
monitorsat  = 1;  %% if == 1 saturation monitors at some points
nome  = 'amostra';
% nome  = 'ref';
et    = 0;
verb  = false;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx  = physicaldim(1);
Ly  = physicaldim(2);
Lz  = physicaldim(3);
nx  = meshg(1);
ny  = meshg(2);
nz  = meshg(3);
well_r= 0.125;          %% well radius
TT    = 300.0;            %% days
nstep = 300;            %% number of time steps for pressure-velocity system
nprint= 25;             %% Number of impressions
ndata = 150;             %% Number of impressions of data
ndt   = 10;
[nprint nprjump] = ajusteImpress(nprint,nstep);
[ndata njump] = ajusteImpress(ndata,nstep);
PRbhp = 0.0;            %% production well pressure
vinj  = 0.5e3/day;      %% Injection rate
patm  = 1.0*atm;        %% Pressure at 0m cote
depth = 1.0e03*meter;   %% depth until the top of reservoir
rhoR  = 2.70e03*kilogram/meter^3;  %% mean density of overload rocks
overburden= 00.0*atm;   %% Load (overburden)
fatk  = milli() * darcy();      %% Factor to permeability
phibeta = 0.146;
phirho  = 0.23;
permbeta= 9.1098e-14;
permrho = 0.597;
permrho = 1.0;
Ebeta   = 1.0225e10;
Erho    = 0.457;
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
%% Read Gaussian fileds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = permbeta * exp(permrho * Y(:,1));
K = [K K K];
phi = phibeta * exp(phirho * Y(:,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rock model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rock = makeRock(G, K, phi);
pv   = poreVolume(G, rock);
mK   = mean(rock.perm/(milli*darcy));
sK   = std(rock.perm/(milli*darcy));
% fprintf('\n==============================================================\n\n')
% fprintf('Mean k_x....: %4.1f mD \t | \t std k_x....: %4.1f mD\n',mK(1),sK(1));
% fprintf('Mean k_y....: %4.1f mD \t | \t std k_y....: %4.1f mD\n',mK(2),sK(2));
% fprintf('Mean k_z....: %4.1f mD \t | \t std k_z....: %4.1f mD\n',mK(3),sK(3));
% fprintf('Mean phi....: %4.2f    \t | \t std phi....: %4.2f   \n',mean(rock.poro),std(rock.poro));
% fprintf('\n==============================================================\n')
clear K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if printa == 1 && heter == 1 
    plot_rock(reverseKlog(rock.perm(:,1),permbeta,permrho),G,'Yn','$\kappa_x$',color,lim,vw,1);
    base=['./figuras/permKx_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    plot_rock(reverseKlog(rock.perm(:,2),permbeta,permrho),G,'Yn','$\kappa_y$',color,lim,vw,2);
    base=['./figuras/permKy_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    plot_rock(reverseKlog(rock.perm(:,3),permbeta,permrho),G,'Yn','$\kappa_z$',color,lim,vw,3);
    base=['./figuras/permKz_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    plot_rock_poro(rock.poro,G,'Yn',1,1,'$\phi$',color,[0 0],vw,14);
    base=['./figuras/phi_' nome];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600', base);
    pause(et); clf; close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FLUID PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fluid = initSimpleFluid('mu' , [ 1, 10]*centi*poise, ...
fluid = initCoreyFluid('mu' , [ 1, 20]*centi*poise, ...
    'rho', [1000, 860]*kilogram/meter^3, ...
    'n'  , [ 2, 2], ...
    'sr' , [0.0, 0.0], ...
    'kwm', [1.0, 1.0]);
% fluid = initSimpleFluid('mu' , [ 1, 1]*centi*poise, ...
%     'rho', [1014, 1014]*kilogram/meter^3, ...
%     'n'  , [ 1, 1]);
[mu,rho] = fluid.properties();
if printa == 1
    s=linspace(0,1,100)';
    kr = fluid.relperm(s);
    two2Dplot(s,[kr(:,1),kr(:,2)],'$s_w$','$kr_\alpha$','$kr_w$','$kr_o$',1);
    pause(et); clf; close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RESERVOIR STATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if grav == 1
    gravity reset on;
    p0 = 0.0;
else
    gravity off;
    p0 = 0.0;                   % initial pressure
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute overburden
if overburden < TOL
    g = norm(gravity);
    overburden = patm + depth * rhoR * g;
end
p_at_topR = overburden;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E     = 10 * giga * Pascal; %% Young's module
nu    = 0.3;               %% Poisson's ratio
alpha = 1;                 %% Biot's coefficient
CompO = 1.0e-15/psia;      %% Oil compressibility
% params = poroParams(mean(rock.poro), true, 'E', mean(E),...
%     'nu', mean(nu), 'alpha', mean(alpha), 'K_f', 1/CompO);
params.B = 1.0;
ptop = overburden * params.B;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute BHPressure
if PRbhp < TOL
    g     = norm(gravity);
    equil = ode23(@(z,p) g.* rho(2), ...
        [0, max(G.nodes.coords(:,3))], patm);
    BHPressure = reshape(deval(equil, max(G.cells.centroids(:,3))),...
        [], 1);  clear equil
else
    BHPressure = PRbhp;
end
% fprintf('\n==============================================================\n')
% fprintf('\n==============================================================\n\n')
% fprintf('Load at the top of the reservoir.........: %4.1f MPa\n',overburden/mega);
% fprintf('Fluid pressure at the top of reservoir...: %4.1f MPa\n',ptop/mega)
% fprintf('BHP (Bottom hole pressure)...............: %4.1f MPa\n',BHPressure/mega);
% fprintf('\n==============================================================\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WELLS five-spot model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[inj_cells,prod1_cells,prod2_cells,prod3_cells,prod4_cells] = ...
    fivespot_wells(G,Lx,Ly,nx,ny);
refdepth = min(G.cells.centroids(:, 3)); % for example...
W = addWell([], G, rock, inj_cells, 'Type', 'rate', 'Comp_i', [1 0],...
     'Val', vinj, 'Radius', well_r*meter, 'Dir', 'z', 'sign',1,...
     'name', '$w_{inj}$', 'refDepth', refdepth); %m^3/sec
% W = addWell([], G, rock, inj_cells, 'Type', 'bhp', 'Comp_i', [1 0],...
%      'Val', BHPressure, 'Radius', well_r*meter, 'Dir', 'z',...
%      'sign',1,'name', '$Winj$'); % Pascal
W = addWell(W, G, rock, prod1_cells, 'Type', 'bhp', 'Comp_i', [0 1],...
     'Val', BHPressure, 'Radius', well_r*meter, 'Dir', 'z', ...
     'name', '$w1$', 'refDepth', refdepth); %Pascal
W = addWell(W, G, rock, prod2_cells, 'Type', 'bhp', 'Comp_i', [0 1],...
     'Val', BHPressure, 'Radius', well_r*meter, 'Dir', 'z', ...
     'name', '$w2$', 'refDepth', refdepth);
W = addWell(W, G, rock, prod3_cells, 'Type', 'bhp', 'Comp_i', [0 1],...
     'Val', BHPressure, 'Radius', well_r*meter, 'Dir', 'z', ...
     'name', '$w3$', 'refDepth', refdepth);
W = addWell(W, G, rock, prod4_cells, 'Type', 'bhp', 'Comp_i', [0 1],...
     'Val', BHPressure, 'Radius', well_r*meter, 'Dir', 'z', ...
     'name', '$w4$', 'refDepth', refdepth);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FLUID SOURCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%src = addSource([] , cells, rates);
%src = addSource(src, cells, rates, 'sat', sat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bc = []; % all outer faces of a grid are assumed as no-flow boundaries
         % unless other conditions are specified explicitly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONSTRUCTION OF LINEAR SYSTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Incompressible Two-Point Pressure Solver
%% Transmissibility
hT = computeTrans(G, rock);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute an initial single-phase pressure solution, from which we estimate
% the final time that corresponds to 1.2 PVI if this flow field remains
% unchanged. With quadratic relperm curves and equal viscosities, the
% multiphase displacement front will propagate at a speed of
% a=1/(2(sqrt(2)-1)) relative to the total velocity.
%% Impose vertical equilibrium
g   = norm(gravity);
pf  = @(z) rho(2) * g * z;
p0  = ptop + pf(abs(depth-G.cells.centroids(:,3)));
s0  = [0.0 1.0];
sol = initState(G, W, p0, s0);
npk = PandSfigures(sol,G,W,printa,vw,nome,et,0,-1,0,0,-1,lim);
T   = 1.2*sum(pv)/sol.wellSol(1).flux;
a   = 1/(2*(sqrt(2)-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize number of time intervals, cell array to hold well solutions,
% and array to hold the oil in place
TT = TT*day;
dt = TT/double(nstep)*ones(1,nstep);
dt = [dt(1).*sort(repmat(2.^-[1:ndt ndt],1,1)) dt(2:end)];
nstep    = numel(dt);
wellSols = cell(nstep+1,1);  
%wellSols{1} = getWellSol(W, sol, fluid);
%solucao{1}  = sol; 
t = 0;
% hwb = waitbar(t,'Simulation ..');
for n=1:nstep
    t = t + dt(n);
%    fprintf(1,'Time step %d/%d <=> %5.4f days\n',n,nstep,(t/day));
    sol  = incompTPFA(sol, G, hT, fluid, 'wells', W, 'verbose', verb);
    sol  = explicitTransport(sol, G, dt(n), rock, fluid,...
        'wells', W, 'verbose', verb, 'dt_factor', 0.75);
    npk  = PandSfigures(sol,G,W,printa,vw,nome,et,n,nprjump,(t/day),npk,ndt,lim);
    wellSols{n} = getWellSol(W, sol, fluid);
    solucao{n} = sol;
end
p_monitores = monitors(G, 1, [Lx/2 Lx/2 (depth+dz*0.5)]);
[pres oilprod] = get_data(G,W,wellSols,solucao,nstep,dt/day,...
    p_monitores,ndt,ndata,njump);


% close(hwb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
telapsed = toc(tstart);
Totaltime= seconds(telapsed);
Totaltime.Format = 'hh:mm:ss';
fprintf('\n =================================================\n')
fprintf(' Total time elapsed.......: %s\n',char(Totaltime))
fprintf(' =================================================\n')%clear
