close all
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath ./tools/
%addpath ../twophaseflow/mrst/
addpath ~/Dropbox/mrst-2021b/
startup
%% INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[expname, prop_method, jump, nStage, num_rockpar, num_datatype, ...
    num_trials, num_select, NC] = finputbox();
[file_ref, file_sample, precision, precision_coarse] = ...
    finputbox2(nStage, num_datatype);
[physical_dim, fine_mesh, coarse_mesh, file_KL, KLM] = ...
    finputbox3(nStage, num_rockpar);
%% READ REFERENCE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataref = load_data(file_ref,num_datatype);
%% READ T matrices from KL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = load_KL(file_KL,num_rockpar,fine_mesh,KLM);
%% Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d     = KLM;
numel = fine_mesh(1) * fine_mesh(2) * fine_mesh(3);
thetan= zeros(d,num_rockpar);
theta = zeros(d,num_rockpar);
Y     = zeros(numel,num_rockpar);
mu    = 0.0;
erro  = zeros(num_trials,1);
cerro = zeros(num_trials,1);
csample = [];
TOL   = 1e-07;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1 : num_trials
    fprintf('\n==================================================\n')
    fprintf('==================================================\n')
    fprintf('\nIteration %d\n',n)
    fprintf('\n==================================================')
    fprintf('\n==================================================\n')
    %% Random fields generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for nk = 1 : num_rockpar
        theta(:,nk) = lhsnorm(0.0,1.0,d);
        Y(:,nk)     = KL(T{nk},theta(:,nk),numel);
    end
    %% Upscaling of rock parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    coarseperm = rockupscaling(Y(:,1),physical_dim, ...
        fine_mesh,coarse_mesh,'perm');
    coarseporo = rockupscaling(Y(:,2),physical_dim,...
        fine_mesh,coarse_mesh,'poro');
    %% Simulation coarse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [cpres cprod] = Simulator([coarseperm coarseporo],...
        physical_dim,coarse_mesh);
    csample{1}   = cpres;
    csample{2}   = cprod;
    clear cpres cprod
    cerro(n) = erromediorel(dataref,csample,num_datatype);
    %% Fine scale avaliation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pres prod] = Simulator([Y(:,1) Y(:,2)],...
        physical_dim,fine_mesh);
    sample{1} = pres;
    sample{2} = prod;
    clear pres prod
    erro(n) = erromediorel(dataref,sample,num_datatype);
    plot(erro(1:n), cerro(1:n),'bo','LineWidth',1,'MarkerSize',6);
    hold on
    pause(0.01)
end
E = [erro cerro];
save('erroCoarseXfine.dat','E','-ascii');
clear;
