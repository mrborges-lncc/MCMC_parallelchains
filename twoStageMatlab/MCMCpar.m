close all
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath ./tools/
%addpath ../twophaseflow/mrst_borges_tools/
addpath ../twophaseflow/mrst/
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
thetan= zeros(d,NC,num_rockpar);
theta = zeros(d,NC,num_rockpar);
Y     = zeros(numel,NC,num_rockpar);
K     = zeros(numel,NC,num_rockpar);
mu    = 0.0;
s2    = 2.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start (First step) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for chain = 1 : NC
    for nk = 1 : num_rockpar
        %% Random fields generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        theta(:,chain,nk) = prop(mu,s2,d);
        Y(:,chain,nk)     = KL(T{nk},theta(:,chain,nk),KLM,numel);
    end
    %% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Simulator(Y(:,chain,:),physical_dim,fine_mesh)
end
%% MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 2 : num_trials
    fprintf('\nIteration %d\n',n)
    if n > num_select
        fprintf('\n======================================================\n')
        fprintf('\nFinished\nNumber of different selected fields achivied\n')
        fprintf('\n======================================================\n')
        break
    end
    for chain = 1 : NC
        x(:,chain) = prop(mu,s2,d)
    end
end
