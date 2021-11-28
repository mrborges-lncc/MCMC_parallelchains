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
thetan= zeros(d,NC,num_rockpar);
theta = zeros(d,NC,num_rockpar);
select_theta = zeros(d,NC,num_rockpar,num_trials);
Y     = zeros(numel,NC,num_rockpar);
mu    = 0.0;
erro  = zeros(num_trials,NC);
csample = [];
TOL   = 1e-07;
if abs(jump) <TOL 
    s2 = (2.38/sqrt(d))^1;
else
    s2 = jump;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start (First step) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n==================================================\n')
fprintf('==================================================\n')
fprintf('\nIteration %d\n',1)
fprintf('\n==================================================')
fprintf('\n==================================================\n')
for chain = 1 : NC
    fprintf('\n==================================================')
    fprintf('\n Iter. %d; chain: %d\n',1,chain);
    fprintf('==================================================')
    for nk = 1 : num_rockpar
        %% Random fields generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        thetan(:,chain,nk) = lhsnorm(0.0,1.0,d);
        Y(:,chain,nk)      = KL(T{nk},thetan(:,chain,nk),numel);
    end
    %% Two-Stage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nStage == 2
        %% Upscaling of rock parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        coarseperm = rockupscaling(Y(:,chain,1),physical_dim, ...
            fine_mesh,coarse_mesh,'perm');
        coarseporo = rockupscaling(Y(:,chain,2),physical_dim,...
            fine_mesh,coarse_mesh,'poro');
        %% Simulation coarse scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [cpres cprod] = Simulator([coarseperm coarseporo],...
            physical_dim,coarse_mesh);
        csamplen{1}   = cpres;
        csamplen{2}   = cprod;
        clear cpres cprod
        CACCEPT = 1;
        coarse_post_ratio = 1.0;
    end
    %% Simulation fine scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pres prod] = Simulator([Y(:,chain,1) Y(:,chain,1)],...
        physical_dim,fine_mesh);
    samplen{1}   = pres;
    samplen{2}   = prod;
    clear pres prod
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    select_theta(:,chain,:,1) = thetan(:,chain,:);
    erro(1,chain) = erromediorel(dataref,samplen,num_datatype);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 2 : num_trials
    fprintf('\n==================================================\n')
    fprintf('==================================================\n')
    fprintf('\nIteration %d\n',n)
    fprintf('\n==================================================')
    fprintf('\n==================================================\n')
    if n > num_select
        fprintf('\n==================================================\n')
        fprintf('\nFinished\nNumber of different selected fields achivied\n')
        fprintf('\n==================================================\n')
        break
    end
    for chain = 1 : NC
        ACCEPT = 0; CACCEPT = 0;
        fprintf('\n==================================================')
        fprintf('\n Iter. %d; chain: %d\n',n,chain);
        fprintf('==================================================')
        %% Random fields generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for nk = 1 : num_rockpar
            theta(:,chain,nk) = prop(thetan(:,chain,nk),mu,s2,d);
            Y(:,chain,nk)     = KL(T{nk},theta(:,chain,nk),numel);
        end
        %% Two-Stage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if nStage == 2
            %% Upscaling of rock parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            coarseperm = rockupscaling(Y(:,chain,1),physical_dim, ...
                fine_mesh,coarse_mesh,'perm');
            coarseporo = rockupscaling(Y(:,chain,2),physical_dim,...
                fine_mesh,coarse_mesh,'poro');
            %% Simulation coarse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [cpres cprod] = Simulator([coarseperm coarseporo],...
                physical_dim,coarse_mesh);
            clear csample
            csample{1}   = cpres;
            csample{2}   = cprod;
            clear cpres cprod
            %% ACCEPTANCE TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [calpha, coarse_post_ratio] = cprob_accept(dataref,csamplen,...
                csample,precision_coarse,num_datatype)
            %% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if rand(1,1) < calpha
                CACCEPT = 1;
                fprintf('###################################\n')
                fprintf('Accepted in coarse scale\nIter. %d; chain: %d\n',n,chain)
                fprintf('###################################')
            else
                CACCEPT = 0;
            end
        else
            CACCEPT = 1;
        end
        %% Fine scale avaliation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if CACCEPT == 1
            %% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [pres prod] = Simulator([Y(:,chain,1) Y(:,chain,2)],...
                physical_dim,fine_mesh);
            sample{1} = pres;
            sample{2} = prod;
            clear pres prod
            %% ACCEPTANCE TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            alpha = prob_accept(dataref,samplen,sample,precision,...
                num_datatype,coarse_post_ratio);
            %% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if rand(1,1) < alpha
                ACCEPT = 1;
            end
            if ACCEPT == 1
                samplen  = sample;
                csamplen = csample;
                clear sample csample
                select_theta(:,chain,:,n) = theta(:,chain,:);
                thetan(:,chain,:) = theta(:,chain,:);
            else
                select_theta(:,chain,:,n) = thetan(:,chain,:);
            end
        else
            select_theta(:,chain,:,n) = thetan(:,chain,:);
        end
        erro(n,chain) = erromediorel(dataref,samplen,num_datatype);
        color = 'r';
        if chain == 1
            color = 'bo';
        end
        plot(1:n, erro(1:n,chain),color,'LineWidth',3);
        hold on
        pause(0.01)
    end
end
