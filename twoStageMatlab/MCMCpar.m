%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Metropolis algorithm
%% Author...: Marcio Borges
%% Date.....: 01/12/2021
%% Revision.: 06/06/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath ./tools/
addpath ../twophaseflow/mrst/
% addpath ~/Dropbox/mrst-2021b/
startup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
   require incomp mimetic coarsegrid upscaling
catch %#ok<CTCH>
   mrstModule add incomp mimetic coarsegrid upscaling
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
homet = './thetas/theta';
homed = './data/data';
homef = './figuras/';
homee = './error/error';
homer = './out/restart';
%% Seed control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seed = 1872;
% rng(seed)
%% INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newinput = false;
[newexp, expname, prop_method, jump, nStage, num_rockpar, ...
    num_datatype, num_trials, NC, freqj, prt, ...
    file_ref, precision, precision_coarse, data_normal,...
    physical_dim, fine_mesh, coarse_mesh, file_KL, KLM] = inputdata(newinput);
%% READ REFERENCE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cut = [1 0];
dataref = load_data(file_ref,num_datatype,cut);
[scalar, dataref] = normalizaREF(dataref,data_normal);
%% READ T matrices from KL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = load_KL(file_KL,num_rockpar,fine_mesh,KLM);
%% Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d        = KLM;
numel    = fine_mesh(1) * fine_mesh(2) * fine_mesh(3);
thetan   = zeros(d,NC,num_rockpar);
theta    = zeros(d,NC,num_rockpar);
% select_theta = zeros(d,NC,num_rockpar,num_trials);
Y        = zeros(numel,NC,num_rockpar);
mu       = 0.0;
erro     = zeros(num_trials,num_datatype,NC);
csample  = cell(num_datatype,NC);
csamplen = cell(num_datatype,NC);
sample   = cell(num_datatype,NC);
samplen  = cell(num_datatype,NC);
jump     = fjump(jump,prop_method,d,NC);
counter  = ones(NC,1);
ccounter = ones(NC,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if newexp
    %% Start (First step) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n==================================================\n')
    fprintf('==================================================\n')
    fprintf('Iteration %d',1)
    fprintf('\n==================================================')
    fprintf('\n==================================================\n')
    for chain = 1 : NC
        fprintf('\n==================================================')
        fprintf('\n Iter. %d \t chain: %d\n',1,chain);
        fprintf('==================================================')
        for nk = 1 : num_rockpar
            %% Random fields generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            thetan(:,chain,nk) = lhsnorm(0.0,1.0,d);
            Y(:,chain,nk)      = KL(T{nk},thetan(:,chain,nk),numel);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Two-Stage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if nStage == 2
            %% Upscaling of rock parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            coarseperm = rockupscaling(Y(:,chain,1),physical_dim, ...
                fine_mesh,coarse_mesh,'perm');
            coarseporo = rockupscaling(Y(:,chain,2),physical_dim,...
                fine_mesh,coarse_mesh,'poro');
            %% Simulation coarse scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [cpres cprod] = Simulator([coarseperm coarseporo],...
                physical_dim,coarse_mesh);
            csamplen{1,chain} = cpres(1+cut(1):end-cut(2),:);
            csamplen{2,chain} = cprod(1+cut(1):end-cut(2),:);
            csamplen(:,chain) = normaliza(csamplen(:,chain),scalar);
            clear cpres cprod
            coarse_post_ratio = 1.0;
        end
        %% Simulation fine scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [pres prod] = Simulator([Y(:,chain,1) Y(:,chain,2)],...
            physical_dim,fine_mesh);
        samplen{1,chain} = pres(1+cut(1):end-cut(2),:);
        samplen{2,chain} = prod(1+cut(1):end-cut(2),:);
        samplen(:,chain) = normaliza(samplen(:,chain),scalar);
        clear pres prod
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        savethetas(thetan,chain,num_rockpar,1,prt,homet,expname);
        savedata(samplen,chain,num_datatype,1,prt,homed,expname,scalar);
%         select_theta(:,chain,:,1) = thetan(:,chain,:);
        erro(1,:,chain) = erromediorel(dataref, samplen(:,chain), chain,...
            num_datatype, homee, expname, prt);
        inicio = 2;
    end
else
    [inicio, csamplen, samplen, precision_coarse, precision, ...
        ccounter, counter] = restart(homer, expname, NC, d, ...
        nStage);
    thetan = loadthetas(thetan,NC,num_rockpar,(inicio-1),homet,expname);
    fprintf('\n||||||||||||||||||||||||||||||||||||||||||||||||||\n')
    fprintf('<><><><><><><><><><><><><><><><><><><><><><><><><>\n')
    fprintf('\nRestart from iteration %d\n',inicio)
    fprintf('\n<><><><><><<><><><>><><><><><><><><><><><><><><><>')
    fprintf('\n||||||||||||||||||||||||||||||||||||||||||||||||||\n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = inicio : num_trials
    fprintf('\n==================================================\n')
    fprintf('==================================================\n')
    fprintf('Iteration %d',n)
    fprintf('\n==================================================')
    fprintf('\n==================================================\n')
%     if n > num_select
%         fprintf('\n==================================================\n')
%         fprintf('\nFinished\nNumber of different selected fields achivied\n')
%         fprintf('\n==================================================\n')
%         break
%     end
    for chain = 1 : NC
        CACCEPT = 0;
        fprintf('\n==================================================')
        fprintf('\n Iter. %d \t chain: %d\n',n,chain);
        fprintf('==================================================')
        %% Random fields generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for nk = 1 : num_rockpar
            theta(:,chain,nk) = prop(prop_method,thetan,chain,nk,jump,...
                d,NC,freqj,n);
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
            csample{1,chain} = cpres(1+cut(1):end-cut(2),:);
            csample{2,chain} = cprod(1+cut(1):end-cut(2),:);
            csample(:,chain) = normaliza(csample(:,chain),scalar);
            clear cpres cprod
            %% ACCEPTANCE TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [calpha, coarse_post_ratio] = cprob_accept(dataref,...
                csamplen(:,chain),csample(:,chain),precision_coarse,num_datatype);
            %% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if rand(1,1) < calpha
                ccounter(chain,1) = ccounter(chain,1) + 1;
                CACCEPT = 1;
                fprintf('##################################################\n')
                fprintf('Accepted in coarse scale (calpha = %f)\nIter. %d; chain: %d\n',calpha,n,chain)
                fprintf('##################################################')
            else
                CACCEPT = 0;
            end
        else
            CACCEPT = 1;
            coarse_post_ratio = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Fine scale avaliation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if CACCEPT == 1
            %% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [pres prod] = Simulator([Y(:,chain,1) Y(:,chain,2)],...
                physical_dim,fine_mesh);
            sample{1,chain} = pres(1+cut(1):end-cut(2),:);
            sample{2,chain} = prod(1+cut(1):end-cut(2),:);
            sample(:,chain) = normaliza(sample(:,chain),scalar);
            clear pres prod
            %% ACCEPTANCE TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            alpha = prob_accept(dataref,samplen(:,chain),sample(:,chain),...
                precision,num_datatype,coarse_post_ratio);
            fprintf('Alpha: %4.2f <===> coarse_post_ratio: %4.2f',alpha,coarse_post_ratio);
            %% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if rand(1,1) < alpha
                samplen(:,chain)  = sample(:,chain);
                csamplen(:,chain) = csample(:,chain);
%                 select_theta(:,chain,:,n) = theta(:,chain,:);
                thetan(:,chain,:) = theta(:,chain,:);
                counter(chain,1)  = counter(chain,1) + 1;
            else
%                 select_theta(:,chain,:,n) = thetan(:,chain,:);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
%             select_theta(:,chain,:,n) = thetan(:,chain,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        erro(n,:,chain) = erromediorel(dataref, samplen(:,chain), chain,...
            num_datatype, homee, expname, prt);
        savethetas(thetan,chain,num_rockpar,n,prt,homet,expname);
        savedata(samplen,chain,num_datatype,n,prt,homed,expname,scalar);
    end
    if mod(n,100) == 0
        for i = 1 : num_rockpar
            plot((1:n)', reshape(erro(1:n,i,:),[n NC]),'LineWidth',3);
            hold on
            pause(0.001)
        end
    end
    saverestart(n,homer,expname,NC,d,nStage,csamplen,samplen,...
        precision_coarse,precision,ccounter,counter);
    fprintf('\n*************************************************\n')
    for i = 1 : NC
        fprintf('Acceptance rate of chain %d:\n => ',i);
        if nStage == 2
            fprintf('coarse scale: %4.2f | ',100*(ccounter(i,1)/n));
%             fprintf('\nfine scale relative to coarse: %4.2f\n',...
%                 100*(counter(i,1)/ccounter(:,1)));
        end
        fprintf('fine scale: %4.2f\n', 100*(counter(i,1)/n));
    end
    fprintf('**************************************************\n')
end
clear
