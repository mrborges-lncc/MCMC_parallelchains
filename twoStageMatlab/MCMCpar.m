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
select_theta = zeros(d,NC,num_rockpar,num_trials);
Y     = zeros(numel,NC,num_rockpar);
K     = zeros(numel,NC,num_rockpar);
mu    = 0.0;
s2    = (2.38/sqrt(d))^2;
postn = zeros(NC,1);
post  = zeros(NC,1);
erro  = zeros(num_trials,NC);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start (First step) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for chain = 1 : NC
    for nk = 1 : num_rockpar
        %% Random fields generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        thetan(:,chain,nk) = prop(mu,s2,d);
        Y(:,chain,nk)      = KL(T{nk},thetan(:,chain,nk),KLM,numel);
    end
    %% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pres prod] = Simulator(Y(:,chain,:),physical_dim,fine_mesh);
    sample{1}   = pres;
    sample{2}   = prod;
    clear pres prod
    %% Likelihood %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    likn          = likelihood(dataref,sample,precision,num_datatype);
    postn(chain)  = posterior(likn,thetan(:,chain,:),KLM)
    ACCEPT = 1;
    if ACCEPT == 1
        for nk = 1 : num_rockpar
             select_theta(:,chain,nk,1) = thetan(:,chain,nk);
        end 
        erro(1,chain) = postn(chain);
    end
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
        ACCEPT = 0;
        %% Random fields generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for nk = 1 : num_rockpar
            theta(:,chain,nk) = thetan(:,chain,nk) + prop(mu,s2,d);
            Y(:,chain,nk)     = KL(T{nk},theta(:,chain,nk),KLM,numel);
        end
        %% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [pres prod] = Simulator(Y(:,chain,:),physical_dim,fine_mesh);
        sample{1}   = pres;
        sample{2}   = prod;
        clear pres prod
        %% Likelihood %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lik          = likelihood(dataref,sample,precision,num_datatype);
        post(chain)  = posterior(lik,theta(:,chain,:),KLM);
        %% ACCEPTANCE TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        alpha = min(1, post(chain)/postn(chain));
        coin  = rand(1,1);
        if coin < alpha
            ACCEPT = 1;
        end
        if ACCEPT == 1
            for nk = 1 : num_rockpar
                select_theta(:,chain,nk,n) = theta(:,chain,nk);
                thetan(:,chain,nk) = theta(:,chain,nk);
            end
            erro(n,chain) = post(chain);
            postn(chain)  = post(chain);
        else
            for nk = 1 : num_rockpar
                select_theta(:,chain,nk,n) = thetan(:,chain,nk);
            end
            erro(n,chain) = postn(chain);
        end
    end
end
