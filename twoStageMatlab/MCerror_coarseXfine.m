close all
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath ./tools/
addpath ../twophaseflow/mrst/
%addpath ~/Dropbox/mrst-2021b/
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
read = 1;
prt  = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if read == 1
    E = load('error/erroCoarseXfine.dat');
    fitfigure(E(:,1), E(:,3), 1);
    name = 'figuras/ErrorFxErrorC1';
%     set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600',name);
    %
    fitfigure(E(:,2), E(:,4), 2);
    name = 'figuras/ErrorFxErrorC2';
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600',name);
    %
    fitfigure(E(:,1) + E(:,2), E(:,3) + E(:,4), 3);
    name = 'figuras/ErrorFxErrorC';
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600',name);
    %
    T = (E(:,1));
    T = 1+(T - mean(T))/(sqrt(var(T)/2));
    
    quisquaredhist(T,'$\mathsf{E_{f}}_{1}$');
    name = 'figuras/logErrorF1';
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600',name);
    %
    T = (E(:,2));
    T = 1+(T - mean(T))/(sqrt(var(T)/2));
    quisquaredhist(T,'$\mathsf{E_{f}}_{1}$');
    name = 'figuras/logErrorF2';
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600',name);
    %
    T = (E(:,3));
    T = 1+(T - mean(T))/(sqrt(var(T)/2));
    quisquaredhist(T,'$\mathsf{E_{c}}_{1}$');
    name = 'figuras/logErrorC1';
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600',name);
%
    T = (E(:,4));
    T = 1+(T - mean(T))/(sqrt(var(T)/2));
    quisquaredhist(T,'$\mathsf{E_{c}}_{2}$');
    name = 'figuras/logErrorC2';
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600',name);
else
    %% INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [newexp, expname, prop_method, jump, nStage, num_rockpar, num_datatype, ...
        num_trials, num_select, NC, freqj, prt] = finputbox();
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
    erro  = zeros(num_trials,num_rockpar);
    cerro = zeros(num_trials,num_rockpar);
    csample = [];
    TOL   = 1e-07;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    chain = 1;
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
        cerro(n,:) = erromedio(dataref, csample, num_datatype);
        %% Fine scale avaliation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [pres prod] = Simulator([Y(:,1) Y(:,2)],...
            physical_dim,fine_mesh);
        sample{1} = pres;
        sample{2} = prod;
        clear pres prod
        erro(n,:) = erromedio(dataref, sample, num_datatype);
    end
    E = [erro cerro];
    save('./error/erroCoarseXfine.dat','E','-ascii');
    fitfigure(erro(1:n,1), cerro(1:n,1),1);
    name = 'figuras/ErrorFxErrorC1';
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r300',name);
    fitfigure(erro(1:n,2), cerro(1:n,2),2);
    name = 'figuras/ErrorFxErrorC2';
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r300',name);
end
clear;