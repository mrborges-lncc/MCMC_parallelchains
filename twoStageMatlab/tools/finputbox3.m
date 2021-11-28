function [physical_dim, fine_mesh, coarse_mesh, file_KL, KLM] = ...
    finputbox3(nStage, num_rockpar)
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dlgtitle = 'Input information to prior generation';
    dims     = [1 80];
    pk = 1;
    dk = 1;
    %% Grid information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prompt{pk}   = 'Domain dimensions (Lx, Ly, Lz):';
    definput{dk} = '500.0 500.0 20.0';
    %% fine mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pk = pk + 1; dk = dk + 1;
    prompt{pk} = 'Fine mesh (nx, ny, nz):';
    definput{dk} = '100 100 1';
    %% coarse mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nStage == 2
        pk = pk + 1; dk = dk + 1;
        prompt{pk} = 'Coarse mesh (cnx, cny, cnz):';
        definput{dk} = '25 25 1';
    end
    %% input data files names
    name = 'Name of KL matrix of prior ';
    for i = 1:num_rockpar
        pk = pk + 1;
        nome   = [name num2str(i,'%d') ':'];
        prompt{pk} = nome;
    end
    for i = 1:num_rockpar
        dk = dk + 1;
        aux = 'prod';
        if i == 1, aux = 'pres'; end
        nome   = ['../gera_KL/MATLAB/out/avet1_500x500x20_100x100x1_l50x50x10_M10000.bin'];
        definput{dk} = nome;
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pk = pk + 1; dk = dk + 1;
    prompt{pk}   = 'Stochastic dimention of KL:';
    definput{dk} = '10000';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    answer   = inputdlg(prompt,dlgtitle,dims,definput);
    ans      = char(answer);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k = 0;
    physical_dim = str2num(strtrim(strjust(ans(k+1,:))));
    fine_mesh    = str2num(strtrim(strjust(ans(k+2,:))));
    k = k + 2;
    coarse_mesh = [];
    if nStage == 2
        coarse_mesh  = str2num(strtrim(strjust(ans(k+1,:))));
        k = k + 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    file_KL = []; % Files with reference data
    for i = 1 : num_rockpar
        k = k + 1;
        file_KL = [file_KL; ans(k,:)];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k = k + 1;
    KLM = str2num(strtrim(strjust(ans(k,:))));
return

