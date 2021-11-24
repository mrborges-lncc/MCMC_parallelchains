function [file_ref,file_sample,precision,precision_coarse] = ...
    finputbox2(nStage, num_datatype)
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dlgtitle = 'Input data set to likelihood';
    dims     = [1 60];
    definput = {};
    prompt = {};
    pk = 0;
    dk = 0;
    %% input data files names
    name = 'Name of input data file (Reference) ';
    for i = 1:num_datatype
        pk = pk + 1;
        nome   = [name num2str(i,'%d') ':'];
        prompt{pk} = nome;
    end
    for i = 1:num_datatype
        dk = dk + 1;
        aux = 'prod';
        if i == 1, aux = 'pres'; end
        nome   = ['../twophaseflow/exp/' aux '/' aux '_ref_0.dat'];
        definput{dk} = nome;
    end    
    %% input data files names
    name = 'Name of input data file (Samples) ';
    for i = 1 : num_datatype
        pk = pk + 1;
        nome   = [name num2str(i,'%d') ':'];
        prompt{pk} = nome;
    end
    for i = 1 : num_datatype
        dk = dk + 1;
        aux = 'prod';
        if i == 1, aux = 'pres'; end
        nome   = ['../twophaseflow/exp/' aux '/' aux '_amostra_0.dat'];
        definput{dk} = nome;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    name = 'Precision for likelihood (\sigma^2) data set n. (';
    if nStage == 1
        for j = 1 : num_datatype
            nome = [name num2str(j,'%d') '):'];
            pk = pk + 1;
            prompt{pk} = nome;
            dk = dk + 1;
            definput{dk} = '2.e-4';
        end
    end
    if nStage == 2
        name = 'Precision for likelihood (\sigma^2) data set n. (';
        name = 'Precision for likelihood (\sigma^2) data set n. (';
        for j = 1 : num_datatype
            nome = [name num2str(j,'%d') ') finescale:'];
            pk = pk + 1;
            prompt{pk} = nome;
            dk = dk + 1;
            definput{dk} = '2.e-4';
            nome = [name num2str(j,'%d') ') coarsescale:'];
            pk = pk + 1;
            prompt{pk} = nome;
            dk = dk + 1;
            definput{dk} = '8.e-4';
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    answer   = inputdlg(prompt,dlgtitle,dims,definput);
    ans      = char(answer);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k = 0;
    file_ref = []; % Files with reference data
    for i = 1 : num_datatype
        k = k + 1;
        file_ref = [file_ref; ans(k,:)];
    end
    file_sample = []; % Files with reference data
    for i = 1 : num_datatype
        k = k + 1;
        file_sample = [file_sample; ans(k,:)];
    end
    precision = [];
    precision_coarse = [];
    if nStage == 1
        for i = 1 : num_datatype
            k = k + 1;
            precision = [precision; str2num(ans(k,:))];
        end
    end
    if nStage == 2
        for i = 1 : num_datatype
            k = k + 1;
            precision = [precision; str2num(ans(k,:))];
            k = k + 1;
            precision_coarse = [precision_coarse; str2num(ans(k,:))];
        end
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return

