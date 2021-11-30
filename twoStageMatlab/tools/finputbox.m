function [expname, prop_method, jump, nStage, num_rockpar, num_datatype, ...
    num_trials, num_select, NC, freqj] = finputbox()
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prompt = {'Name of experiment:', 'Type of proposal (DE, RW, CW):', ...
        'Size of jump:', 'Frequence update jumping:','Number of stages:', ...
        'Number of rock parameters (permeability, porosity, ...)',...
        'Number of data types (production, pressure, saturation, etc.)',...
        'Number of trials:',...
        'Number of different selected fields:','Number of chains:'};
    dlgtitle = 'Input data';
    dims     = [1 60];
    definput = {'RW','RW', '0', '10', '2', '2', '2', '150', '1000','2'};
    answer   = inputdlg(prompt,dlgtitle,dims,definput);
    ans      = char(answer);
    expname      = strtrim(strjust(ans(1,:)));
    prop_method  = strtrim(strjust(ans(2,:)));
    jump         = str2num(ans(3,:));
    freqj        = str2num(ans(4,:));
    nStage       = str2num(ans(5,:));
    num_rockpar  = str2num(ans(6,:));
    num_datatype = str2num(ans(7,:));
    num_trials   = str2num(ans(8,:));
    num_select   = str2num(ans(9,:));
    NC           = str2num(ans(10,:));
return

