function [newexp, expname, prop_method, jump, nStage, num_rockpar, num_datatype, ...
    num_trials, num_select, NC, freqj, prt] = finputbox()
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prompt = {'New experiment: (true or false)','Name of experiment:', 'Type of proposal (DE, RW, CW):', ...
        'Size of jump:', 'Frequence update jumping:','Number of stages:', ...
        'Number of rock parameters (permeability, porosity, ...)',...
        'Number of data types (production, pressure, saturation, etc.)',...
        'Number of trials:',...
        'Number of different selected fields:','Number of chains:',...
        'Save files? (1 == yes)'};
    dlgtitle = 'Input data';
    dims     = [1 60];
    definput = {'false','RW','RW', '0', '10', '2', '2', '2', '150', '1000', '2', '1'};
    answer   = inputdlg(prompt,dlgtitle,dims,definput);
    ans      = char(answer);
    k        = 1;
    newexp   = strtrim(strjust(ans(k+0,:)));
    expname      = strtrim(strjust(ans(k+1,:)));
    prop_method  = strtrim(strjust(ans(k+2,:)));
    jump         = str2num(ans(k+3,:));
    freqj        = str2num(ans(k+4,:));
    nStage       = str2num(ans(k+5,:));
    num_rockpar  = str2num(ans(k+6,:));
    num_datatype = str2num(ans(k+7,:));
    num_trials   = str2num(ans(k+8,:));
    num_select   = str2num(ans(k+9,:));
    NC           = str2num(ans(k+10,:));
    prt          = str2num(ans(k+11,:));
    if newexp(1:4) == 'true'
        newexp = true;
    else
        if newexp(1:4) == 'fals'
            newexp = false;
        end
    end
return

