function [newexp, expname, prop_method, jump, nStage, num_rockpar, num_datatype, ...
    num_trials, NC, freqj, prt] = finputbox()
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prompt = {'Starting a new experiment? (true or false)',...
        'Name of experiment:', 'Type of proposal (options: DE, RW, CW):', ...
        'Size of jump (if == 0 => 2.38/sqrt(d)):',...
        'Jump update frequency:','Number of stages:', ...
        'Number of rock parameters (permeability, porosity, ...)',...
        'Number of data types (production, pressure, saturation, etc.)',...
        'Number of trials:','Number of chains:',...
        'Save files? (1 == yes)'};
    dlgtitle = 'Input data';
    dims     = [1 60];
    definput = {'true','nCW','CW', '0', '10', '1', '2', '2', '10', '2', '1'};
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
    NC           = str2num(ans(k+9,:));
    prt          = str2num(ans(k+10,:));
    if newexp(1:4) == 'true'
        newexp = true;
    else
        if newexp(1:4) == 'fals'
            newexp = false;
        end
    end
return

