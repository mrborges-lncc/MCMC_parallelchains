function [n, csample, sample, cprec, prec, ccount, count] = ...
    restart(home,name,NC,d,stage)
    TOL  = 1.0e-7;
    nome = [home name '_stage' num2str(stage,'%d') '_NC' ...
        num2str(NC,'%d') '_d' num2str(d,'%d') '.txt'];
    if exist(nome, 'file') == 0
        error('\nError file %s does not exist.\n',nome);
    end
    [fileID, errmsg] = fopen(nome,'r');
    A = fscanf(fileID,'%g');
    fclose(fileID);
    l = (length(A) - 4) / 2;
    n = A(1);
    if abs(A(2) - stage) > TOL, error('Wrong number of stage'); end
    if abs(A(3) - d)     > TOL, error('Wrong stochastic dimension'); end
    if abs(A(4) - NC)    > TOL, error('Wrong number of chains'); end
    cprec = [];
    prec = [];
    k = 4;
    for i = 1 : l
        k = k + 1;
        cprec = [cprec; A(k)];
        k = k + 1;
        prec  = [prec; A(k)];
    end
    %
    nome = [home name '_stage' num2str(stage,'%d') '_NC' ...
        num2str(NC,'%d') '_d' num2str(d,'%d') ...
        '_samplen.mat'];
    aux = load(nome);
    sample = aux.sample
    nome = [home name '_stage' num2str(stage,'%d') '_NC' ...
        num2str(NC,'%d') '_d' num2str(d,'%d') ...
        '_csamplen.mat'];
    aux = load(nome);
    csample = aux.csample;
    %
    nome = [home name '_stage' num2str(stage,'%d') '_NC' ...
        num2str(NC,'%d') '_d' num2str(d,'%d') ...
        '_ccounter.mat'];
    aux = load(nome);
    ccount = aux.ccounter;
    nome = [home name '_stage' num2str(stage,'%d') '_NC' ...
        num2str(NC,'%d') '_d' num2str(d,'%d') ...
        '_counter.mat'];
    aux = load(nome);
    count = aux.counter;
end