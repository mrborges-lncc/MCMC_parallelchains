function saverestart(n,home,name,NC,d,stage,csample,sample,cprec,prec)
    nome = [home name '_stage' num2str(stage,'%d') '_NC' ...
        num2str(NC,'%d') '_d' num2str(d,'%d') '.txt'];
    [fileID, errmsg] = fopen(nome,'w');
    fprintf(fileID,'%d %d %d %d',n+1,stage,d,NC);
    for i = 1 : length(prec)
        fprintf(fileID,' %g %g',cprec(i),prec(i));
    end
    fclose(fileID);
    nome = [home name '_stage' num2str(stage,'%d') '_NC' ...
        num2str(NC,'%d') '_d' num2str(d,'%d') ...
        '_samplen.mat'];
    save(nome,'sample');
    nome = [home name '_stage' num2str(stage,'%d') '_NC' ...
        num2str(NC,'%d') '_d' num2str(d,'%d') ...
        '_csamplen.mat'];
    save(nome,'csample');
end