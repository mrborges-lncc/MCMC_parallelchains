function saverestart(n,home,name,NC,d,stage,csample,sample)
    nome = [home name '_restart_stage' num2str(stage,'%d') '_NC' ...
        num2str(NC,'%d') '_d' num2str(d,'%d') '.txt'];
    [fileID, errmsg] = fopen(nome,'w');
    fprintf(fileID,'%d %d %d %d',n,stage,d,NC);
    fclose(fileID);
end