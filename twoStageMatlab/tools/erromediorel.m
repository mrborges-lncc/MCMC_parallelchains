function [out] = erromediorel(cref,csample,chain,numdatatype,home,name,prt)
    out   = zeros(1,numdatatype);
    for k = 1 : numdatatype
        ref    = cref{k};
        sample = csample{k};
        out(k) = ((norm(ref(:,2:end) - sample(:,2:end))^2)/...
            (norm(ref(:,2:end))^2));
    end
    if prt == 1
        nome = [home name '_chain' num2str(chain,'%d') '.dat'];
        if exist(nome, 'file') == 0
            [fileID, errmsg]  = fopen(nome,'w');
        else
            [fileID, errmsg]  = fopen(nome,'a');
        end
        for k = 1 : numdatatype
            fprintf(fileID , '%f ' ,out(k));
        end
        fprintf(fileID , '\n' );
        fclose(fileID);
    end
    return
end