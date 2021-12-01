function savedata(samplen,chain,nk,iter,prt,home,name)
    if prt == 1
        for i = 1 : nk
            nome = [home name '_chain' num2str(chain,'%d') ...
                '_v' num2str(i,'%d') '_' num2str(iter,'%d') '.dat'];
            t = samplen{i,chain};
            save(nome,'t');
        end
    end
return
end