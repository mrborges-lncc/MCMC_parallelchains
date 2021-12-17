function savedata(samplen,chain,nk,iter,prt,home,name,scal)
    if prt == 1
        for i = 1 : nk
            nome = [home name '_chain' num2str(chain,'%d') ...
                '_v' num2str(i,'%d') '_' num2str(iter,'%d') '.dat'];
            t = samplen{i,chain};
            t(:,2:end) = t(:,2:end) * (scal(i,2) - scal(i,1)) + scal(i,1);
            save(nome,'t');
        end
    end
return
end