function savethetas(thetan,chain,nk,iter,prt,home,name)
    if prt == 1
        for i = 1 : nk
            nome = [home name '_chain' num2str(chain,'%d') '_t' ...
                num2str(i,'%d') '_' num2str(iter,'%d') '.dat'];
            t = thetan(:,chain,i);
            save(nome,'t');
        end
    end
return
end