function thetan = loadthetas(thetan,nc,nk,iter,home,name)
    for i = 1 : nc
        for j = 1 : nk
            nome = [home name '_chain' num2str(i,'%d') '_t' ...
                num2str(j,'%d') '_' num2str(iter,'%d') '.dat'];
            t = load(nome,'-mat');
            thetan(:,i,j) = t.t;
        end
    end
    return
end