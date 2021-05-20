function np = save_data(sol,G,W,wellsol,printa,nome,t,n,njump,ndt,rep,smon,pmon,dmon)
    if printa == 1 && mod(n-ndt,njump) == 0 && (n > ndt || n == 0)
        nomeprod = ['./prod/prod_' nome '_' num2str(rep,'%d') '.dat'];
        nomewcut = ['./prod/wcut_' nome '_' num2str(rep,'%d') '.dat'];
        nomepres = ['./pres/presinj_' nome '_' num2str(rep,'%d') '.dat'];
        nomesmon = ['./conc/sw_' nome '_' num2str(rep,'%d') '.dat'];
        nomepmon = ['./pres/pres_' nome '_' num2str(rep,'%d') '.dat'];
        nomedmon = ['./disp/disp_' nome '_' num2str(rep,'%d') '.dat'];
        nw  = numel(W);
        winj= []; wbhp = [];
        for i=1:nw
            if W(i).type(1:1) == 'r' && W(i).val >= 0, winj = [winj i]; end
            if W(i).type(1:1) == 'r' && W(i).val < 0 , wbhp = [wbhp i]; end
            if W(i).type(1:1) == 'b', wbhp = [wbhp i]; end
        end
        water_cut = t;
        oil_prod  = t;
        bhpress   = t;
        for i = wbhp
            water_cut = [water_cut -wellsol(i).qWs*day];
            oil_prod  = [oil_prod  -wellsol(i).qOs*day];
        end
        for i = winj
            cells   = W(i).cells;
            [maxcoord pos] = max(G.cells.centroids(cells,3));
            bhpress = [bhpress sol.pressure(cells(pos))/mega];
        end
        satw = [t sol.s(smon,1)'];
        pres = [t sol.pressure(pmon,1)'/mega];
        if isempty(dmon)
            disp = [];
        else
            disp = [t sol.uu(dmon,3)'];
        end
        %% files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if n == 0
            if isfile(nomeprod)
                mode = 'w';
            else
                mode = 'a';
            end
        else
            mode = 'a';
        end
        [fileIDprod, errmsg] = fopen(nomeprod,mode);
        if length(errmsg) ~= 0
            error(errmsg);
        end
        fprintf(fileIDprod,'%8.5e ',oil_prod);
        fprintf(fileIDprod,'\n');
        fclose(fileIDprod);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if n == 0
            if isfile(nomewcut)
                mode = 'w';
            else
                mode = 'a';
            end
        else
            mode = 'a';
        end
        [fileIDwcut, errmsg] = fopen(nomewcut,mode);
        if length(errmsg) ~= 0
            error(errmsg);
        end
        fprintf(fileIDwcut,'%8.5e ',water_cut);
        fprintf(fileIDwcut,'\n');
        fclose(fileIDwcut);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if n == 0
            if isfile(nomepres)
                mode = 'w';
            else
                mode = 'a';
            end
        else
            mode = 'a';
        end
        [fileIDpres, errmsg] = fopen(nomepres,mode);
        if length(errmsg) ~= 0
            error(errmsg);
        end
        fprintf(fileIDpres,'%8.5e ',bhpress);
        fprintf(fileIDpres,'\n');
        fclose(fileIDpres);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if n == 0
            if isfile(nomesmon)
                mode = 'w';
            else
                mode = 'a';
            end
        else
            mode = 'a';
        end
        [fileIDsmon, errmsg] = fopen(nomesmon,mode);
        if length(errmsg) ~= 0
            error(errmsg);
        end
        fprintf(fileIDsmon,'%8.5e ',satw);
        fprintf(fileIDsmon,'\n');
        fclose(fileIDsmon);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if n == 0
            if isfile(nomepmon)
                mode = 'w';
            else
                mode = 'a';
            end
        else
            mode = 'a';
        end
        [fileIDpmon, errmsg] = fopen(nomepmon,mode);
        if length(errmsg) ~= 0
            error(errmsg);
        end
        fprintf(fileIDpmon,'%8.5e ',pres);
        fprintf(fileIDpmon,'\n');
        fclose(fileIDpmon);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if n == 0
            if isfile(nomedmon)
                mode = 'w';
            else
                mode = 'a';
            end
        else
            mode = 'a';
        end
        [fileIDdmon, errmsg] = fopen(nomedmon,mode);
        if length(errmsg) ~= 0
            error(errmsg);
        end
        fprintf(fileIDdmon,'%8.5e ',disp);
        fprintf(fileIDdmon,'\n');
        fclose(fileIDdmon);
    end
end

