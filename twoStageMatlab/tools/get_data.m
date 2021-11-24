function  get_data(G,W,wellsolution,soluc,nstep,dt,pmon)
    nw  = numel(W);
    winj= []; wbhp = [];
    for i=1:nw
        if W(i).type(1:1) == 'r' && W(i).val >= 0, winj = [winj i]; end
        if W(i).type(1:1) == 'r' && W(i).val < 0 , wbhp = [wbhp i]; end
        if W(i).type(1:1) == 'b', wbhp = [wbhp i]; end
    end
    water_cut = [];
    oil_prod  = [];
    bhpress   = [];
    t = 0;
    for n = 1 : nstep+1
        wellsol = wellsolution{n};
        sol     = soluc{n}
        for i = wbhp
            water_cut = [water_cut; t -wellsol(i).qWs*day];
            oil_prod  = [oil_prod; t  -wellsol(i).qOs*day];
        end
        for i = winj
            cells   = W(i).cells;
            [maxcoord pos] = max(G.cells.centroids(cells,3));
            bhpress = [bhpress sol.pressure(cells(pos))/mega];
        end
        satw = [t sol.s(smon,1)'];
        pres = [t sol.pressure(pmon,1)'/mega];
        t = t + dt(n)
    end
    
end

