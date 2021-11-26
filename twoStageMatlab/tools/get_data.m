function [pres oil_prod] = get_data(G,W,wellsolution,soluc,nstep,dt,...
    pmon,ndt,ndata,njump)
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
    pres      = [];
    t = 0;
    oil_prod  = [oil_prod; t zeros(1,length(wbhp))];
    pres      = [pres; t 0.0];
    for n = 1 : nstep
        t = t + dt(n);
        if mod(n-ndt,njump) == 0 && (n > ndt || n == 0)
            %fprintf('\n%f\n',t)
            wellsol = wellsolution{n};
            sol     = soluc{n};
            prod    = [];
            for i = wbhp
                prod  = [prod  -wellsol(i).qOs*day];
            end
            oil_prod  = [oil_prod; t prod];
            for i = winj
                cells   = W(i).cells;
                [maxcoord pos] = max(G.cells.centroids(cells,3));
                bhpress = [bhpress sol.pressure(cells(pos))/mega];
            end
            pres = [pres; t sol.pressure(pmon,1)'/mega];
        end
    end
end

