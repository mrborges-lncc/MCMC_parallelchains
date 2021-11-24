function [points] = monitors(G, monitor, local);
    points = [];
    if monitor == 1
        coords = local;
        nc = size(coords,1);
        for k=1:nc
            p = coords(k,:);
            d = 1e20;
            for i=1:G.cells.num
                if(norm(p - G.cells.centroids(i,:)) < d)
                    ponto = i;
                    d = norm(p - G.cells.centroids(i,:));
                end
            end
            points = [points ponto];
        end
    end
end

