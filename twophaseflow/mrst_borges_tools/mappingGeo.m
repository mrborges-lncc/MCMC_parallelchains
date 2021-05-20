function [perm] = mappingGeo(geo,g,k,dx,dy,dz)
    crock = makeRock(geo, k, 1.0);
    rock  = makeRock(g, [0.0 0.0 0.0], 1.0);
    cells = 1:g.cells.num;
    ncell = [];
    for i=1:geo.cells.num
        cr = geo.cells.centroids(i,:) + [dx/2 dy/2 dz/2];
        cl = geo.cells.centroids(i,:) - [dx/2 dy/2 dz/2];
        elems = setdiff(cells,ncell);
        for j=1:length(elems)
            k = elems(j);
            c = g.cells.centroids(k,:);
            if(c(1) < cr(1) && c(1) >= cl(1))
                if(c(2) < cr(2) && c(2) >= cl(2))
                    if(c(3) < cr(3) && c(3) >= cl(3))
                        rock.perm(k,:) = crock.perm(i,:);
                        ncell = [ncell k];
                    end
                end
            end
        end
    end
    perm = rock.perm;
end

