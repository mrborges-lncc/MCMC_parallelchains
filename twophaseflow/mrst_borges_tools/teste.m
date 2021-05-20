function [perm] = teste(geo,g,k,dx,dy,dz)
    crock = makeRock(geo, k, 1.0);
    rock  = makeRock(g, [0.0 0.0 0.0], 1.0);
    for j=1:g.cells.num
        c = g.cells.centroids(j,:);
        for i=1:geo.cells.num
            cr = geo.cells.centroids(i,:) + [dx/2 dy/2 dz/2];
            cl = geo.cells.centroids(i,:) - [dx/2 dy/2 dz/2];
            if(c(1) < cr(1) && c(1) >= cl(1))
                if(c(2) < cr(2) && c(2) >= cl(2))
                    if(c(3) < cr(3) && c(3) >= cl(3))
                        rock.perm(j,:) = crock.perm(i,:);
                        break
                    end
                end
            end
        end
    end
    perm = rock.perm;
end

