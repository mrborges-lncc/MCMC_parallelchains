function [ncells1,ncells2] = quarterFS(pwell1,pwell2,G,Lx,Ly,nx,ny,flag)
    TOL = 1.0e-05;
    ncells1 = []; ncells2 = [];
    c = G.cells.centroids;
    if flag == 1
        for i = 1:G.cells.num
            aux=c(i,1:2);
            if abs(aux(1)-pwell1(1)) < TOL
                if abs(aux(2)-pwell1(2)) < TOL
                    ncells1 = [ncells1; i];
                end
            end
            if abs(aux(1)-pwell2(1)) < TOL
                if abs(aux(2)-pwell2(2)) < TOL
                    ncells2 = [ncells2; i];
                end
            end   
        end
    else
        min1 = 1.0e32; min2 = 1.0e32;
        for i = 1:G.cells.num
            aux=c(i,1:2);
            d1 = norm(aux-pwell1);
            d2 = norm(aux-pwell2);
            if d1 < min1
                min1 = d1;
                n1 = i;
            end
            if d2 < min2
                min2 = d2;
                n2 = i;
            end
        end
        p1 = c(n1,1:2);
        p2 = c(n2,1:2);
        for i = 1:G.cells.num
            aux=c(i,1:2);
            if norm(p1-aux)<TOL
                ncells1 = [ncells1; i];
            end
            if norm(p2-aux)<TOL
                ncells2 = [ncells2; i];
            end
        end
    end
end

