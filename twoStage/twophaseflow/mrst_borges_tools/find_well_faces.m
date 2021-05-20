function [wfaces bfaces] = find_well_faces(G,raio,pos,r)
    TOL = 1e-8;
    wfaces = [];
    bfaces= [];
    f = boundaryFaces(G); %% all exterior faces
    f = f(abs(G.faces.normals(f,3))<TOL); %% exterior vertical faces
    fsz = numel(f);
    TOL = 0.5 * (raio-r);
    for i = 1:fsz
        c = G.faces.centroids(f(i),1:2).';
        d = norm(pos-c,2);
        if d < TOL
            wfaces = [wfaces f(i)];
        else
            bfaces = [bfaces f(i)];
        end
    end
end

