function [bc,faces,pressure] = bc_fivespotDirichlet(G,BHP,rho,g,init_sat)
    TOL= 1.0e-07;
    faces = [];
    % well 1
    f  = find(abs(G.faces.centroids(:,1) - min(G.nodes.coords(:,1))) < TOL);
    d  = min(G.faces.centroids(f,2));
    y  = find(abs(d - G.faces.centroids(f,2)) < TOL);
    faces = [faces; f(y)];
    f  = find(abs(G.faces.centroids(:,2) - min(G.nodes.coords(:,2))) < TOL);
    d  = min(G.faces.centroids(f,1));
    y  = find(abs(d - G.faces.centroids(f,1)) < TOL);
    faces = [faces; f(y)];
    % well 2
    f  = find(abs(G.faces.centroids(:,1) - max(G.nodes.coords(:,1))) < TOL);
    d  = min(G.faces.centroids(f,2));
    y  = find(abs(d - G.faces.centroids(f,2)) < TOL);
    faces = [faces; f(y)];
    f  = find(abs(G.faces.centroids(:,2) - min(G.nodes.coords(:,2))) < TOL);
    d  = max(G.faces.centroids(f,1));
    y  = find(abs(d - G.faces.centroids(f,1)) < TOL);
    faces = [faces; f(y)];
    % well 3
    f  = find(abs(G.faces.centroids(:,1) - max(G.nodes.coords(:,1))) < TOL);
    d  = max(G.faces.centroids(f,2));
    y  = find(abs(d - G.faces.centroids(f,2)) < TOL);
    faces = [faces; f(y)];
    f  = find(abs(G.faces.centroids(:,2) - max(G.nodes.coords(:,2))) < TOL);
    d  = max(G.faces.centroids(f,1));
    y  = find(abs(d - G.faces.centroids(f,1)) < TOL);
    faces = [faces; f(y)];
    % well 4
    f  = find(abs(G.faces.centroids(:,1) - min(G.nodes.coords(:,1))) < TOL);
    d  = max(G.faces.centroids(f,2));
    y  = find(abs(d - G.faces.centroids(f,2)) < TOL);
    faces = [faces; f(y)];
    f  = find(abs(G.faces.centroids(:,2) - max(G.nodes.coords(:,2))) < TOL);
    d  = min(G.faces.centroids(f,1));
    y  = find(abs(d - G.faces.centroids(f,1)) < TOL);
    faces = [faces; f(y)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    equil = ode23(@(z,p) g.* rho(p), ...
        [max(G.faces.centroids(faces,3)), 0], BHP);
    pressure = reshape(deval(equil, G.faces.centroids(faces,3)), ...
        [], 1);  clear equil
    bc.face  = faces;
    bc.type  = repmat({'pressure'}, 1, numel(faces));
    bc.value = pressure;  
    bc.sat   = ones(numel(faces), 1) * init_sat;
end

