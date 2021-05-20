function [top_faces, bottom_nodes, top_nodes, right_nodes, left_nodes,...
    front_nodes, back_nodes, fr_corner_nodes, fl_corner_nodes, ...
    br_corner_nodes, bl_corner_nodes, tr_corner_nodes, tl_corner_nodes, ...
    tf_corner_nodes, tb_corner_nodes, bor_corner_nodes, ...
    bol_corner_nodes, bof_corner_nodes, bob_corner_nodes, ...
    trf_corner_nodes, trb_corner_nodes, tlf_corner_nodes, ...
    tlb_corner_nodes, borf_corner_nodes, borb_corner_nodes, ...
    bolf_corner_nodes, bolb_corner_nodes, bottom_innermost] = ...
    find_boudary_nodes(G)
    TOL = 1.0e-07;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find bottom nodes
    bottom_nodes = find(abs(G.nodes.coords(:,3) - max(G.nodes.coords(:,3))) < TOL);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find top nodes
    top_nodes    = find(abs(G.nodes.coords(:,3) - min(G.nodes.coords(:,3))) < TOL);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find top faces
    top_faces    = find(abs(G.faces.centroids(:,3) - min(G.faces.centroids(:,3))) < TOL);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find right nodes
    right_nodes  = find(abs(G.nodes.coords(:,1) - max(G.nodes.coords(:,1))) < TOL);
    right_nodes  = setdiff(right_nodes, top_nodes);
    right_nodes  = setdiff(right_nodes, bottom_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find left nodes
    left_nodes   = find(abs(G.nodes.coords(:,1) - min(G.nodes.coords(:,1))) < TOL);
    left_nodes   = setdiff(left_nodes, top_nodes);
    left_nodes   = setdiff(left_nodes, bottom_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find front nodes
    front_nodes  = find(abs(G.nodes.coords(:,2) - min(G.nodes.coords(:,2))) < TOL);
    front_nodes  = setdiff(front_nodes, top_nodes);
    front_nodes  = setdiff(front_nodes, bottom_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find back nodes
    back_nodes  = find(abs(G.nodes.coords(:,2) - max(G.nodes.coords(:,2))) < TOL);
    back_nodes  = setdiff(back_nodes, top_nodes);
    back_nodes  = setdiff(back_nodes, bottom_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner FR
    fr_corner_nodes = intersect(front_nodes,right_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner FL
    fl_corner_nodes = intersect(front_nodes,left_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner BR
    br_corner_nodes = intersect(back_nodes,right_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner BL
    bl_corner_nodes = intersect(back_nodes,left_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner TR
    tr_corner_nodes = find(abs(G.nodes.coords(top_nodes,1) - ...
        max(G.nodes.coords(top_nodes,1))) < TOL);
    tr_corner_nodes = top_nodes(tr_corner_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner TL
    tl_corner_nodes = find(abs(G.nodes.coords(top_nodes,1) - ...
        min(G.nodes.coords(top_nodes,1))) < TOL);
    tl_corner_nodes = top_nodes(tl_corner_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner TB
    tb_corner_nodes = find(abs(G.nodes.coords(top_nodes,2) - ...
        max(G.nodes.coords(top_nodes,2))) < TOL);
    tb_corner_nodes = top_nodes(tb_corner_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner TF
    tf_corner_nodes = find(abs(G.nodes.coords(top_nodes,2) - ...
        min(G.nodes.coords(top_nodes,2))) < TOL);
    tf_corner_nodes = top_nodes(tf_corner_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner BoR
    bor_corner_nodes = find(abs(G.nodes.coords(bottom_nodes,1) - ...
        max(G.nodes.coords(bottom_nodes,1))) < TOL);
    bor_corner_nodes = bottom_nodes(bor_corner_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner BoL
    bol_corner_nodes = find(abs(G.nodes.coords(bottom_nodes,1) - ...
        min(G.nodes.coords(bottom_nodes,1))) < TOL);
    bol_corner_nodes = bottom_nodes(bol_corner_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner BoB
    bob_corner_nodes = find(abs(G.nodes.coords(bottom_nodes,2) - ...
        max(G.nodes.coords(bottom_nodes,2))) < TOL);
    bob_corner_nodes = bottom_nodes(bob_corner_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner BoF
    bof_corner_nodes = find(abs(G.nodes.coords(bottom_nodes,2) - ...
        min(G.nodes.coords(bottom_nodes,2))) < TOL);
    bof_corner_nodes = bottom_nodes(bof_corner_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner TRB
    trb_corner_nodes = intersect(tr_corner_nodes,tb_corner_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner TRF
    trf_corner_nodes = intersect(tr_corner_nodes,tf_corner_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner TLB
    tlb_corner_nodes = intersect(tl_corner_nodes,tb_corner_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner TLF
    tlf_corner_nodes = intersect(tl_corner_nodes,tf_corner_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner BoRB
    borb_corner_nodes = intersect(bor_corner_nodes,bob_corner_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner BoRF
    borf_corner_nodes = intersect(bor_corner_nodes,bof_corner_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner BoLB
    bolb_corner_nodes = intersect(bol_corner_nodes,bob_corner_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % corner BoLF
    bolf_corner_nodes = intersect(bol_corner_nodes,bof_corner_nodes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % adjust (- corners)
    right_nodes = setdiff(right_nodes,fr_corner_nodes);
    right_nodes = setdiff(right_nodes,br_corner_nodes);
    left_nodes  = setdiff(left_nodes,fl_corner_nodes);
    left_nodes  = setdiff(left_nodes,bl_corner_nodes);
    front_nodes = setdiff(front_nodes,fr_corner_nodes);
    front_nodes = setdiff(front_nodes,fl_corner_nodes);
    back_nodes  = setdiff(back_nodes,br_corner_nodes);
    back_nodes  = setdiff(back_nodes,bl_corner_nodes);
    top_nodes   = setdiff(top_nodes,tr_corner_nodes);
    top_nodes   = setdiff(top_nodes,tl_corner_nodes);
    top_nodes   = setdiff(top_nodes,tb_corner_nodes);
    top_nodes   = setdiff(top_nodes,tf_corner_nodes);
    bottom_nodes= setdiff(bottom_nodes,bor_corner_nodes);
    bottom_nodes= setdiff(bottom_nodes,bol_corner_nodes);
    bottom_nodes= setdiff(bottom_nodes,bob_corner_nodes);
    bottom_nodes= setdiff(bottom_nodes,bof_corner_nodes);
    tr_corner_nodes = setdiff(tr_corner_nodes,trb_corner_nodes);
    tr_corner_nodes = setdiff(tr_corner_nodes,trf_corner_nodes);
    tb_corner_nodes = setdiff(tb_corner_nodes,trb_corner_nodes);
    tb_corner_nodes = setdiff(tb_corner_nodes,tlb_corner_nodes);
    tl_corner_nodes = setdiff(tl_corner_nodes,tlb_corner_nodes);
    tl_corner_nodes = setdiff(tl_corner_nodes,tlf_corner_nodes);
    tf_corner_nodes = setdiff(tf_corner_nodes,trf_corner_nodes);
    tf_corner_nodes = setdiff(tf_corner_nodes,tlf_corner_nodes);
    bor_corner_nodes = setdiff(bor_corner_nodes,borb_corner_nodes);
    bor_corner_nodes = setdiff(bor_corner_nodes,borf_corner_nodes);
    bob_corner_nodes = setdiff(bob_corner_nodes,borb_corner_nodes);
    bob_corner_nodes = setdiff(bob_corner_nodes,bolb_corner_nodes);
    bol_corner_nodes = setdiff(bol_corner_nodes,bolb_corner_nodes);
    bol_corner_nodes = setdiff(bol_corner_nodes,bolf_corner_nodes);
    bof_corner_nodes = setdiff(bof_corner_nodes,borf_corner_nodes);
    bof_corner_nodes = setdiff(bof_corner_nodes,bolf_corner_nodes);
    % inner most bottom nodes
    c = [(max(G.nodes.coords(:,1)) + min(G.nodes.coords(:,1)))*0.5 ...
        (max(G.nodes.coords(:,2)) + min(G.nodes.coords(:,2)))*0.5];
    tol = 1e32;
    for i=1:numel(bottom_nodes)
        dist = norm(G.nodes.coords(bottom_nodes(i), 1:2)-c);
        if dist < tol
            tol = dist;
            bottom_innermost = bottom_nodes(i);
        end
    end
    bottom_nodes = setdiff(bottom_nodes,bottom_innermost);
end

