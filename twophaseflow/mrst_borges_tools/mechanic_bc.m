function el_bc = mechanic_bc(optcase, G, overburden, top_faces, ...
    bottom_nodes, top_nodes, right_nodes, left_nodes, front_nodes, ...
    back_nodes, fr_corner_nodes, fl_corner_nodes, br_corner_nodes, ...
    bl_corner_nodes, tr_corner_nodes, tl_corner_nodes, ...
    tf_corner_nodes, tb_corner_nodes, bor_corner_nodes, ...
    bol_corner_nodes, bof_corner_nodes, bob_corner_nodes, ...
    trf_corner_nodes, trb_corner_nodes, tlf_corner_nodes, ...
    tlb_corner_nodes, borf_corner_nodes, borb_corner_nodes, ...
    bolf_corner_nodes, bolb_corner_nodes, bottom_innermost)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch optcase
        %==================================================================
        case 'load + roller'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % zero vertical displacement for bottom nodes and zero 
            % displacement for  innermost nodes to anchor the problem
            Nbo = numel(bottom_nodes);
            el_bc.disp_bc.nodes = bottom_nodes;
            el_bc.disp_bc.uu    = repmat([0, 0, 0], Nbo, 1);
            el_bc.disp_bc.mask  = repmat([true, true, true], Nbo, 1);
            % bottom_innermost
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; bottom_innermost];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], 1, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, true, true], 1, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % roller boundary conditions
            % right
            Nr = numel(right_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; right_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Nr, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, false, false], Nr, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % left
            Nl = numel(left_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; left_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Nl, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, false, false], Nl, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % front
            Nf = numel(front_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; front_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Nf, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([false, true, false], Nf, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % back
            Nb = numel(back_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; back_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Nb, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([false, true, false], Nb, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % fr
            Nfr = numel(fr_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; fr_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Nfr, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, true, false], Nfr, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % fl
            Nfl = numel(fl_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; fl_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Nfl, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, true, false], Nfl, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % br
            Nbr = numel(br_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; br_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Nbr, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, true, false], Nbr, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % bl
            Nbl = numel(bl_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; bl_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Nbl, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, true, false], Nbl, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % tr
            Ntr = numel(tr_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; tr_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Ntr, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, false, false], Ntr, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % tl
            Ntl = numel(tl_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; tl_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Ntl, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, false, false], Ntl, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % tb
            Ntb = numel(tb_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; tb_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Ntb, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([false, true, false], Ntb, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % tf
            Ntf = numel(tf_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; tf_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Ntf, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([false, true, false], Ntf, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % bor
            Nbor = numel(bor_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; bor_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Nbor, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, false, true], Nbor, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % bol
            Nbol = numel(bol_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; bol_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Nbol, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, false, true], Nbol, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % bof
            Nbof = numel(bof_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; bof_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Nbof, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([false, true, true], Nbof, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % bob
            Nbob = numel(bob_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; bob_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Nbob, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([false, true, true], Nbob, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TRF
            Ntrf = numel(trf_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; trf_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Ntrf, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, true, false], Ntrf, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TLF
            Ntlf = numel(tlf_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; tlf_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Ntlf, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, true, false], Ntlf, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TRB
            Ntrb = numel(trb_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; trb_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Ntrb, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, true, false], Ntrb, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TLB
            Ntlb = numel(trf_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; tlb_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Ntlb, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, true, false], Ntlb, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % BoRB
            Nborb = numel(borb_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; borb_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Nborb, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, true, true], Nborb, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % BoLB
            Nbolb = numel(bolb_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; bolb_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Nbolb, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, true, true], Nbolb, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % BoRF
            Nborf = numel(borf_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; borf_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Nborf, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, true, true], Nborf, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % BoLF
            Nbolf = numel(bolf_corner_nodes);
            el_bc.disp_bc.nodes = [el_bc.disp_bc.nodes; bolf_corner_nodes];
            el_bc.disp_bc.uu    = [el_bc.disp_bc.uu;   repmat([0, 0, 0], Nbolf, 1)];
            el_bc.disp_bc.mask  = [el_bc.disp_bc.mask; repmat([true, true, true], Nbolf, 1)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % force applied at top
            top_force = overburden;
            signcoef = (G.faces.neighbors(top_faces, 1) == 0) - ...
                (G.faces.neighbors(top_faces, 2) == 0);
            n = bsxfun(@times, G.faces.normals(top_faces, :), signcoef./ ...
                G.faces.areas(top_faces));
            force = bsxfun(@times, n, top_force);
            el_bc.force_bc.faces = top_faces;
            el_bc.force_bc.force = force;
            %el_bc.force_bc.force = top_force * G.faces.normals(top_faces,:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %==================================================================    
        case 'bottom fixed'
            nx = G.cartDims(1);
            ny = G.cartDims(2);
            nz = G.cartDims(3);

            % Find the bottom nodes. On these nodes, we impose zero displacement

            c = zeros(nx*ny*nz, 1);
            c(G.cells.indexMap) = (1 : numel(G.cells.indexMap))';
            bottomcells = c(nx*ny*(nz - 1) +  (1 : (nx*ny))');
            faces = G.cells.faces(mcolon(G.cells.facePos(bottomcells), ...
                G.cells.facePos(bottomcells + 1) - 1), :);
            bottomfaces = faces( faces(:, 2) == 6  , 1);
            indbottom_nodes = mcolon(G.faces.nodePos(bottomfaces), ...
                G.faces.nodePos(bottomfaces + 1) - 1);
            bottom_nodes = G.faces.nodes(indbottom_nodes);
            isbottom_node = false(G.nodes.num, 1);
            isbottom_node(bottom_nodes) = true;
            bcnodes = find(isbottom_node);

            nn = numel(bcnodes);
            u = zeros(nn, 3);
            m = ones(nn, 3);
            disp_bc = struct('nodes', bcnodes, 'uu', u, 'mask', m);

            % Find outer faces that are not at the bottom. On these faces, we impose
            % a given pressure.

            is_outerface1 = (G.faces.neighbors(:, 1) == 0);
            is_outerface1(bottomfaces) = false;
            is_outerface2 = G.faces.neighbors(:, 2) == 0;
            is_outerface2(bottomfaces) = false;

            is_outerface = is_outerface1 | is_outerface2;

            outer_faces = find(is_outerface);

            outer_pressure = overburden;
            signcoef = (G.faces.neighbors(outer_faces, 1) == 0) - ...
                (G.faces.neighbors(outer_faces, 2) == 0);
            n = bsxfun(@times, G.faces.normals(outer_faces, :), signcoef./ ...
                       G.faces.areas(outer_faces));
            force = bsxfun(@times, n, outer_pressure);

            force_bc = struct('faces', outer_faces, 'force', force);
            
            el_bc = struct('disp_bc' , disp_bc, ...
                'force_bc', force_bc);

        %==================================================================
        case 'no displacement'
            ind = (G.faces.neighbors(:, 1) == 0 | G.faces.neighbors(:, 2) == 0);
            ind = find(ind);
            nodesind = mcolon(G.faces.nodePos(ind), G.faces.nodePos(ind + 1) - 1);
            nodes = G.faces.nodes(nodesind);
            bcnodes = zeros(G.nodes.num);
            bcnodes(nodes) = 1;
            bcnodes = find(bcnodes == 1);
            nn = numel(bcnodes);
            u = zeros(nn, 3);
            m = ones(nn, 3);
            disp_bc = struct('nodes', bcnodes, 'uu', u, 'mask', m);
            force_bc = [];
            
            el_bc = struct('disp_bc' , disp_bc, ...
                'force_bc', force_bc);
        %==================================================================    
        %==================================================================
        otherwise
            error('bc_cases not recognized')
    end