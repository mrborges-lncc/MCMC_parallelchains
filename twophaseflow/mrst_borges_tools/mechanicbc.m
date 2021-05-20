el_bc = mechanicbc(optcase, overburden, top_faces, bottom_nodes, ...
    top_nodes, right_nodes, left_nodes, front_nodes, back_nodes, ...
    fr_corner_nodes, fl_corner_nodes, br_corner_nodes, ...
    bl_corner_nodes, bottom_innermost)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %switch optcase
        %==================================================================
        %case 'load + roller'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % zero vertical displacement for bottom nodes and zero 
            % displacement for  innermost nodes to anchor the problem
            Nbo = numel(bottom_nodes);
            el_bc.disp_bc.nodes = bottom_nodes;
            el_bc.disp_bc.uu    = repmat([0, 0, 0], Nbo, 1);
            el_bc.disp_bc.mask  = repmat([false, false, true], Nbo, 1);
            el_bc.disp_bc.mask(bottom_innermost, 1:2) = true;
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            % force applied at top
            top_force = overburden;
            el_bc.force_bc.faces = top_faces;
            el_bc.force_bc.force = repmat([0, 0, top_force], numel(top_faces), 1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %==================================================================
        %otherwise
            %error('bc_cases not recognized')
    %end