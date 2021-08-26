function [model,states,initState0,wellSols,W,dt] = runMechanicSchedule(G,...
    rock,E,nu,overburden,patm,depth,rhoR,grav,vinj,vprod,BHP,well_r,Tf,Tinj, ...
    Tinjup,Tinjstop,nstep,compW,compO,compR,ndt,printa,pdirichlet,varargin)
%% Example: Poroelasticity simulation applied to the Norne case.
% 
% The simulation options are gathered in the opt structure. If opt=[] the
% simulation is run with the default options defined below
%
% **  Summary of the options ** 
%
% option 'norne_case' :
%
%     * 'full'       : 7392 cells
%     * 'mini Norne' :  605 cells
%
% option 'bc_case' :
%
%     * 'no displacement' : All nodes belonging to external faces have displacement
%                           equal to zero
%     * 'bottom fixed'    : The nodes that belong to the bottom have zero
%                           displacement, while a given pressure is imposed on
%                           the external faces that are not bottom faces.
%
% option 'method' :
%
%     * 'fully coupled'          : The mechanical and flow equations are solved fully couplde.
%     * 'fixed stress splitting' : The mechanical and flow equations are solved
%                                  sequentially using a fixed stress splitting
%
% option 'fluid_model' :
%
%     * 'blackoil'  : blackoil model is used for the fluid (gas is injected, see
%                     schedule below)
%     * 'oil water' : Two phase oil-water
%     * 'water'     : water model is used for the fluid

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui

    % setup default option values
    opt = struct('norne_case'         , 'full'          , ...
                 'bc_case'            , 'bottom fixed'  , ...
                 'method'             , 'fully coupled' , ...
                 'fluid_model'        , 'water'         , ...
                 'nonlinearTolerance' , 1e-6            , ...
                 'splittingTolerance' , 1e-3            , ...
                 'verbose'            , false           , ...
                 'splittingVerbose'   , false);
    opt = merge_options(opt, varargin{:});

    % overwrite the default options by the given option
    % optvals = cellfun(@(x) opt.(x), fieldnames(opt), 'uniformoutput', false);
    % optlist = reshape(vertcat(fieldnames(opt)', optvals'), [], 1);
    % opt = merge_options(default_opt, optlist{:});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Gravity
    % The gravity in this option affects only the fluid behavior
    if grav == 1
        gravity reset on;
    else
        gravity off;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params = poroParams(mean(rock.poro), true, 'E', mean(E),...
        'nu', mean(nu), 'alpha', mean(rock.alpha), 'K_f', 1/compO);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute overburden
    TOL = 1.0e-07;
    if overburden < TOL
        g    = norm(gravity);
        overburden = patm + depth * rhoR * g;
    end
    ptop = overburden * params.gamma;
    fprintf('\n==============================================================\n\n')
    fprintf('Load at the top of the reservoir.........: %4.1f MPa\n',overburden/mega);
    fprintf('Fluid pressure at the top of reservoir...: %4.1f MPa\n',ptop/mega)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Setup fluid parameters from SPE1
    switch opt.fluid_model
        case 'blackoil'
            pth = getDatasetPath('spe1');
            fn  = fullfile(pth, 'BENCH_SPE1.DATA');
            deck = readEclipseDeck(fn);
            deck = convertDeckUnits(deck);
            fluid = initDeckADIFluid(deck);
            if isfield(fluid, 'pcOW')
                fluid = rmfield(fluid, 'pcOW');
            end
            if isfield(fluid, 'pcOG')
                fluid = rmfield(fluid, 'pcOG');
            end
            % Setup quadratic relative permeabilities, since SPE1 relperm are a bit rough.
            fluid.krW  = @(s) s.^2;
            fluid.krG  = @(s) s.^2;
            fluid.krOW = @(s) s.^2;
            fluid.krOG = @(s) s.^2;
            pRef = deck.PROPS.PVTW(1);
            
        case {'oil water'}
            fluid = initSimpleADIFluid('phases', 'WO', ...
                'mu', [1, 20]*centi*poise, ...
                'n',  [2, 2], 'rho', [1000.0, 860.0]*kilogram/meter^3,...
                'c', [compW, compO], 'cR', mean(compR), ...
                'pRef', patm);
            funrho = @(p) fluid.rhoOS * fluid.bO(p);
            if printa == 1
                p = linspace(patm, overburden*1.2,100)';
                two2Dplot(p/barsa, [fluid.rhoWS * fluid.bW(p) fluid.rhoOS * fluid.bO(p)],...
                    '$p\ (bar)$','$\rho_{\alpha}$','water','oil',1);
                pause(3); clf; close all
            end        
            
        case {'water'}
            fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise,...
                'rho', 1000.0*kilogram/meter^3, 'c', compW, ...
                'cR', mean(compR), 'pRef', patm);
            funrho = @(p) fluid.rhoWS * fluid.bW(p);
        
        otherwise
            error('fluid_model  not recognized.');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Setup boundary conditions for mechanics (no displacement)
    [top_faces, bottom_nodes, top_nodes, right_nodes, left_nodes, ...
        front_nodes, back_nodes, fr_corner_nodes, fl_corner_nodes, ...
        br_corner_nodes, bl_corner_nodes, tr_corner_nodes, ...
        tl_corner_nodes, tf_corner_nodes, tb_corner_nodes, ...
        bor_corner_nodes, bol_corner_nodes, bof_corner_nodes, ...
        bob_corner_nodes, trf_corner_nodes, trb_corner_nodes, ...
        tlf_corner_nodes, tlb_corner_nodes, borf_corner_nodes, ...
        borb_corner_nodes, bolf_corner_nodes, bolb_corner_nodes, ...
        bottom_innermost] = find_boudary_nodes(G);
    
    el_bc = mechanic_bc(opt.bc_case, G, overburden, top_faces, ...
        bottom_nodes, top_nodes, right_nodes, left_nodes, front_nodes, ...
        back_nodes, fr_corner_nodes, fl_corner_nodes, br_corner_nodes, ...
        bl_corner_nodes, tr_corner_nodes, tl_corner_nodes, ...
        tf_corner_nodes, tb_corner_nodes, bor_corner_nodes, ...
        bol_corner_nodes, bof_corner_nodes, bob_corner_nodes, ...
        trf_corner_nodes, trb_corner_nodes, tlf_corner_nodes, ...
        tlb_corner_nodes, borf_corner_nodes, borb_corner_nodes, ...
        bolf_corner_nodes, bolb_corner_nodes, bottom_innermost);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Setup load for mechanics
    % In this example we do not impose any volumetric force
    %loadfun = @(x) (0*x);
    % in the default case, we ignore gravity, but we need to specify a load
    % function nevertheless.
    loadfun = @(x) repmat(rhoR * gravity(), size(x, 1), 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Gather all the mechanical parameters in a struct
    mech = struct('E', E, 'nu', nu, 'el_bc', el_bc, 'load', loadfun);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute BHPressure
    switch opt.fluid_model
        case {'oil water'}
            if BHP < TOL
                pRef  = patm;
                g     = norm(gravity);
                rhoO  = @(p) fluid.rhoOS * fluid.bO(p);
                equil = ode23(@(z,p) g.* rhoO(p), ...
                    [0, max(G.nodes.coords(:,3))], pRef);
                BHPressure = reshape(deval(equil, ...
                    max(G.nodes.coords(:,3))), [], 1);
                BHPressure = reshape(deval(equil, ...
                    max(G.cells.centroids(:,3))), [], 1);  clear equil
            else
                BHPressure = BHP;
            end
        case{'water'}
            if BHP < TOL
                pRef  = patm;
                g     = norm(gravity);
                rhoWR = fluid.bW(pRef)*fluid.rhoWS;
                rhoW  = @(p) fluid.rhoWS * fluid.bW(p);
                equil = ode23(@(z,p) g.* rhoW(p), ...
                    [0, max(G.nodes.coords(:,3))], pRef);
                BHPressure = reshape(deval(equil, ...
                    max(G.cells.centroids(:,3))), [], 1);  clear equil
            else
                BHPressure = BHP;
            end
    end
    fprintf('BHP (Bottom hole pressure)...............: %4.1f MPa\n',BHPressure/mega);
    fprintf('\n==============================================================\n')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Setup model

    modeltype = [opt.method, ' and ', opt.fluid_model];
    fullycoupledOptions = {'verbose', opt.verbose};
    splittingOptions = {'splittingTolerance', opt.splittingTolerance, ...
                        'splittingVerbose', opt.splittingVerbose};
    switch modeltype

      case 'fully coupled and blackoil'
        model = MechBlackOilModel(G, rock, fluid, mech, ...
            fullycoupledOptions{:});

      case 'fixed stress splitting and blackoil'
        model = MechFluidFixedStressSplitModel(G, rock, fluid, mech, ...
                                               'fluidModelType', 'blackoil', ...
                                               splittingOptions{:});

      case 'fully coupled and oil water'
        model = MechOilWaterModel(G, rock, fluid, mech,...
            fullycoupledOptions{:});

      case 'fixed stress splitting and oil water'
        model = MechFluidFixedStressSplitModel(G, rock, fluid, mech, ...
                                               'fluidModelType', 'oil water', ...
                                               splittingOptions{:});

      case 'fully coupled and water'
        model = MechWaterModel(G, rock, fluid, mech, fullycoupledOptions{: });

      case 'fixed stress splitting and water'
        model = MechFluidFixedStressSplitModel(G, rock, fluid, mech, ...
                                               'fluidModelType', 'water', ...
                                               splittingOptions{:});

      otherwise
        error('modeltype not recognized.');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Setup wells
    L  = max(G.nodes.coords);
    Lx = L(1); Ly = L(2); Lz = L(3);
    nx = G.cartDims(1); ny = G.cartDims(2); nz = G.cartDims(3);
    [inj_cells,prod1_cells,prod2_cells,prod3_cells,prod4_cells] = ...
        fivespot_wells(G,Lx,Ly,nx,ny);

    refdepth = min(G.cells.centroids(:, 3)); % for example...
%     refdepth = min(G.cells.centroids(:, 3)); % for example...
    W = []; bc = [];
%    
    W = addWell(W, G, rock, inj_cells,...
        'Type', 'rate', 'Val', vinj,...
        'Comp_i', [1, 0, 0],...
        'Radius', well_r, 'Dir', 'z',...
        'Sign',1,'Name', '$w_{inj}$',...
        'refDepth', refdepth); %m^3/sec
    if abs(vprod) < TOL
        if pdirichlet
            [bc, bcfaces, bcpressures] = bc_fivespotDirichlet(G, BHPressure, ...
                funrho, norm(gravity), init_sat);
    %         prodcells = [prod1_cells,prod2_cells,prod3_cells,prod4_cells];
    %         equil = ode23(@(z,p) g.* funrho(p), ...
    %             [max(G.cells.centroids(:,3)), 0], BHP);
    %         initState.pressure(prodcells) = reshape(deval(equil, ...
    %             G.cells.centroids(prodcells,3)), [], 1);  clear equil        
        else
            W = addWell(W, G, rock, prod1_cells,...
                'Type', 'bhp', 'Comp_i', [0, 1, 0],...
                'Val', BHPressure, 'Radius', well_r*meter,...
                'Dir', 'z', 'name', '$w1$',...
                'refDepth', refdepth); %Pascal
            W = addWell(W, G, rock, prod2_cells,...
                'Type', 'bhp', 'Comp_i', [0, 1, 0],...
                'Val', BHPressure, 'Radius', well_r*meter,...
                'Dir', 'z', 'name', '$w2$',...
                'refDepth', refdepth); %Pascal
            W = addWell(W, G, rock, prod3_cells,...
                'Type', 'bhp', 'Comp_i', [0, 1, 0],...
                'Val', BHPressure, 'Radius', well_r*meter,...
                'Dir', 'z', 'name', '$w3$',...
                'refDepth', refdepth); %Pascal
            W = addWell(W, G, rock, prod4_cells,...
                'Type', 'bhp', 'Comp_i', [0, 1, 0],...
                'Val', BHPressure, 'Radius', well_r*meter,...
                'Dir', 'z', 'name', '$w4$',...
                'refDepth', refdepth); %Pascal
        end
    else
        W = addWell(W, G, rock, prod1_cells,...
            'Type', 'rate', 'Comp_i', [0, 1, 0],...
            'Val', -vprod, 'Radius', well_r*meter,...
            'Dir', 'z', 'name', '$w1$',...
            'refDepth', refdepth); %Pascal
        W = addWell(W, G, rock, prod2_cells,...
            'Type', 'rate', 'Comp_i', [0, 1, 0],...
            'Val', -vprod, 'Radius', well_r*meter,...
            'Dir', 'z', 'name', '$w2$',...
            'refDepth', refdepth); %Pascal
        W = addWell(W, G, rock, prod3_cells,...
            'Type', 'rate', 'Comp_i', [0, 1, 0],...
            'Val', -vprod, 'Radius', well_r*meter,...
            'Dir', 'z', 'name', '$w3$',...
            'refDepth', refdepth); %Pascal
        W = addWell(W, G, rock, prod4_cells,...
            'Type', 'rate', 'Comp_i', [0, 1, 0],...
            'Val', -vprod, 'Radius', well_r*meter,...
            'Dir', 'z', 'name', '$w4$',...
            'refDepth', refdepth); %Pascal
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch opt.fluid_model
        case 'blackoil'
            error('falta ajustar');
            W(1).compi = [0, 0, 1];
            W(2).compi = [0, 1, 0];
        case 'oil water'
            for i = 1:numel(W)
                tipo = W(i).type;
                if tipo(1:3) == 'bhp'
                    W(i).compi = [0 1];
                else
                    W(i).compi = [1 0];
                end
            end
        case 'water'
            for i = 1:numel(W)
                W(i).compi = [1];
            end
        otherwise
            error('fluid_model not recognized.')
    end

    facilityModel = FacilityModel(model.fluidModel);
    facilityModel = facilityModel.setupWells(W);
    model.FacilityModel = facilityModel;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Setup schedule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Tf       = Tf*day;
    Tinj     = Tinj*day;
    Tinjup   = Tinjup*day;
    Tinjstop = Tinjstop*day;
    [schedule1,schedule2,schedule3,schedule4,nstep1,nstep2,...
    nstep3,nstep4,dt] = injectionScheme(bc,W,vinj,...
        Tf,Tinj,Tinjup,Tinjstop,nstep,ndt);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Setup initial state
    clear initState;
    
    switch opt.fluid_model
        case 'blackoil'
            initState.pressure = pRef*ones(G.cells.num, 1);
            init_sat = [0, 1, 0];
            initState.rs  = 0.5*fluid.rsSat(initState.pressure);
        case 'oil water'
            init_sat = [0, 1];
            rhoO  = @(p) fluid.rhoOS * fluid.bO(p);
            equil = ode23(@(z,p) g.* rhoO(p), ...
                [min(G.nodes.coords(:,3)), max(G.nodes.coords(:,3))], ptop);
            initState.pressure = reshape(deval(equil, ...
                    G.cells.centroids(:,3)), [], 1);  clear equil

        case 'water'
            init_sat = [1];
            rhoW  = @(p) fluid.rhoWS * fluid.bW(p);
            equil = ode23(@(z,p) g.* rhoW(p), ...
                [min(G.nodes.coords(:,3)), max(G.nodes.coords(:,3))], ptop);
            initState.pressure = reshape(deval(equil, ...
                    G.cells.centroids(:,3)), [], 1);  clear equil
        otherwise
            error('fluid_model not recognized.')
    end
    initState.xd = zeros(nnz(~model.mechModel.operators.isdirdofs), 1);
    initState.s  = ones(G.cells.num, 1)*init_sat;
    initState    = computeInitDisp(model, initState, zeros(G.nodes.num,3));
    initState    = addDerivedQuantities(model.mechModel, initState);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    initState0   = initState;
    states       = []; wellSols = [];
    solver       = NonLinearSolver('maxIterations', 100);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tempo = 0.0;
    %% First period
    aux = [];
    if nstep1 > 0
        [a p1] = max(G.nodes.coords(prod1_cells,3));
        [a p2] = min(G.nodes.coords(prod2_cells,3));
        [a p3] = min(G.nodes.coords(prod3_cells,3));
        [a p4] = min(G.nodes.coords(prod4_cells,3));
        bhp_pos= [prod1_cells(p1); prod2_cells(p2); prod3_cells(p3);...
            prod4_cells(p4)];
        BHP_0 = initState.pressure(bhp_pos);
        BHP   = [BHPressure; BHPressure; BHPressure; BHPressure];
        tt    = sum(schedule1.step.val);
        dt1   = schedule1.step.val;
        W(1).val = 0.0;
        for i=1:nstep1
            t = sum(dt1(1:i));
            v = t * (BHP - BHP_0) / tt + BHP_0;
            for j=2:5
                W(j).val = v(j-1);
            end
            schedule1.step.val     = dt1(i);
            schedule1.step.control = ones(numel(schedule1.step.val), 1);
            schedule1.control      = struct('W', W, 'bc', bc);
            [wellSols1, states1, schedulereport] = simulateScheduleAD(...
                initState, model, schedule1, 'nonlinearsolver', solver);
            states   = [states  ; states1];
            wellSols = [wellSols; wellSols1];
            initState= states{end};
            clear wellSols1 states1
            tempo = tempo + dt1(i)/day;
            fprintf('\nTIME...........................: %5.2f days\n',tempo)
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Second period
    if nstep2 > 0
        tt = sum(schedule2.step.val);
        dt2= schedule2.step.val;
        for i=1:nstep2
            t = sum(dt2(1:i));
            W(1).val = vinj*t/tt;
            schedule2.step.val     = dt2(i);
            schedule2.step.control = ones(numel(schedule2.step.val), 1);
            schedule2.control      = struct('W', W, 'bc', bc);
            [wellSols2, states2, schedulereport] = simulateScheduleAD(...
                initState, model, schedule2, 'nonlinearsolver', solver);
            states   = [states  ; states2];
            wellSols = [wellSols; wellSols2];
            initState= states{end};
            tempo = tempo + dt2(i)/day;
            fprintf('\nTIME...........................: %5.2f days\n',tempo)
        end
        clear wellSols2 states2
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Tird period
    if nstep3 > 0
        [wellSols3, states3, schedulereport] = simulateScheduleAD(...,
            initState, model, schedule3, 'nonlinearsolver', solver);
        states   = [states  ; states3];
        wellSols = [wellSols; wellSols3];
        initState= states{end};
        clear wellSols3 states3
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Fourth period
    if nstep4 > 0
        [wellSols4, states4, schedulereport] = simulateScheduleAD(...,
            initState, model, schedule4, 'nonlinearsolver', solver);
        states   = [states  ; states4];
        wellSols = [wellSols; wellSols4];
        clear wellSols4 states4
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
