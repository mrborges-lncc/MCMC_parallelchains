function [description, options, state0, model, schedule, plotOptions, connectionIndices] ...
    = ensemble_base_problem_simple_gpsnet_model(varargin)
% Creates a GPSNET model with two-phase flow between two injectors and two
% producers.
%
% SYNOPSIS:
%   [description, options, state0, model, schedule, plotOptions, connectionIndices] ...
%        = ensemble_base_problem_simple_gpsnet_model('pn1', pv1, ...)
%
% DESCRIPTION:
%   
%   Very simple example that generates a 3D reservoir with a two-phase flow
%   problem. It has two injectors and two producers, of which one of the
%   injectors are not in a corner. May be used as a stand-alone example 
%   definition, or to construct an instance of `MRSTExample` as 
%   example = MRSTExample('ensemble_base_problem_3D_reservoir');
%   At the time of writing, the main purpose is to use this example for a
%   base in an example ensemble simulation.
%
% OPTIONAL PARAMETERS:
%   This example currently does not take any optional inputs
%
% RETURNS:
%   description - One-line example description, displayed in list-examples,
%                 and the only input argument if the function is called as
%                 description = my_example_wog()
%
%   options     - A struct of the optional input parameters, with defaults
%                 for all arguments that were not passed as optional
%                 parameters. Returned for convenient access to the example
%                 configuration.
%
%   state0, model, schedule - Initial state, model, and simulation schedule
%                             that can be passed to `simulateScheduleAD`
%
%   plotOptions - Cell array on the form {'pn1', pv1, ...} with arguments
%                 that can be used in any of the following ways
%                   - set(myAxis, 'pn1, vn1, ...)
%                   - figure('pn1', vn1, ...)
%                   - plotToolbar(G, state, 'pn1', vn1, ...)
%                 In addition to the standard optional parameters of
%                 `figure`, {'Size', [width, height]} can also be provided,
%                 which `MRSTExample` interprets as
%                 [pos(1:2), [width, height]], where
%                 pos = get(0, 'DefaultFigurePosition')
%
% SEE ALSO:
%   `MRSTExample`, `listExamples`, `exampleSuiteTutorial`

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

    % Each example must start with the description and options, followed by
    % an nargout check that returns if we only asked for the description
    % and options
    description = 'two-phase flow in fixed GPSNET reservoir connecting two injectors and two producers';
    
    % Each example can have any number of optional input arguments, and
    % must return a (possibly empy) options struct
    options = struct('deleteOldResults', false, ...
                     'cellsPerConnection', 10, ...
                     'gpsnetPoro', 0.2, ...
                     'gpsnetPerm', 1000*milli*darcy, ...
                     'plotNetwork', false, ...
                     'fullExampleName', 'ensemble_base_problem_3d_reservoir');
    options = merge_options(options, varargin{:});
    
    if nargout <= 2, return; end
    
    
    % Define any module dependencies for the example. The following are
    % typically always needed
    require ad-core ad-props ad-blackoil incomp dd-models ensemble diagnostics
    
    %% Define model
    % Create full 3D reservoir based on a full 3D example reservoir 
    fullExample = MRSTExample(options.fullExampleName);
    
    % Simulate the example so that we can do flow diagnostics
    fullProblem = fullExample.getPackedSimulationProblem();
    if options.deleteOldResults
        clearPackedSimulatorOutput(fullProblem, 'prompt', false);
    end
    
    [ok, status] = simulatePackedProblem(fullProblem);
    [fullWellSols, fullStates, reports] = getPackedSimulatorOutput(fullProblem);
    
    introAnimations = false;
    introPlots = false;
    
    if  introAnimations
        % Animate water saturation
        figure
        title('Water saturation');
        plotGrid(fullProblem.SimulatorSetup.model.G, 'FaceAlpha', 0, 'EdgeAlpha', 0.1);
        plotWell(fullProblem.SimulatorSetup.model.G, ...
                 fullProblem.SimulatorSetup.schedule.control.W);
        view(30, 50);
        pause(1);
        hs = []; % handle for saturation plot, empty initially
        for i = 1:size(fullProblem.SimulatorSetup.schedule.step.val, 1)
            hs = plotCellData(fullProblem.SimulatorSetup.model.G, ...
                              fullStates{i}.s(:,1), fullStates{i}.s(:,1) > 0.1);
            drawnow, pause(0.5);
        end
    end
    
    if introPlots
        plotWellSols(fullWellSols)
    end

    
    
    %% Create data-driven graph of well connections
    gpsnet = WellPairNetwork(fullProblem.SimulatorSetup.model, ...
                         fullProblem.SimulatorSetup.schedule, ...
                         fullStates, ...
                         fullProblem.SimulatorSetup.state0, ...
                         fullWellSols);
    gpsnet = gpsnet.filter_wps(1*stb/day);
    
    if options.plotNetwork
        figure;
        gpsnet.plotWellPairConnections();
    end
    
    %% Create the grid and rock representing the gpsnet model
    
    % Compute the required reservoir volume needed to fit all the oil in
    % the original model
    totalReservoirVolume = sum(fullProblem.SimulatorSetup.model.operators.pv)/options.gpsnetPoro;
    
    % Find number of well-pair connections
    numConnections = size(gpsnet.Graph.Edges, 1);
    
    % Create 2D grid of 1D connections
    length = (totalReservoirVolume*25)^(1/3);
    G = cartGrid([options.cellsPerConnection, 1, numConnections], ...
                 [length, length/5, length/5]*meter^3);
    G = computeGeometry(G);
    
    rock = makeRock(G, options.gpsnetPerm, options.gpsnetPoro);
    
    %% Create fluid and physical model
    fluid = fullProblem.SimulatorSetup.model.fluid;
    
    gravity off
    model = GenericBlackOilModel(G, rock, fluid);
    model.gas = false;
    model.OutputStateFunctions = {};
    
    [model, W, connectionIndices] = createDDmodel_1(model, ...
                                                    options.cellsPerConnection, ...
                                                    gpsnet.Graph, ...
                                                    fullProblem.SimulatorSetup.schedule.control.W);

    %% Define the suitable initial state and schedule
    
    state0 = initState(G, W, ...
                       fullProblem.SimulatorSetup.state0.pressure(1), ...
                       fullProblem.SimulatorSetup.state0.s(1,:));
                   
    schedule = simpleSchedule(fullProblem.SimulatorSetup.schedule.step.val, ...
                              'W', W);
    

    %% Plot options
    % plotOptions are only by MRSTExample. In case of empty plotOptions,
    % MRSTExample will attempt to set reasonable defaults
    plotOptions = {};
    
    %% Output connectionIndices
    %  To be able to run optimization algorithms or history matching on
    %  this model, we need to output the connectionIndices (lists of grid
    %  indices to identify each connection in G).
    %  Since this function is called from inside `MRSTExample`, we can't
    %  add new output variables, but has to attach it to some of those we
    %  have.
    %  Here, we will add the connectionIndices to the options output
    %  parameter.
    
    options.connectionIndices = connectionIndices;
    
end
