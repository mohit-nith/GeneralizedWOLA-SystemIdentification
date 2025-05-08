% Generalized_WOLA_AEC_System.m
% -------------------------------------------------------------------------
% This is the main controller class for AEC simulations using the Generalized (or Per-Tone) WOLA framework.
% It integrates following stages:
%   1. AECParameters, PTWOLAConfig  - Parameter and configuration setup
%   2. RIRGenerator                 - Room Impulse Response (RIR) generation
%   3. SignalGenerator              - Far-end, near-end, and mic signal generation
%   4. WindowDesigner               - Analysis and synthesis window design
%   5. AECSimulation                - AEC simulation execution
%   6. ResultsManager               - Results visualization and saving
%
% Usage:
%   % Run full simulation with default settings
%   system = Generalized_WOLA_AEC_System();            % Initialize with default parameters
%   system = system.runAll();                          % Run full simulation and save results
%
%   % Run with custom parameters and save flag
%   params = AECParameters();
%   config = PTWOLAConfig();
%   system = Generalized_WOLA_AEC_System(params, config);
%   system = system.runAll(true);                     % true/false to enable/disable saving
%
% Inputs:
%   - params (optional)             – Custom AECParameters object
%   - config (optional)             – Custom PTWOLAConfig object
%   - SaveResults (optional)        – Boolean flag to save results (default: true)
%
% Outputs:
%   - Object properties updated with RIR, generated signals, windows, and results
%
% Dependencies:
%   - RIRGenerator, SignalGenerator, WindowDesigner, AECSimulation, ResultsManager
%
% Author: Mohit Sharma
% Date: 01/2025
% -------------------------------------------------------------------------
classdef Generalized_WOLA_AEC_System

    properties
        Params          % AECParameters object containing simulation settings
        PTWOLAConfig    % PTWOLAConfig object defining WOLA filter bank and adaptive algorithm parameters
        RIR             % Room Impulse Response (generated based on config)
        Signals         % Structure with far-end, near-end, and mic signals
        AnaWin          % Analysis window used in WOLA filter bank
        SynWin          % Synthesis window used in WOLA filter bank
        Results         % Structure with performance metrics like ERLE, ResidualError, etc.
    end

    % Static standalone utility methods
    methods (Static)
        function rir = generateRIR(Params)
            % Generate RIR using defined method in AECParameters
            rirGen = RIRGenerator(Params);
            rir = rirGen.generate();
        end

        function signals = generateSignals(Params,rir)
            % Generate far-end, near-end, echo, and mic signals
            sigGen = SignalGenerator(Params, rir);
            signals = sigGen.generate();
        end

        function [anaWin, synWin] = generateWindows(PTWOLAConfig)
            % Design analysis and synthesis windows based on settings in PTWOLAConfig
            winDesigner = WindowDesigner(PTWOLAConfig);
            [anaWin, synWin] = winDesigner.generate();
        end

        function results = RunSimulation(signals, PTWOLAConfig, anaWin, synWin)
            % Run the generalized WOLA-based AEC and return simulation results
            RunSim = AECSimulation(signals, PTWOLAConfig, anaWin, synWin);
            results = RunSim.run();
        end
    end

    % Dynamic methods
    methods
        function obj = Generalized_WOLA_AEC_System(params, config)
            % Constructor – Initialize with default or custom parameter/config objects
            if nargin < 1
                obj.Params = AECParameters();       % Use default parameters if not provided
            else
                obj.Params = params;
            end

            if nargin < 2
                obj.PTWOLAConfig = PTWOLAConfig();  % Use default configuration if not provided
            else
                obj.PTWOLAConfig = config;
            end
        end

        function saveAndPlot(obj, results, rir, anaWin, synWin, SaveResults)
            % Generate performance metrics plot and save results (if enabled)
            if ~exist('SaveResults', 'var')
                SaveResults = true;         % Save results if flag is not provided
            end
            ResultsManager.plotResults(results, obj.PTWOLAConfig, obj.Params, SaveResults);
            if SaveResults
                ResultsManager.saveResults(results, rir, anaWin, synWin, obj.PTWOLAConfig, obj.Params);
            end
        end

        function obj = runAll(obj, SaveResults)
            % Run full simulation pipeline: RIR → signals → windows → AEC → results
            if nargin < 2
                SaveResults = true; % Default saves results
            end

            % Generate Room Impulse Response (RIR)
            obj.RIR = Generalized_WOLA_AEC_System.generateRIR(obj.Params);

            % Generate Far-End, Near-End, and Mic Signals
            obj.Signals = Generalized_WOLA_AEC_System.generateSignals(obj.Params, obj.RIR);

            % Design Analysis and Synthesis Windows
            [obj.AnaWin, obj.SynWin] = Generalized_WOLA_AEC_System.generateWindows(obj.PTWOLAConfig);

            % Step 4: Run the AEC Simulation and Generate Results
            obj.Results = Generalized_WOLA_AEC_System.RunSimulation(obj.Signals, obj.PTWOLAConfig, obj.AnaWin, obj.SynWin);

            % Step 5: Save and Plot Results
            obj.saveAndPlot(obj.Results, obj.RIR, obj.AnaWin, obj.SynWin, SaveResults);

            disp('Full simulation complete.');
        end
    end


end