% run_WOLA_AEC_Simulation.m
% -------------------------------------------------------------------------
% This script implements the Generalized (or Per-Tone) WOLA- based Acoustic Echo Cancellation (AEC) simulation. It can:
%   - Run the entire simulation pipeline in a single step (for typical use)
%   - Step through each stage 
%
% Author: Mohit Sharma
% Date: 01/2025
% -------------------------------------------------------------------------

clc; clear;
addpath(genpath(pwd));

% Simulation Control Flags
RunFullSimulation = false;      % true → full pipeline, false → step-by-step (debug) mode
SaveResults       = true;       % true → save plots and .mat result files to Params.OutputDir. % false → skip saving

% Instantiate Main Controller Object
sys = Generalized_WOLA_AEC_System();
if RunFullSimulation
    %% Executes the full pipeline (RIR → Signals → Windows → Simulation → Save/Plot)
    sys = sys.runAll(SaveResults); 

else
    %% Run each stage independently

    % Generate Room Impulse Response (RIR)
    sys.RIR = Generalized_WOLA_AEC_System.generateRIR(sys.Params);

    % Generate Far-End, Near-End, and Mic Signals
    sys.Signals = Generalized_WOLA_AEC_System.generateSignals(sys.Params, sys.RIR);

    % Design Analysis and Synthesis Windows
    [sys.AnaWin, sys.SynWin] = Generalized_WOLA_AEC_System.generateWindows(sys.PTWOLAConfig);

    % Run the AEC Simulation and and Generate Results
    sys.Results = Generalized_WOLA_AEC_System.RunSimulation(sys.Signals, sys.PTWOLAConfig, sys.AnaWin, sys.SynWin);

    % Save and Plot Results
    sys.saveAndPlot(sys.Results, sys.RIR, sys.AnaWin, sys.SynWin,SaveResults);
end