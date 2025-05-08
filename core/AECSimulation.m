% AECSimulation.m
% -------------------------------------------------------------------------
% This class executes the Acoustic Echo Cancellation (AEC) simulation using a Generalized Weighted OverLap-Add (WOLA) framework.
% It supports two implementation approaches:
%   - 'swdft' – Full-complexity Generalized WOLA filtering using standard Sliding DFT (SWDFT)
%   - 'pt-wola' – Efficient Generalized WOLA using a Per-Tone structure with cross-term and difference-term as inputs
%
% Usage:
%   sim = AECSimulation(signals, ptwolaConfig, anaMic, synMic);    % Initialize
%   results = sim.run();                                           % Run simulation
%
% Inputs:
%   signals         – Structure containing y (mic), far_end, far_end_speaker
%   ptwolaConfig    – Instance of PTWOLAConfig defining algorithm and filter settings
%   anaWin          – Analysis window
%   synWin          – Synthesis window
%
% Outputs:
%   results         – Structure containing:
%       * ERLE              – Echo Return Loss Enhancement
%       * ResidualError     – Residual echo signal(s)
%       * FilterEstimates   – Estimated subband filter coefficients
%
% Dependencies:
%   - Generalized_WOLA_AEC_Processor.run.m for adaptive filtering
%   - padcat (for padding variable-length ERLE arrays into matrix form)
%   - Requires prior generation of windows and input signals
%
% Author: Mohit Sharma
% Date: 01/2025
% -------------------------------------------------------------------------
classdef AECSimulation
    properties
        Signals         % Structure containing far-end (echo), near-end, mic signals
        PTWOLAConfig    % PTWOLA Configuration
        AnaWin          % Analysis window
        SynWin          % Synthesis window
    end

    methods
        function obj = AECSimulation(signals, ptwolaCfg, anaWin, synWin)
            % Constructor to initialize the simulation environment
            obj.Signals = signals;
            obj.PTWOLAConfig = ptwolaCfg;
            obj.AnaWin = anaWin;
            obj.SynWin = synWin;
        end

        function results = run(obj)
            % Runs the AEC simulation loop and returns result metrics (structure with ERLE, ResidualError, FilterEstimates (optional))
            % For readability extract structures from object
            sig = obj.Signals;
            PT_WOLA = obj.PTWOLAConfig;
            anaWindow = obj.AnaWin;
            synWindow = obj.SynWin;

            switch lower(PT_WOLA.implementation)

                %%  Case 1: PTWOLA implementation
                case 'pt-wola'
                    results = obj.runPTWOLA(sig, PT_WOLA, anaWindow, synWindow);
                                fprintf('Simulation complete for Generalized WOLA with %s implementation | Config: %s \n', upper(PT_WOLA.implementation), strrep(PT_WOLA.Plot_constant, '_', ' '));

                    %%  Case 2: SWDFT (Standard generalized WOLA implementation with sliding DFT)
                case 'swdft'
                    results = obj.runSWDFT(sig, PT_WOLA, anaWindow, synWindow);
                                fprintf('Simulation complete for Generalized WOLA with %s implementation \n', upper(PT_WOLA.implementation));

                otherwise
                    error('Unknown implementation type: %s', PT_WOLA.implementation);
            end

        end

        function results = runPTWOLA(~, sig, PT_WOLA, anaWindow, synWindow)
            % Runs generalized WOLA with low-complexity PTWOLA implementation using cross and difference terms.
            % Supports both Constant-Difference-Terms and Constant-Cross-Terms modes
            if strcmpi(PT_WOLA.Plot_constant, 'Const_Diff_Terms')
                % Define valid number of cross terms in cross_length_matrix for each row (for each PTWOLA filter length)
                matrix = PT_WOLA.cross_length_matrix;
                validPerRow = sum(matrix >= 0, 2);
                secondaryMatrix = PT_WOLA.diff_filter_length;
            else
                % Define valid number of PTEQ filter taps in diff_filter_length matrix for each row (for each cross-term length)
                matrix = PT_WOLA.diff_filter_length;
                validPerRow = sum(matrix >= 1, 2);
                secondaryMatrix = PT_WOLA.cross_length_matrix;
            end


            % Preallocate memory to store ERLE vs valid crossterms (or difference terms) for each difference term filter (or each cross-term filter) length in each cell
            ERLE_cell = cell(size(matrix, 1), 1);
            Res_cell = cell(size(matrix, 1), 1);

            % Process each row of matrix as a configuration group
            numRows = size(matrix,1);
            for rowIdx = 1:numRows
                % Preallocate memory to a temporary row vector which stores ERLE values (and Residual errors) within nested for loop.
                ERLE_row = NaN(1, validPerRow(rowIdx));
                num_frames = floor(length(sig.far_end) / PT_WOLA.N) - 1;
                Residual_row = NaN(num_frames, validPerRow(rowIdx));
                % Extract valid (non-negative) entries in the current row
                validEntries = matrix(rowIdx, matrix(rowIdx,:) >= 0 | matrix(rowIdx,:) >= 1);

                % Loop over each valid configuration column for the current row
                for colIdx = 1:length(validEntries)
                    % Dynamically allocate difference term length (nf) and single sided cross-term length (R) based on Const_Diff_Terms or Const_Cross_Terms
                    nf = strcmpi(PT_WOLA.Plot_constant, 'Const_Diff_Terms') * secondaryMatrix(rowIdx) + ...
                        strcmpi(PT_WOLA.Plot_constant, 'Const_Cross_Terms') * validEntries(colIdx);

                    R  = strcmpi(PT_WOLA.Plot_constant, 'Const_Diff_Terms') * validEntries(colIdx) + ...
                        strcmpi(PT_WOLA.Plot_constant, 'Const_Cross_Terms') * secondaryMatrix(rowIdx);

                    % Run the adaptive filtering (or LSQR) function for current (nf, R) configuration
                    [erle, residual, ~] = Generalized_WOLA_AEC_Processor.run(sig.y, sig.far_end_speaker, sig.far_end, nf, R, ...
                        PT_WOLA.M, PT_WOLA.analysis_window, PT_WOLA.Synthesis_window_def, PT_WOLA.implementation, PT_WOLA.algo, anaWindow, synWindow);

                    % Store the computed ERLE and residual error for this config
                    ERLE_row(colIdx) = erle;
                    Residual_row(:, colIdx) = residual;
                end

                ERLE_cell{rowIdx} = ERLE_row;
                Residual_cell{rowIdx} = Residual_row;
            end

            % Pad cell array into matrix
            results.ERLE = padcat(ERLE_cell{:});
            results.ResidualError = Residual_cell;
        end

        function results = runSWDFT(~, sig, PT_WOLA, anaWindow, synWindow)
            % Runs full complexity Generalized WOLA with sliding DFT terms
            % List of all filter lengths
            nfList = PT_WOLA.Total_filter;
            ERLE = NaN(length(nfList), 1);
            num_frames = floor(length(sig.far_end) / PT_WOLA.N) - 1;
            Residual = NaN(num_frames, length(nfList));
            F_est = cell(length(nfList), 1);

            for i = 1:length(nfList)
                nf = nfList(i);

                % Run the adaptive filtering (or LSQR) function for current (nf, R) configuration
                [erle, resErr, F] = Generalized_WOLA_AEC_Processor.run(sig.y, sig.far_end_speaker, sig.far_end, nf, 0, ...
                    PT_WOLA.M, PT_WOLA.analysis_window, PT_WOLA.Synthesis_window_def, PT_WOLA.implementation, PT_WOLA.algo, anaWindow, synWindow);

                % Store the computed ERLE and residual error for this config
                ERLE(i) = erle;
                Residual(:, i) = resErr;
                F_est{i} = F;
            end

            % Collect output in structure
            results.ERLE = ERLE;
            results.ResidualError = Residual;
            results.FilterEstimates = F_est;
        end
    end
end