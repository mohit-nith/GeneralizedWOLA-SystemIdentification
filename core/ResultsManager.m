% ResultsManager.m
% -------------------------------------------------------------------------
% This class handles visualization and saving of results from Generalized WOLA-based Acoustic Echo Cancellation (AEC) simulations.
%
% Function: plotResults
%   - Plot ERLE (Echo Return Loss Enhancement) curves for:
%       * 'SWDFT' – Single line plot of ERLE vs filter length
%       * 'PT-WOLA' – Matrix plots for varying difference-terms and cross-terms configurations
%
% Function: saveResults
%   - Saving ERLE plots to the output directory defined in AECParameters object
%   - Saving simulation data as structure:
%       * Results (ERLE, residuals)
%       * Room impulse response
%       * Analysis and synthesis windows
%       * Full PT-WOLA configuration
%
% Usage:
%   ResultsManager.plotResults(results, config, params);                        % Plot results and save them
%   ResultsManager.saveResults(results, rir, anaWin, synWin, config, params);   % Save simulation data
%
% Inputs:
%   - results   - Struct from AECSimulation containing performance metrics like ERLE, etc
%   - config:   - PTWOLAConfig object used during simulation
%   - params:   - AECParameters object (for output directory path)
%   - rir:      - Room impulse response vector
%   - anaWin:   - Analysis window
%   - synWin:   - Synthesis window
%
% Outputs:
%   - Saved plot image (if enabled)
%   - Saved .mat file containing complete simulation state
%
% Dependencies:
%   - None
%
% Author: Mohit Sharma
% Date: 01/2025
% -------------------------------------------------------------------------
classdef ResultsManager

    methods (Static)
        %% Display and optionally save the ERLE plot
        function plotResults(resultsStruct, config, params, SaveResults)
            % Default save results
            if nargin < 4
                SaveResults = true;
            end
            % Plot ERLE matrix from resultsStruct
            ERLE = resultsStruct.ERLE;

            % Results for SWDFT (vector)
            figure;
            if strcmpi(config.implementation, 'SWDFT')
                plot(config.Total_filter, ERLE, '-o', 'LineWidth', 1.5, 'MarkerIndices', 1:3:length(config.Total_filter));
                xlabel('Filter Length');
                ylabel('ERLE (dB)');
                title('ERLE vs Filter Length (SWDFT)');
                xlim([1, config.Total_filter(end) + 1]);
                xticks(config.Total_filter(1:2:end));
                ylim([min(ERLE, [], 'all')-2, max(ERLE, [], 'all')+2]);
                grid on;

                % Results for PTWOLA implementation (matrix)
            elseif strcmpi(config.implementation, 'PT-WOLA')
                hold on;
                Colors = lines(size(ERLE, 1));   % Auto color selection
                MarkerType = {'o','s','^','d','x','+','v','*','>','<'};  % Cycle markers

                switch lower(config.Plot_constant)

                    %  Case 1: Constant Cross-Terms, Varying nf 
                    case 'const_cross_terms'
                        for row = 1:size(ERLE, 1)
                            y_vals = ERLE(row, :);

                            % Skip rows that are entirely NaN
                            if all(isnan(y_vals))
                                continue;
                            end

                            % Move NaNs to beginning for alignment
                            nan_mask = isnan(y_vals);
                            y_vals = [y_vals(nan_mask), y_vals(~nan_mask)];

                            legend_label = sprintf('R = %d', config.cross_length_matrix(row));
                            title_adjustment = 'nf = T - 2R';
                            plot(config.Total_filter, y_vals, ...
                                'Color', Colors(row,:), ...
                                'LineWidth', 1.5, ...
                                'Marker', MarkerType{mod(row-1,length(MarkerType))+1}, ...
                                'DisplayName', legend_label, ...
                                'MarkerIndices', 1:3:length(config.Total_filter));
                        end
                        xlabel('Total filter length T');

                        %  Case 2: Constant Filter Taps, Varying R 
                    case 'const_diff_terms'
                        for row = 1:size(ERLE, 1)
                            y_vals = ERLE(row, :);

                            % Skip rows that are entirely NaN
                            if all(isnan(y_vals))
                                continue;
                            end

                            % Move NaNs to beginning for alignment
                            nan_mask = isnan(y_vals);
                            y_vals = [y_vals(nan_mask), y_vals(~nan_mask)];

                            legend_label = sprintf('nf = %d', config.diff_filter_length(row));
                            title_adjustment = 'R = (T - nf)/2';
                            plot(config.Total_filter, y_vals, ...
                                'Color', Colors(row,:), ...
                                'LineWidth', 1.5, ...
                                'Marker', MarkerType{mod(row-1,length(MarkerType))+1}, ...
                                'DisplayName', legend_label, ...
                                'MarkerIndices', 1:3:length(config.Total_filter));
                        end
                        xlabel('Total filter length T');

                    otherwise
                        warning('Unknown Plot_constant type: %s', config.Plot_constant);
                end

                %  Common plot details 
                ylabel('ERLE (dB)');
                title(sprintf('ERLE vs Total Filter Length | Config : %s | %s', strrep(config.Plot_constant, '_', ' '), title_adjustment));
                legend('NumColumns', 2, 'Location', 'northwest');
                xlim([1, config.Total_filter(end) + 1]);
                xticks(config.Total_filter(1:2:end));
                ylim([min(ERLE, [], 'all')-2, max(ERLE, [], 'all')+2]);
                grid on;
                hold off;

            else
                warning('Unsupported implementation for ERLE plot: %s', config.implementation);
                return;
            end

            % Save plot if SaveResults = true
            if SaveResults
                %  Ensure directory exists 
                if ~exist(params.OutputDir, 'dir')
                    mkdir(params.OutputDir);
                end

                fname = sprintf('ERLE_plot_%s_%s_%s.png', config.implementation, lower(config.analysis_window), lower(config.Synthesis_window_def(1:4)));
                saveas(gcf, fullfile(params.OutputDir, fname));
                fprintf('ERLE plot saved: %s\n', fname);
            end
        end

        %% Optionally save the simulation data
        function saveResults(resultsStruct, rir, anaWin, synWin, config, params)
            % Save all simulation data into a single structure
            SimulationData.results = resultsStruct;
            SimulationData.rir = rir;
            SimulationData.anaWin = anaWin;
            SimulationData.synWin = synWin;
            SimulationData.config = config;

            % Save results to file with name format: RIR length -> analysis window -> synthesis window generation method)
            fname = sprintf('RIR_L_%d_%s_%s',length(rir), lower(config.analysis_window), lower(config.Synthesis_window_def(1:4)));
            fullPath = fullfile(params.OutputDir, [fname '.mat']);

            %  Save file 
            save(fullPath, 'SimulationData', '-v7.3');
            fprintf('Results saved: %s\n', fullPath);
        end


    end
end