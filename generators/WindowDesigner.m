% WindowDesigner.m
% -------------------------------------------------------------------------
% This class generates analysis and synthesis windows for WOLA-based Acoustic Echo Cancellation (AEC) systems, based on the
% configuration provided in a PTWOLAConfig object.
%
% The window generation includes:
%   - Analysis window options:
%       * 'Rectangular' – Flat window scaled by sqrt(0.5)
%       * 'Cosine' – Standard periodic Hann window
%       * 'Sqrt-Hann' – Square root of a periodic Hann window
%   - Synthesis window design options:
%       * 'Distortion_minimizing' – Uses optimization for perfect reconstruction with minimum signal distortion. Method is based on M.
%       Sharma and M. Moonen, "Prototype filter design for weighted overlap-add filter bank based sub-band adaptive filtering
%       applications," 2023 31st European Signal Processing Conference (EUSIPCO), Helsinki, Finland, 2023, pp. 366-370.
%       * 'Norm_minimizing' – Conventional method. Ensures Constant OverLap-Add (COLA) property and norm minimizing window design.
%
% Usage:
%   winGen = WindowDesigner(ptwolaConfig);        % Create window designer object
%   [anaWin, synWin] = winGen.generate();         % Generate windows
%
% Inputs:
%   ptwolaConfig – Instance of PTWOLAConfig defining window parameters
%
% Outputs:
%   anaWin – Analysis window (column vector of length M)
%   synWin – Synthesis window (column vector of length M)
%
% Dependencies:
%   - For distortion-minimizing synthesis: Requires `PTEQ_WOLA_FB_window_generation_fcn_opt_tool.m` and MatLab optimization toolbox
%
% Author: Mohit
% Date: 01/2025
% -------------------------------------------------------------------------
classdef WindowDesigner
    properties
        PTWOLAConfig   % Instance of PTEQConfig class
    end

    methods
        function obj = WindowDesigner(ptwolaConfig)
            % Constructor with PTEQConfig as input
            obj.PTWOLAConfig = ptwolaConfig;
        end

        function [anaWin, synWin] = generate(obj)
            % Generate the analysis (anaWin) and synthesis (synWin) windows For readability extract structures from object
            M = obj.PTWOLAConfig.M;
            N = obj.PTWOLAConfig.N;
            ana_type = lower(obj.PTWOLAConfig.analysis_window);
            syn_def  = lower(obj.PTWOLAConfig.Synthesis_window_def);

            %%  Generate Analysis Window
            switch ana_type
                case 'rectangular'
                    anaWin = sqrt(0.5) * ones(M, 1);

                case 'cosine'
                    anaWin = hann(M, 'periodic');           % Equivalent to 0.5*(1-cos(2*pi*(0:N-1)/(N))).';

                case 'sqrt-hann'
                    anaWin = sqrt(hann(M, 'periodic'));     % Equivalent to sqrt(0.5*(1-cos(2*pi*(0:N-1)/(N-1)))).';

                otherwise
                    error('Unknown analysis window type: "%s"', ana_type);
            end

            %%  Generate Synthesis Window
            switch syn_def
                case 'distortion_minimizing'
                    if strcmpi(ana_type, 'rectangular')
                        % Analytical solution for rectangular case
                        synWin = [zeros(M/2,1); sqrt(2)*ones(M/2,1)];
                    else
                        % Optimized distortion minimizing synthesis window for a given analysis window satisfying PR. Needs optimization toolbox
                        synWin = design_wola_synthesis_window(anaWin(:));
                    end

                case 'norm_minimizing'
                    % For COLA-compatiblity setups- Works when both the analysis and synthesis windows are same.
                    if ~strcmpi(ana_type, 'cosine')
                        [isCola, medianSum, ~] = iscola(anaWin, N, 'WOLA');
                        if ~isCola
                            error('Analysis window does not satisfy COLA for N = %d', N);
                        end
                        % Scaling to ensure output is not scaled version of input in PR. Eg: when using sqrt(hann), the constant is 1.08 not 1.
                        synWin = anaWin / medianSum;
                    else
                        % For cosine analysis window the synthesis window shape is not same as analysis window shape
                        synWin = design_wola_synthesis_window(anaWin(:),syn_def);
                    end

                otherwise
                    warning('Unknown synthesis window def. Using empty.');
                    synWin = [];
            end

            fprintf('Windows created: Analysis - %s | Synthesis design criterion - %s\n', ...
                obj.PTWOLAConfig.analysis_window, ...
                strrep(obj.PTWOLAConfig.Synthesis_window_def, '_', ' '));
        end
    end
end