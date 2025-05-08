% PTWOLAConfig.m
% -------------------------------------------------------------------------
% This class defines the configuration structure for Generalized Weighted OverLap-Add (WOLA) processing using either a standard Sliding
% DFT (SWDFT)-based implementation or a low complexity Per-Tone WOLA (PT-WOLA) approach with cross-term or difference-term as inputs.
%
% The configuration defines parameters for the generalized WOLA based subband adaptive filtering process, including:
%   - Implementation method: 'SWDFT' or 'PT-WOLA'
%   - Total filter lengths to sweep over for evaluating performance
%   - Plotting mode: For PT-WOLA, controls how filter taps are allocated
%     between difference terms and single-sided cross terms. Useful for analyzing the individual roles of each term:
%       * 'Const_Cross_Terms' – Keeps cross-term filter taps constant, resulting in: Difference terms = Total taps - 2 × Cross terms
%       * 'Const_Diff_Terms' – Keeps difference-term filter taps constant, resulting in: Cross terms = (Total taps - Difference terms)/2
%   - Hop size and DFT sizes for WOLA (N: hop size, M: DFT size)
%   - Algorithm type for subband system identification:
%       'LSQR', 'RLS', or 'NLMS'
%   - Analysis and synthesis windows.
%       * Analysis window shape:
%           'Rectangular', 'Cosine', or 'Sqrt-Hann'
%       * Synthesis window design criterion: Norm minimizing or Distortion minimizing or custom:
%           'Norm_minimizing', 'Distortion_minimizing', or 'Saved_windows'
%
% If the 'PT-WOLA' implementation is selected, the class internally computes either the filter length for subband difference terms and
% filter length for the number of crossband terms based on the chosen plotting mode (filter taps or cross terms). These are stored in
% `diff_filter_length` and `cross_length_matrix`.
%
% Usage:
%   config = PTWOLAConfig();    % Initializes the configuration with default values
%
% Dependencies:
%   None. This class is used by AECParameters.m and associated generalized WOLA processing scripts.
%
% Author: Mohit
% Date: 01/2025
% -------------------------------------------------------------------------
classdef PTWOLAConfig
    properties
        implementation            % 'SWDFT' or 'PT-WOLA'
        Plot_constant             % If implementation -> 'PT-WOLA' then choose 'Const_Cross_Terms' or 'Const_Diff_Terms'
        N                         % Block shift (hop size) in WOLA
        M                         % Window length (usually 2*N)
        algo                      % 'LSQR', 'RLS', 'NLMS'.
        analysis_window           % 'Rectangular', 'Cosine', 'Sqrt-Hann'
        Synthesis_window_def      % 'Distortion_minimizing', 'Norm_minimizing' or any other 'Saved_windows'
        Total_filter              % List of total filter lengths to sweep over
        diff_filter_length        % Derived lengths of difference term filter if total filter and cross-term filter lengths are given. Plot_constant = 'Const_Cross_Terms'.
        cross_length_matrix       % Derived lengths of cross-term filter if total filter and difference-term filter lengths are given. Plot_constant = 'Const_Diff_Terms'.
    end

    methods
        function obj = PTWOLAConfig()
            % Default settings
            obj.implementation = 'SWDFT';                     % 'SWDFT' or 'PT-WOLA'
            obj.Plot_constant = 'Const_Cross_Terms';            % If implementation -> 'PT-WOLA' then choose 'Const_Cross_Terms' or 'Const_Diff_Terms'
            obj.N = 512;
            obj.M = 2 * obj.N;
            obj.algo = 'LSQR';                                  % 'LSQR', 'RLS', or 'NLMS'
            obj.analysis_window = 'Sqrt-Hann';                  % 'Rectangular', 'Cosine', or 'Sqrt-Hann'
            obj.Synthesis_window_def = 'Distortion_minimizing'; % 'Norm_minimizing', 'Distortion_minimizing', or 'Saved_windows'
            obj.Total_filter = 1:4:64;

            % Compute cross term and filter length matrices based on implementation and Plot_constant
            obj = obj.computeCrossLengthMatrix();
        end

        %% Function definition for computeCrossLengthMatrix
        function obj = computeCrossLengthMatrix(obj, R_vals)

            switch obj.implementation
                case 'PT-WOLA'
                    if strcmpi(obj.Plot_constant, 'Const_Diff_Terms')
                        % Subband filter length for PTWOLA with diff terms (column vector)
                        obj.diff_filter_length = reshape([1:2:10, 20:10:60],[],1);
                        % Find crossband filter lengths (R) corresponding to PTEQ filter lengths (nf) and total required filter lengths (T) T = nf+2*R
                        for nf_idx = 1:length(obj.diff_filter_length)
                            % Number of cross-terms for corresponding PTWOLA difference terms
                            obj.cross_length_matrix(nf_idx,:) = ...
                                floor((obj.Total_filter - obj.diff_filter_length(nf_idx)) / 2);
                        end

                    else

                        % For the case of constant difference terms
                        % Default cross-term values if not provided
                        if nargin < 2
                            R_vals = [0, 1, 3, 5, 10, 20:floor(40/8):45, 51];   % Number of one sided cross-terms
                        end
                        % Find PTEQ filter lengths (nf) noncorroding to Crossband filter lengths (R_val) and total required filter lengths (T) T = nf+2*R_val
                        obj.cross_length_matrix = R_vals;
                        for idx = 1:length(R_vals)
                            obj.diff_filter_length(idx,:) = obj.Total_filter - 2 * R_vals(idx);
                        end
                    end
                otherwise
                    % For SWDFT and other types: no cross terms needed
                    obj.cross_length_matrix = [];
                    obj.diff_filter_length = obj.Total_filter;
            end
        end
    end
end