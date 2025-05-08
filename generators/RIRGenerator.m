% RIRGenerator.m
% -------------------------------------------------------------------------
% This class handles Room Impulse Response (RIR) generation based on the configuration specified in the AECParameters object.
% Supports multiple RIR modeling methods to simulate acoustic environment for echo cancellation.
%
% The supported RIR generation methods are:
%   - 'rim' – Generates impulse response of the (randomised) image method, proposed in De Sena et al. "On the modeling of rectangular
%       geometries in  room acoustic simulations." IEEE/ACM Transactions on Audio, Speech and Language Processing (TASLP) 23.4 (2015): 774-786.
%   - 'windowed_noise' – Generates an artificial RIR using exponentially decaying white noise.
%   - 'load_rir' – Load a precomputed or measured RIR from a file.
%
% Usage:
%   rirGen = RIRGenerator(params);    % Create object with given parameters
%   rir = rirGen.generate();          % Generate and return the RIR
%
% Inputs:
%   params – An instance of the AECParameters class containing RIR settings
%
% Outputs:
%   rir – A vector representing the generated or loaded room impulse response
%
% Dependencies:
%   - For 'rim' method: Requires `rim()` function in the path
%   - For 'load_rir': Define the file name with a variable `RIR` in folder 'auxiliary_data/rir/'
%
% Author: Mohit
% Date: 01/2025
% -------------------------------------------------------------------------
classdef RIRGenerator
    properties
        Parameters  % AECParameters
    end

    methods
        function obj = RIRGenerator(params)
            % Constructor with params as input
            obj.Parameters = params;
        end

        function rir = generate(obj)
            % Generate the RIR based on the method specified in Parameters.RIR
            Fs = obj.Parameters.Fs;
            RIR_length = obj.Parameters.RIR_length;
            method = lower(obj.Parameters.RIR);

            switch method
                case 'rim'
                    % Randomized Image Method (RIM) based RIR ("On the modeling of rectangular geometries in  room acoustic simulations.")
                    RIM.Lx = 5; RIM.Ly = 5; RIM.Lz = 3;           % Room dimensions
                    RIM.beta = 0.05 * ones(2, 3);                 % Reflection coefficients
                    RIM.rir_length = RIR_length / Fs;             % RIR length in seconds
                    RIM.rand_dist = 0.08;                         % Random offset for positions
                    RIM.Tw = 40 / Fs;                             % Fractional delay filter length
                    RIM.Fc = 0.9 * (Fs / 2);                      % Fractional delay cutoff frequency
                    RIM.c = 343;                                  % Speed of sound
                    RIM.T60 = 0.16 * RIM.Lx * RIM.Ly * RIM.Lz / ...
                        (2 * (RIM.Lx*RIM.Ly + RIM.Lz*RIM.Ly + RIM.Lx*RIM.Lz) * ...
                        (1 - mean(RIM.beta, 'all')^2));           % Reverberation time estimate
                    RIM.lpk_pos = [2.5; 1; 1.5];                  % Loudspeaker position
                    RIM.mic_pos = [2; 2.5; 1.5];                  % Microphone position

                    % Generate RIR calling rim() function
                    rir = rim(RIM.mic_pos, RIM.lpk_pos, [RIM.Lx, RIM.Ly, RIM.Lz], RIM.beta, RIM.rir_length, Fs, RIM.rand_dist, ...
                        RIM.Tw, RIM.Fc, RIM.c);

                case 'windowed_noise'
                    % Exponentially decaying windowed white noise RIR
                    if RIR_length>100
                        beta = 50;              % Taper length to ensure smooth declining RIR tail (for long RIRs)
                        alpha_RIR = 0.006;      % Exponential decay rate
                    else
                        beta=0;
                        alpha_RIR = 0.004;      % Exponential decay rate
                    end
                    % Generate exponentially decaying noise
                    rir = wgn(RIR_length, 1, 0) .* exp(-alpha_RIR * (0:RIR_length-1).');
                    % Apply cosine window taper to the tail of the RIR
                    cos_window = 0.5 * (1 + cos(pi * (1/2:beta) / beta));
                    rir(end - beta + 1:end) = rir(end - beta + 1:end) .* cos_window';
                    rir = rir / norm(rir);      %  Normalize the RIR

                case 'load_rir'
                    % Load measured or generated RIR from a file
                    if isfile('auxiliary_data/rir/RIR_HA_smooth.mat')
                        data = load('RIR_HA_smooth.mat');
                        if isfield(data, 'RIR')
                            rir = data.RIR;
                        else
                            error('RIR variable not found in file.');
                        end
                    else
                        error('File RIR_HA_smooth.mat not found.');
                    end

                otherwise
                    error('Unsupported RIR method: "%s"', obj.Parameters.RIR);
            end

            fprintf('RIR generated using method: %s\n', lower(method));
        end
    end
end