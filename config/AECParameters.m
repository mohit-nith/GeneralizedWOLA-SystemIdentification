% AECParameters.m
% -------------------------------------------------------------------------
% This class defines the parameter configuration for various Acoustic Echo Cancellation (AEC) algorithms used in a generalized
% Weighted OverLap-Add (WOLA) processing framework.
%
% The class defines setup for:
%   - SignalGenerator.m: Signal generation settings (e.g., sampling rate, SNR, echo conditions, Continuous Double Talk (CDT) flags)
%   - RIRGenerator.m: Room Impulse Response (RIR) modeling options
%   - Configuration structures for specific AEC methods:
%       * Generalized (or Per-Tone) WOLA
%       * FDAF (Frequency-Domain Adaptive Filter)
%       * PB-WOLA (Partitioned Block WOLA)
%       * Crossband WOLA
%
% Usage:
%   params = AECParameters();     % Creates an object with default values
%
% Dependencies:
%   PTWOLAConfig.m should be present in the path.
%
% Author: Mohit
% Date: 01/2025
% -------------------------------------------------------------------------

classdef AECParameters
    properties
        % Static parameters
        Fs = 16000                  % Sampling frequency
        sig_time = 20               % Duration in seconds
        SNR_awgn = inf              % SNR of additive white noise
        CDT = false                 % Continuous Double Talk flag
        EBR = 10                    % Echo-to-background ratio in dB (required if CDT = true)
        silence_period = 0          % Duration before near-end starts
        near_type = 'simulated'     % 'simulated' or 'measured'
        far_type = 'simulated'      % 'simulated' or 'measured'
        load_setup = false          % Load saved signal setup containing all observed signals

        RIR = 'windowed_noise'      % RIR generation method - 'RIM', 'windowed_noise' or any saved RIR 'load_RIR'
        RIR_length = 64             % RIR length in samples

        % Dynamic (dependent) parameters
        OutputDir                   % Output path to save results
        PTWOLA                      % Generalized (or Per-tone) WOLA
        % FDAF                      % To be added: FDAF (frequency domain block LMS implementation)
        % PB_WOLA                   % To be added: Partition block WOLA (PB-WOLA) based AEC (without crossterms)  || PB_WOLA filter taps = PTEQ_WOLA filter taps
        % Crossband                 % To be added: Crossband WOLA (System Identification in the STFT Domain With Crossband Filtering, by Y. Avargel and I. Cohen)
    end

    methods
        function obj = AECParameters()
            % Set OutputDir relative to project root (Assuming its one level above the folder of this class)
            thisFile = mfilename('fullpath');     % Full path to AECParameters.m
            configDir = fileparts(thisFile);      % /.../root/ThisFolder
            rootDir = fileparts(configDir);       % /.../root

            obj.OutputDir = fullfile(rootDir, 'results', 'Generalized_WOLA');   % fullfile: Output directory path format immune to Linux or Windows

             % Define PT-WOLA configuration from PTWOLAConfig.m
            obj.PTWOLA = PTWOLAConfig();     
            
            % % FDAF config
            % obj.FDAF.N = obj.RIR_length;            % In FDAF (50% overlap), the block shift is denoted as N and is equal to the size of the RIR in FDAF (careful with PB-FDAF).
            % obj.FDAF.M = 2 * obj.RIR_length;        % Window size is twice the size of RIR in FDAF
            % obj.FDAF.lambda = 0.9;                  % Forgetiing factor to recursively calculate the estimated power of the loudspeaker signal.
            % obj.FDAF.mu = 0.05;                     % Step size for the gradient-descent based algorithms
            % obj.FDAF.alpha = 1e-10;                 % Regularisation term in case denominator of normalisation becomes zero (NLMS)
            % obj.FDAF.AR_order = 0;                  % Order of PEM AR filter
            % obj.FDAF.AR_hop = obj.RIR_length;       % Hopping window size for AR calculation
            % 
            % % PB_WOLA config
            % obj.PB_WOLA.N = obj.RIR_length;         % N denotes the block shift in WOLA (i.e., overlap or downsampling ratio)
            % obj.PB_WOLA.M = 2 * obj.PTWOLA.N;       % Window size (or the DFT size)
            % obj.PB_WOLA.run = false;                % Flag for running PB_WOLA based AEC
            % obj.PB_WOLA.algo = 'RLS';               % 'RLS' or 'NLMS' based adaptive algorithm
            % 
            % % Crossband WOLA config
            % obj.Crossband.N = obj.RIR_length;       % N denotes the block shift in WOLA (i.e., overlap or downsampling ratio)
            % obj.Crossband.M = 2 * obj.PTWOLA.N;     % Window size (or the DFT size)
            % obj.Crossband.run = false;              % Flag for running Crossband WOLA based AEC
            % obj.Crossband.algo = 'RLS';             % 'RLS' or 'NLMS' based adaptive algorithm
            % obj.Crossband.frames = 2;               % Decides total number of frames to be used (like PB_WOLA) (=1 implies no previous frames)
            % obj.Crossband.filter_length = 2;        % Decides number of tones used on each side of the center tones (total filter coeff = frames*(2*filter_length+1))
        end
    end
end