% SignalGenerator.m
% -------------------------------------------------------------------------
% This class generates the set of signals required for simulating Acoustic Echo Cancellation (AEC), including far-end, near-end, echo,
% microphone, and noise signals. The behavior is defined by the configuration provided in an AECParameters object and a Room Impulse Response (RIR).
%
% The generated signals include:
%   - Far-end speaker signal:
%       * 'simulated' – White Gaussian noise shaped with an AR model
%       * 'measured' – Real speech signal resampled to the desired Fs
%   - Near-end signal (only if CDT is enabled):
%       * 'simulated' – AR-filtered noise
%       * 'measured' – Pre-recorded speech signal resampled to the desired Fs
%   - Echo signal – Convolution of the far-end speaker signal with the provided RIR
%   - Microphone signal – Sum of echo, near-end, and added background noise
%
% Usage:
%   sigGen = SignalGenerator(params, rir);     % Create object
%   signals = sigGen.generate();               % Generate signals
%
% Inputs:
%   params – AECParameters object defining simulation settings
%   rir    – Room impulse response vector used for echo modeling
%
% Outputs: signals – Struct containing all relevant simulation signals
%   far_end_speaker – Original far-end excitation signal (played by loudspeaker)
%   far_end         – Echo signal after RIR convolution
%   near_end        – Near-end signal (0 if CDT is false)
%   y_clean         – Clean mic signal (echo + near-end)
%   y               – Microphone signal with added AWGN
%   noise           – Background noise component
%
% Dependencies:
%   - For measured signals, data files must exist in 'auxiliary_data/input_signals/'
%
% Author: Mohit
% Date: 01/2025
% -------------------------------------------------------------------------
classdef SignalGenerator
    properties
        Parameters   % AECParameters instance
        RIR          % Room impulse response (vector)
    end

    methods
        function obj = SignalGenerator(params, rir)
            % Constructor  with params and rir as input
            obj.Parameters = params;
            obj.RIR = rir;
        end

        function signals = generate(obj)
            % Generate far-end, near-end (if CBT flag is true), echo, and mic signal.
            Fs = obj.Parameters.Fs;
            T_sig = obj.Parameters.sig_time;
            sig_len = Fs * T_sig;

            %% Generate Far-End Signal 
            switch lower(obj.Parameters.far_type)
                case 'simulated'
                    % White Gaussian as excitation, shaped with AR model
                    far_end_excitation = randn(sig_len, 1);
                    b = fir1(1024, 0.3);                     % Generate a low pass filter with normalised passband freq.
                    [ar_model, ~] = lpc(b, 0);               % AR_order = 0 ⇒ white noise
                    far_end_speaker = filter(1, ar_model, far_end_excitation);

                case 'measured'
                    % Load real speech and resample
                    data = load('auxiliary_data/input_signals/speech_female_LPP');
                    % Resample the  measured  audio signal to match the required sampling frequency in simulation
                    [P, Q] = rat(Fs / data.Fs);                 % Returns approx num (P) and deno (Q): (fractional) output of a rational number
                    resampled = resample(data.sig_2, P, Q);
                    far_end_speaker = resampled(301:sig_len+300, 1);    % To ignore some silent samples (= 0) in the begining

                otherwise
                    error('Unknown far_type: "%s"', obj.Parameters.far_type);
            end

            %%  Generate Near-End Signal 
            if obj.Parameters.CDT
                switch lower(obj.Parameters.near_type)
                    case 'simulated'
                        near_end_excitation = randn(sig_len, 1);
                        b = fir1(1024, 0.3);                   % LPF
                        [ar_model, ~] = lpc(b, 5);             % AR_order = 5 ⇒ speech-like
                        near_end = filter(1, ar_model, near_end_excitation);
                        near_end(1:Fs * obj.Parameters.silence_period) = 0;

                    case 'measured'
                        data = load('auxiliary_data/input_signals/speech_female_north_wind');
                        % Resample the  measured  audio signal to match the required sampling frequency in simulation
                        [P, Q] = rat(Fs / data.Fs);                     % Returns approx num (P) and deno (Q): (fractional) output of a rational number
                        resampled = resample(data.sig_1, P, Q);
                        near_end = resampled(301:sig_len+300, 1);       % To ignore some silent samples (= 0) in the begining
                        near_end(1:Fs * obj.Parameters.silence_period) = 0;

                    otherwise
                        error('Unknown near_type: "%s"', obj.Parameters.near_type);
                end
            else
                near_end = zeros(sig_len, 1);
            end

            %%  Generate Echo Signal 
            echo = fftfilt(obj.RIR, far_end_speaker);

            %%  EBR Adjustment if CDT flag is true 
            if obj.Parameters.CDT
                scale = db2pow(obj.Parameters.EBR) / (var(echo) / var(near_end));   % Scaling factor to ensure EBR
                far_end_speaker = sqrt(scale) * far_end_speaker;                    % Rescaled far end speaker signal to ensure EBR is acheived (before filtering with RIR)
                echo = fftfilt(obj.RIR, far_end_speaker);                           % Echo signal based on rescalled white noise signal such that the desired EBR is acheived
            end

            %%  Generate Microphone Signal 
            y_clean = echo + near_end;
            y = awgn(y_clean, obj.Parameters.SNR_awgn, 'measured');                 % Add noise to the microphone signal
            noise = y - y_clean;

            %%  Create a structure of relevant signals  
            signals.far_end_speaker = far_end_speaker;
            signals.far_end = echo;
            signals.near_end = near_end;
            signals.y_clean = y_clean;
            signals.y = y;
            signals.noise = noise;

            %%  Logging 
            if obj.Parameters.CDT
                actualEBR = var(echo) / var(near_end);
                fprintf('Signal generated with CDT | EBR = %.2f dB\n', pow2db(actualEBR));
            else
                fprintf('Signal generated without near-end (single-talk)\n');
            end
        end
    end
end