% Generalized_WOLA_AEC_Processor
% -------------------------------------------------------------------------
% Core engine that performs adaptive filtering using Generalized WOLA. Supports:
%   - SWDFT or PT-WOLA implementation
%   - LSQR (batch) or RLS (block / online)
%
% Used by AECSimulation.m for executing subband adaptive AEC.
% Inputs:
%   y                       – Microphone signal (near + echo + noise)
%   far_end_exctn           – Loudspeaker excitation (pre-RIR)
%   far_end                 – Echo signal (far-end filtered with RIR)
%   nf                      – Number of difference-term filter taps
%   R                       – Number of cross terms (on one side) filter taps
%   M                       – DFT-size or window length (assumed same)
%   analysis_window         – Type of analysis window (e.g., 'Sqrt-Hann')
%   Synthesis_window_def    – Synthesis window type: 'Norm_minimizing', etc.
%   WOLA_Implementation     – Specifies 'pt-wola' or 'swdft' implementation for Generalized WOLA
%   AdaptAlgo               – 'LSQR' or 'RLS'
%   anaWin                  – (Optional) analysis window vector
%   synWin                  – (Optional) synthesis window vector
%
% Outputs:
%   ERLE_mean               – Mean Echo Return Loss Enhancement (dB)
%   Res_error_mean          – Mean power of residual error in each frequency bin (dB)
%   F_est                   – Estimated subband filter coefficients
%
% Author: Mohit
% Last updated: 01/2025
% -------------------------------------------------------------------------
classdef Generalized_WOLA_AEC_Processor
    methods (Static)
        function [ERLE_mean, Res_error_mean, F_est] = run(y,far_end_exctn,far_end,nf,R,M,...
                                                analysis_window,Synthesis_window_def,WOLA_Implementation,AdaptAlgo,anaWin,synWin)
%% Parameters and Setup
nf_total = nf + 2 * R;          % Total number of subband filter taps (diff + cross terms)
OL = 0.5;                       % Overlap factor (Not all OL guarantee perfect reconstruction). OL=2 been shown to work well in practice
N = round((1-OL)*M);            % Redefine the value of hop size based on the overlapping factor (in case (1-OL)*M is non-integer)
scaling_fact = 2*(1-OL);        % If windows are designed, should be 1 irrespective of OL (window design algorithm takes care of scaling)
pad_length = round(1 / (1 - OL) - 1) * N;

% Zero padding to account for delay/overlap
far_end = [zeros(pad_length, 1); far_end];
y = [zeros(pad_length, 1); y];
far_end_exctn = [zeros(pad_length + nf - 1, 1); far_end_exctn];

num_frames = floor((length(far_end) - M) / N);

% Generate analysis and synthesis windows if not provided
if isempty(anaWin) || isempty(synWin)
    % Generate a structure to pass to WindowDesigner class
    cfg.M = M;
    cfg.N = N;  % Assuming 50% overlap
    cfg.analysis_window = analysis_window;
    cfg.Synthesis_window_def = Synthesis_window_def;

    % Call WindowDesigner class to generate windows
    winGen = WindowDesigner(cfg);
    [anaWin, synWin] = winGen.generate();
end
ana_loud = anaWin;  % In case analysis window for microphone and loudspeaker are different, define the window for loudspeaker here

% Memory Preallocation
Y_STFT      = zeros(M, num_frames);
Echo_STFT   = zeros(M, num_frames);             % STFT domain echo signal  
Far_SWDFT   = zeros(M, nf_total, num_frames);   % STFT domain Far end excitation-signal before RIR 


% Preallocate adaptive filter output and define hyperparameters
if strcmpi(AdaptAlgo, 'LSQR')
    F_est = zeros(M, nf_total);
elseif strcmpi(AdaptAlgo, 'RLS')
    F_est = zeros(M, nf_total, num_frames + 1);
    delta =1e-12;                                   % Hyperparameter. Value depends on SNR: 1e-12 for cos, norm_min. nf=3
    P = repmat(eye(nf_total)/delta, [1, 1, M]);     % Initial estimate of the inverse (of) correlation matrix
    forget_factor = 1;                              % Forgetting factor
end


%% Starting the loop for the Generalized WOLA based AEC

% Stage 1: Analysis Filter Bank
for frame_idx = 1:num_frames
    % Extract current block
    currBlk = (frame_idx - 1) * N + (1:M);
    y_block = y(currBlk);                                           % Microphone signal (far_end + noise + near_end)
    u_block = far_end_exctn(currBlk(1):(currBlk(end) + nf - 1));    % Loudspeaker signal before RIR
    echo_block = far_end(currBlk);                                  % Echo signal: Loudspeaker signal after RIR


    % Analysis Filter Bank: Microphone and Far end 
    Y_STFT(:,frame_idx) = (1 / sqrt(M)) * fft(anaWin .* y_block);
    Echo_STFT(:,frame_idx) = (1 / sqrt(M)) * fft(anaWin .* echo_block);

    % Analysis Filter Bank: Loudspeaker
    switch lower(WOLA_Implementation)
        case 'pt-wola'
            Far_SWDFT(:, :, frame_idx) = computePTWOLAMatrix(u_block, M, nf, R);
        case 'swdft'
            Far_SWDFT(:, :, frame_idx) = computeSWDFTMatrix(u_block, ana_loud, M, nf);
        otherwise
            error('Unsupported Generalized WOLA Implementation: %s', WOLA_Implementation);
    end
end
    
% Stage 2: Subband Adaptive Filtering
switch lower(AdaptAlgo)
    % LSQR (batch processing)
    case 'lsqr'
        for freq_idx = 1:M/2+1
            U = permute(Far_SWDFT(freq_idx, :, :), [3, 2, 1]);    % Reshape Far_SWDFT to [num_frames x nf_total] matrix
            d = Y_STFT(freq_idx, :).';                                  % Desired signal vector for freq_idx [frames x 1]
            % Solve the least-squares system
            F_est(freq_idx, :) = conj(U \ d);
        end
        % Reconstruct conjugate-symmetric part of filter
        F_est(M/2+2:end, :) = conj(flipud(F_est(2:M/2, :)));

    case 'rls'
        % RLS (block processing)
        for frame_idx = 1:num_frames
            [F_est(:, :, frame_idx + 1), P] =   updateRLS(Y_STFT(:, frame_idx), squeeze(Far_SWDFT(:, :, frame_idx)), ...
                                                F_est(:, :, frame_idx), P, M, forget_factor);
        end
    otherwise
        error('Unsupported algorithm: %s', AdaptAlgo);
end


% Subband Echo Estimation
Echo_est_STFT = conj(squeeze(sum(F_est .* conj(Far_SWDFT), 2)));
Subband_res_err  = Echo_STFT - Echo_est_STFT;

% Stage 3: Reconstruction with Synthesis Filter Bank
echo_recon      = reconstructTD_SynthesisFB(Echo_est_STFT, synWin, N, M, num_frames);
residual_recon  = reconstructTD_SynthesisFB(Subband_res_err, synWin, N, M, num_frames);

% Performance Metrics
[Res_error_mean, ERLE_mean] = computeERLE(far_end, residual_recon, N, pad_length, scaling_fact, num_frames);

        end
    end
end