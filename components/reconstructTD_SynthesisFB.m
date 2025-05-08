% reconstructTD_SynthesisFB.m
% -------------------------------------------------------------------------
% This function reconstructs the time-domain signal from its STFT-domain representation using the synthesis filter bank in WOLA
% framework (assumes 50% overlap).
%
% The reconstruction involves:
%   - Inverse DFT (IDFT) of each complex STFT frame
%   - Application of a synthesis window to each frame
%   - Overlapping and adding the windowed frames at hop size (N) intervals to form the time-domain signal
%
% Usage:
%   y_TD = reconstructTD_SynthesisFB(STFT_data, syn_win, N, M, num_frames);
%
% Inputs:
%   STFT_data   – Complex STFT-domain matrix of size (M × num_frames), each column is a frequency-domain frame
%   syn_win     – Synthesis window of length M
%   N           – Hop size between successive frames (e.g., M/2 for 50% overlap)
%   M           – DFT size (assumed equal to the window length)
%   num_frames  – Total number of STFT frames
%
% Output:
%   y_TD        – Reconstructed time-domain signal of length (N × num_frames + M)
%
% Author: Mohit
% Date: 01/2025
% -------------------------------------------------------------------------
function y_TD = reconstructTD_SynthesisFB(STFT_data, syn_win, N, M, num_frames)

% Preallocate output signal buffer
y_TD = zeros(N * num_frames + M, 1);

% Loop over each frame in the STFT-domain data
for frame_idx = 1:num_frames
    % Compute sample indices for the current output frame position
    currBlk = (frame_idx - 1) * N + (1:M);

    % Compute inverse DFT of the current frame
    ISTFT_signal = sqrt(M) * ifft(STFT_data(:, frame_idx), 'symmetric');


    % Apply synthesis window and accumulate via overlap-add
    y_TD(currBlk) = y_TD(currBlk) + syn_win .* ISTFT_signal;
end