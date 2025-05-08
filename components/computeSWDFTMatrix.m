% computeSWDFTMatrix.m
% -------------------------------------------------------------------------
% This function computes the STFT-domain excitation matrix for generalized WOLA processing using the full-complexity Sliding Windowed
% DFT (SWDFT) approach. The matrix is from a window applied to each sliding segment of the far-end signal, followed by a normalized DFT operation.
%
% The excitation matrix contains:
%   - nf windowed DFT vectors corresponding to nf sliding blocks of the current far-end signal frame.
%   - Each column of the excitation matrix corresponds to the M-point DFT of one sliding windowed segment.
%
% Usage:
%   SWDFT_mat = computeSWDFTMatrix(u_block, ana_loud, M, nf);
%
% Inputs:
%   u_block   – Current unwindowed segment of the far-end signal frame of size ([M + nf − 1] × 1)
%   ana_loud  – Analysis window of length M applied to each sliding segment
%   M         – DFT size (assumed even)
%   nf        – Number of adaptive filter taps (= the number of sliding DFT blocks)
%
% Output:
%   SWDFT_mat – STFT-domain excitation matrix of size (M × nf), where each column is a DFT of a windowed segment
%
% Author: Mohit
% Date: 01/2025
% -------------------------------------------------------------------------
function SWDFT_mat = computeSWDFTMatrix(u_block, ana_loud, M, nf)

% Preallocate excitation matrix: [total DFT bins × number of filter taps]
SWDFT_mat = zeros(M, nf);
for sliding_idx = 1:nf
    % Compute the offset to select current sliding window samples
    offset = nf - sliding_idx;

    % Apply analysis window to the selected current sliding block of samples
    u_win = ana_loud .* u_block(offset + (1:M));

    % Compute normalized DFT and store in corresponding column
    SWDFT_mat(:, sliding_idx) = (1 / sqrt(M)) * fft(u_win);
end
end