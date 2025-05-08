% computePTWOLAMatrix.m
% -------------------------------------------------------------------------
% This function computes a STFT-domain excitation matrix for low-complexity Per-Tone Weighted Overlap-Add (PT-WOLA) implementation of
% generalized WOLA processing. The excitation matrix is constructed from unwindowed frames of the far-end signal (i.e., prior to room
% impulse response convolution), and includes both frequency-domain cross-terms and time-domain difference terms for each DFT bin.
%
% The excitation matrix includes:
%   - Frequency-domain cross-terms: 
%       * For each DFT bin, 2*R neighboring bins within a ±R range, selected to represent spectral aliasing of windowing operation.
%   - Time-domain difference terms:
%       * nf-1 (real) difference terms, computed from current unwindowed far-end signal frame.
%       * Difference terms are identical across all DFT bins.
%
% Usage:
%   SWDFT_mat = computePTWOLAMatrix(u_block, M, nf, R);
%
% Inputs:
%   u_block – Current unwindowed far-end signal frame of size ([M + nf − 1] × 1)
%   M       – DFT size (assumed even)
%   nf      – Number of filter taps associated with the difference-terms in the adaptive filter
%   R       – Number of filter taps associated with the single-sided cross-terms in the adaptive filter
%
% Output:
%   SWDFT_mat – Excitation matrix of size (M × [nf + 2*R]), where each row (DFT bin) contains:
%               2*R cross-terms, nf-1 difference terms, and 1 current frame DFT term.
%
% Author: Mohit
% Date: 01/2025
% -------------------------------------------------------------------------
function SWDFT_mat = computePTWOLAMatrix(u_block, M, nf, R)

% Total number of excitation terms per DFT bin: nf difference + 2*R cross terms
nf_total = nf + 2 * R;

% Preallocate output excitation matrix: [#DFT bins × #excitation terms]
SWDFT_mat = zeros(M, nf_total);

% Compute normalized DFT of the current unwindowed frame (length M)
unwindowed_DFT = (1 / sqrt(M)) * fft(u_block(nf : nf + M - 1));

% Repetition ensures circular access when selecting ±R neighboring bins
DFT_vec_rep = repmat(1:M, 1, 3);

%% Generate frequency-domain cross-terms
for dft_idx = 1 : M / 2 + 1
    % Select indices for ±R neighbors centered at current bin
    bins_idx = M + dft_idx + (-R : R);    
    % Assign neighboring DFT bin values (cross-terms) to corresponding row
    SWDFT_mat(dft_idx, 1 : 2 * R + 1) = unwindowed_DFT(DFT_vec_rep(bins_idx));       
end

% Fill mirrored bins
SWDFT_mat(M/2+2:end, :) = conj(flipud(SWDFT_mat(2:M/2, :)));

%% Generate difference terms
Diff_mat = [-eye(nf-1), zeros(nf-1, M - nf + 1), eye(nf-1)];
Diff_terms = flipud(Diff_mat * u_block(:));

% Repeat difference terms across all frequency bins (rows)
SWDFT_mat(:, 2 * R + 2:end) = repmat(Diff_terms(:).', M, 1);
end