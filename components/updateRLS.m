% updateRLS.m
% -------------------------------------------------------------------------
% This function performs frame-wise update of subband adaptive filter coefficients in STFT domain using recursive least squares (RLS)
% algorithm. The RLS algorithm is applied independently to each frequency bin and uses an exponentially weighted error criterion to
% recursively minimize the a priori error between the observed microphone signal and the adaptive filter output.
%
% The RLS update includes:
%   - Computation of the RLS (or Kalman) gain vector per frequency bin
%   - Subband filter coefficient update using a priori error
%   - Recursive update of the inverse correlation matrix (P_prev matrix) using the matrix inversion lemma
%
% Usage:
%   [F_next, P_updated] = updateRLS(Y_bin, U_bin, F_prev, P_prev, M, forget_factor);
%
% Inputs:
%   Y_bin         – STFT-domain microphone signal vector (M × 1)
%   U_bin         – Excitation signal STFT matrix (M × nf_total)
%   F_prev        – Previous subband filter coefficient estimates (M × nf_total)
%   P_prev        – Prior estimate of inverse correlation matrices (per bin) (M × nf_total × nf_total)
%   M             – DFT size
%   forget_factor – Forgetting factor for RLS (typically close to 1). Controls memory depth and adaptation speed
%
% Outputs:
%   F_next     – Updated subband filter coefficient estimates (M × nf_total)
%   P_updated  – Posterior (updated) inverse correlation matrices (per bin) (M × nf_total × nf_total)
%
% Author: Mohit
% Date: 01/2025
% -------------------------------------------------------------------------
function [F_next, P_updated] = updateRLS(Y_bin, U_bin, F_prev, P_prev, M, forget_factor)

% Initialize updated filter and correlation matrix with previous values
F_next = F_prev;
P_updated = P_prev;

% Loop over DFT bins: DC to Nyquist
for freq_idx = 1:M/2+1
    % Extract current excitation vector for current frequency bin (nf_total x 1)
    u = squeeze(U_bin(freq_idx, :)).';

    % Extract desired signal (microphone) DFT coefficient at current frequency bin
    d = conj(Y_bin(freq_idx));

    % Project current excitation vector onto inverse correlation matrix
    Pi = squeeze(P_prev(freq_idx, :, :)) * u;       % Intermediate vector for gain computation

    % Compute RLS (or Kalman) gain vector for current frequency bin
    k = Pi / (forget_factor + u' * Pi);

    % Compute a priori error at current frequency bin
    e = d - u' * F_prev(freq_idx, :).';

    % Update filter coefficients for current frequency bin
    F_next(freq_idx, :) = F_prev(freq_idx, :) + (k * conj(e)).';

    % Update inverse correlation matrix using matrix inversion lemma for current frequency bin
    P_updated(freq_idx, :, :) = (1/forget_factor) * (P_prev(freq_idx, :, :) - k * u' * P_prev(freq_idx, :, :));
end

% Reconstruct the conjugate-symmetric half of DFT bins  (Since all time-domain signals are real-valued)
F_next(M/2+2:end, :) = conj(flipud(F_next(2:M/2, :)));
end