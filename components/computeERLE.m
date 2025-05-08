% computeERLE.m
% -------------------------------------------------------------------------
% Computes frame-wise and average ERLE (Echo Return Loss Enhancement) from
% clean echo and residual echo signals in the time domain.
%
% Inputs:
%   clean_echo   – Echo signal before cancellation (far-end passed through RIR)
%   res_echo     – Residual echo signal after cancellation
%   N            – Hop size 
%   pad_length   – Number of initial samples to discard (zero-padding offset)
%   num_frames   – Total number of STFT frames used by WOLA filter bank
%
% Outputs:
%   res_error_db – Frame-wise ERLE values for convergence behavior  (in dB)
%   ERLE_mean    – Mean ERLE across the signal duration (in dB)
% -------------------------------------------------------------------------
function [ERLE_frame, ERLE_mean] = computeERLE(Clean_echo, Res_echo, N, pad_length, scaling_fact, num_frames)
% Remove zero-padding from the beginning
clean_echo = Clean_echo(pad_length + 1:end);
res_echo   = scaling_fact * real(Res_echo(pad_length + 1:end));

ERLE_frame = zeros(num_frames, 1);

for frame_idx = 1:num_frames
    currBlk = (frame_idx - 1) * N + (1:N);
    num = sum(clean_echo(currBlk).^2);
    den = sum(res_echo(currBlk).^2 + eps);  % eps to avoid division by 0
    ERLE_frame(frame_idx) = 10 * log10(num / den);
end

% Mean ERLE over all samples (not frame-averaged)
ERLE_mean = 10 * log10(sum(clean_echo.^2 + eps) / sum(res_echo.^2 + eps));
end