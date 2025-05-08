% demo_design_wola_synthesis_window.m
% -------------------------------------------------------------------------
% This script demonstrates the use of design_wola_synthesis_window.m. It creates an analysis prototype window and then designs a
% synthesis window for a WOLA filter bank under different synthesis window design strategies: 
% 'distortion_minimizing' and 'norm_minimizing'.
% The resulting synthesis windows are plotted for comparison.
% NOTE: Assumes 50% overlap
% Author: Mohit Sharma 
% -------------------------------------------------------------------------

clear; clc; close all;
%% Parameters
M = 1024;                               % Length of the analysis window
% Create a prototype analysis window
anaWin = sqrt(hann(M, 'periodic'));     % OR cosine window: hann(M, 'periodic');, rectangular window: sqrt(0.5) * ones(M, 1)

%% Design Synthesis Windows
% Mode: 'distortion_minimizing'
syn_win_distort = design_wola_synthesis_window(anaWin, 'distortion_minimizing');

% Mode: 'norm_minimizing'
syn_win_norm = design_wola_synthesis_window(anaWin, 'norm_minimizing');

%% Plot the Results
figure;
subplot(3,1,1);
plot(anaWin, 'LineWidth', 2);
xlim([1 M]);
title('Analysis Window');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(syn_win_distort, 'r', 'LineWidth', 2);
title('Synthesis Window (Distortion Minimizing)');
xlabel('Sample Index');
ylabel('Amplitude');
xlim([1 M]);
grid on;

subplot(3,1,3);
plot(syn_win_norm, 'b', 'LineWidth', 2);
title('Synthesis Window (Norm Minimizing)');
xlabel('Sample Index');
ylabel('Amplitude');
xlim([1 M]);
grid on;

sgtitle('WOLA Synthesis Window Design Comparison');

%% Evaluate Perfect Reconstruction (PR) Constraint
D = M / 2;
PR_error_distort = abs(anaWin(1:D) .* syn_win_distort(1:D) + anaWin(D+1:end) .* syn_win_distort(D+1:end) - 1);
PR_error_norm = abs(anaWin(1:D) .* syn_win_norm(1:D) + anaWin(D+1:end) .* syn_win_norm(D+1:end) - 1);
                  
fprintf('Max PR error (Distortion Minimizing): %e\n', max(PR_error_distort));
fprintf('Max PR error (Norm Minimizing): %e\n', max(PR_error_norm));