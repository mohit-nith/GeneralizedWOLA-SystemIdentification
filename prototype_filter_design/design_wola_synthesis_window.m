% design_wola_synthesis_window.m
% -------------------------------------------------------------------------
% Designs a prototype synthesis window for a Weighted Overlap-Add (WOLA) filter bank using convex optimization based on the framework
% presented in:
%       M. Sharma and M. Moonen, "Prototype filter design for weighted overlap-add filter bank based sub-band adaptive filtering
%           applications," 2023 31st European Signal Processing Conference (EUSIPCO), Helsinki, Finland, 2023, pp. 366-370.
%
% Synthesis window design options:
%   * 'Distortion_minimizing'   – Synthesis window design with minimum signal distortion.
%   * 'Norm_minimizing'         – Conventional norm minimizing synthesis window design
%
% Inputs:
%   ana_win : [M x 1] analysis prototype window
%   syn_def : (optional) string, with synthesis window design strategy
%             'distortion_minimizing' (default) – Eq. 21 in paper
%             'norm_minimizing'                 – minimize norm only
%
% Output:
%   syn_win : [M x 1] synthesis window (i.e., prototype synthesis filter for DFT-modulated FB)
%
% NOTE:
% This implementation currently supports only distortion minimization (g₀) via real-valued operations. Alias term support will be added in future.
% Author: Mohit
% Date: 08/2023
% -------------------------------------------------------------------------
function synWin = design_wola_synthesis_window(anaWin, syn_def)

% Define default behavior
if nargin < 2
    syn_def = 'distortion_minimizing';
end

% Fix input formatting
anaWin = anaWin(:);

% Setup Parameters
M = length(anaWin);
D = M / 2;      % Assumes 50% overlap
dist_len = D;
numAlias = 0;   % Aliasing terms to consider (g_1, g_2, ... g_numAlias)
                %Note:  for numAlias>0: Error: In MATLAB Optimization variables cannot be combined with complex data. (use CVX to define
                % problem). Will be updated.

if numAlias > 0
    warning(['Aliasing terms are not supported in this version of the function due to complex-valued ' ...
        'optimization variable limitations in MatLab. Matlab error: Error using optim.problemdef.OptimizationNumeric' ...
        ' Optimization variables cannot be combined with complex data. Use CVX for full aliasing-aware design.']);
end

% Weight selection based on syn_def
switch lower(syn_def)
    case 'distortion_minimizing'    % w_distort should be >> w_alias for distortion_minimizing. Ideal w_distort = 1
        w_alias = 0;
        w_distort = 1 - w_alias;
        w_regularize = 120;         % Weight for regularization term
    case 'norm_minimizing'
        w_alias = 0;
        w_distort = 0;
        w_regularize = 1;           % Weight for regularization term
    otherwise
        error('Unsupported syn_def "%s". Use "distortion_minimizing" or "norm_minimizing".', syn_def);
end

% Define Optimization Variable
f = optimvar('f', M, 1, 'LowerBound', -Inf);  % Synthesis window

% Generate modulated analysis filters h_n and their Toeplitz convolution matrices H_n
Hn = cell(numAlias + 1, 1);   % H0 = distortion, Hn (n>0) = alias
g = cell(numAlias + 1, 1);    % g{n} = H_n * f

for shift_idx = 0:numAlias
    hMod = anaWin(:).' .* exp(1i * 2 * pi * (0:M-1) * shift_idx / D);
    Hn{shift_idx+1} = toeplitz([hMod(:); zeros(M - 1, 1)], [hMod(1); zeros(M - 1, 1)]);
    g{shift_idx+1} = Hn{shift_idx+1} * f;  % ERROR for numAlias>0: In MATLAB Optimization variables cannot be combined with complex data.
    % Workaround: Use CVX solver https://cvxr.com/cvx/ for optimization problem (to be added)
end


% Objective function: Combine distortion, aliasing, and regularization
distortion_term = norm(g{1}(1:dist_len)) + norm(g{1}(2*D:2*D+dist_len-1) - D);

alias_term = 0;
for i = 2:numAlias+1
    alias_term = alias_term + norm(g{i}(2*D+1:3*D));
end

regularization_term = norm(f);

% Initialize optimization problem
prob = optimproblem('ObjectiveSense','min');

% Define objective function
prob.Objective  =   w_distort * distortion_term + ...
    w_regularize * regularization_term + ...
    w_alias * alias_term;

% Define Perfect Reconstruction Constraints
ceq = optimexpr(1, D);
for win_idx = 1:D
    ceq(win_idx) =  anaWin(win_idx)*f(win_idx) + ...
        anaWin(D+win_idx)*f(D+win_idx) - 1;
end
prob.Constraints.window_constraints = ceq == 0;

% Solve using  optimization toolbox
x0.f = ones(M, 1);  % or any other initial guess of size N×1
opts = optimoptions('fmincon', 'Display', 'off');
sol = solve(prob, x0, 'Options', opts);

% Extract output
synWin = sol.f;

end