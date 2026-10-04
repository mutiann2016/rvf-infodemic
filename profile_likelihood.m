function profile_likelihood
% Compute 95% profile-likelihood intervals for the three estimated
% parameters sigma1, BetaA, phiF, using the reduced objective of
% Section 2.3 of the manuscript.
%
% For each parameter p_j, we grid-fix p_j over a range around its
% fitted value, re-optimise the other two parameters by fmincon, and
% record the minimum objective value J_j(v) at each grid point.
% The 95% interval is the set of v for which
%     J_j(v) - J_opt <= chi2inv(0.95, 1) = 3.8415.
%
% The output is saved to 'RVF_profile_output.txt' and printed to the
% Command Window.

clc; close all;

fprintf('\n=== Profile likelihood for the reduced fit ===\n');
fprintf('Started: %s\n\n', datestr(now));

% ---------- Load literature parameters and set globals ----------
Pardef

global Z0 tt
global LambdaH LambdaL LambdaV muH muL muV muM
global gammaH deltaH gammaL epslonL epslonV RoV
global sigma0 Phi0 Phi1 BetaHV BetaLV BetaVL Kappa2 Kappa4 mu2
global theta BetaHL BetaR bmax Kappa1 Kappa3 mu1 mu3 mu4 gamma Alpha
global sigma1 BetaA phiF

% Fixed parameters (Table 3 of the manuscript)
theta  = 0.47;
BetaHL = 0.05;
BetaR  = 0.33;
bmax   = 182;
Kappa1 = 0.29;
Kappa3 = 0.21;
mu1    = 1.7e-4;
mu3    = 0.38;
mu4    = 0.26;
gamma  = 0.031;
Alpha  = 0.63;

% Fitted parameters (Table 5 of the manuscript)
q_opt = [0.954; 1.159; 1.365];

% ---------- Initial condition and time grid (matches main.m) ----------
Astar = 0.617;
S0H0  = LambdaH / (muH + q_opt(2) * Astar);
S1H0  = (q_opt(2) * Astar / muH) * S0H0;
SL0   = LambdaL / muL;
EL0   = 1;
IL0   = 0;
RL0   = 0;
A0    = Astar;
R0v   = 0;
SV0   = LambdaV / muV;
EV0   = 0;
IV0   = 0;
MV0   = 0;
Z0 = [SL0; EL0; IL0; RL0; S0H0; S1H0; 0; 0; A0; R0v; SV0; EV0; IV0; MV0];

T    = 0.5;
Nmax = 500;
tt   = 0 : T/Nmax : T;

% ---------- Bounds (same as main.m) ----------
LB = [1e-3; 1e-2; 0];
UB = [5;    2;    3];

% ---------- Optimiser options ----------
OPTIONS = optimset('fmincon');
OPTIONS = optimset(OPTIONS, 'Display','off', ...
    'MaxFunEvals', 3000, 'MaxIter', 300, ...
    'TolFun', 1e-5, 'TolX', 1e-5);

% ---------- Evaluate objective at the optimum ----------
J_opt = objective(q_opt);
fprintf('Objective at fitted optimum: J_opt = %.6g\n', J_opt);
fprintf('Threshold for 95%% profile interval: J_opt + 3.8415 = %.6g\n\n', ...
        J_opt + chi2inv(0.95, 1));

threshold = J_opt + chi2inv(0.95, 1);

% ---------- Parameter names and grid half-widths ----------
names     = {'sigma1','BetaA','phiF'};
halfwidth = [0.30, 0.30, 0.30];   % search +/- 30% of the fitted value

% ---------- Results storage ----------
results = struct();

for j = 1:3
    fprintf('--- Profiling parameter %d: %s ---\n', j, names{j});

    p_opt  = q_opt(j);
    lo_v   = max(LB(j), p_opt * (1 - halfwidth(j)));
    hi_v   = min(UB(j), p_opt * (1 + halfwidth(j)));

    grid_v = linspace(lo_v, hi_v, 41);
    J_grid = zeros(size(grid_v));

    for k = 1:numel(grid_v)
        v  = grid_v(k);
        q0 = q_opt;
        q0(j) = v;

        % Optimise the other two parameters with p_j fixed at v
        idx_free = setdiff(1:3, j);
        objfix   = @(qfree) objective_set_free(q_opt, j, v, qfree, idx_free);

        try
            [~, Jk] = fmincon(objfix, q_opt(idx_free), [],[],[],[], ...
                              LB(idx_free), UB(idx_free), [], OPTIONS);
        catch
            Jk = Inf;
        end
        J_grid(k) = Jk;
    end

    % Find the 95% interval
    in_interval = J_grid <= threshold;
    if any(in_interval)
        v_lo = min(grid_v(in_interval));
        v_hi = max(grid_v(in_interval));
    else
        v_lo = NaN;
        v_hi = NaN;
    end

    fprintf('  Fitted value        : %.6g\n', p_opt);
    fprintf('  95%% profile interval: [%.4g, %.4g]\n\n', v_lo, v_hi);

    results(j).name    = names{j};
    results(j).value   = p_opt;
    results(j).lo      = v_lo;
    results(j).hi      = v_hi;
    results(j).grid_v  = grid_v;
    results(j).grid_J  = J_grid;

    % Plot the profile
    fig = figure('Units','pixels','Position',[100+j*40 100 800 500]);
    plot(grid_v, J_grid, 'b-', 'LineWidth', 1.6); hold on;
    yline(threshold, 'r--', 'LineWidth', 1.3);
    yline(J_opt, 'k:', 'LineWidth', 1.1);
    if ~isnan(v_lo)
        xline(v_lo, 'g--', 'LineWidth', 1.3);
        xline(v_hi, 'g--', 'LineWidth', 1.3);
    end
    xlabel(sprintf('%s', names{j}), 'Interpreter', 'none');
    ylabel('Objective J(u)');
    title(sprintf('Profile likelihood for %s', names{j}));
    grid on;
    legend('Profile J(v)','95% threshold','J_{opt}','Interval bounds', ...
           'Location', 'best');

    pngfile = fullfile(pwd, sprintf('profile_%s.png', names{j}));
    saveas(fig, pngfile);
    fprintf('  Saved profile plot to %s\n\n', pngfile);
end

% ---------- Save the results to text file ----------
outfile = fullfile(pwd, 'RVF_profile_output.txt');
fid = fopen(outfile, 'w');
fprintf(fid, 'Profile-likelihood 95%% intervals for the reduced fit\n');
fprintf(fid, 'Run at %s\n\n', datestr(now));
fprintf(fid, 'J_opt      = %.6g\n', J_opt);
fprintf(fid, 'Threshold  = %.6g\n\n', threshold);
fprintf(fid, '%-10s %12s %12s %12s\n', ...
        'Parameter', 'Value', 'Lower 95%', 'Upper 95%');
for j = 1:3
    fprintf(fid, '%-10s %12.4g %12.4g %12.4g\n', ...
            results(j).name, results(j).value, results(j).lo, results(j).hi);
end
fclose(fid);

fprintf('=== Profile likelihood complete ===\n');
fprintf('Results saved to: %s\n', outfile);
fprintf('Finished: %s\n\n', datestr(now));

end

% =============================================================
%  Helper: set one parameter to a fixed value and the other two
%  to the free vector, then evaluate the objective.
% =============================================================
function J = objective_set_free(q_opt, j_fixed, v_fixed, qfree, idx_free)
    q = q_opt;
    q(j_fixed) = v_fixed;
    q(idx_free) = qfree;
    J = objective(q);
end