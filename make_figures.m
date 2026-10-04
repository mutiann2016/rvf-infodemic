function make_figures
% Regenerate every numerical figure referenced in the manuscript.
% Saves PNG files to a 'figs/' subfolder next to the working directory.
% Uses subplot for MATLAB R2018b and earlier compatibility.
% Shading is drawn first so the data lines appear on top.

clc; close all;

% ---------- Output folder ----------
outdir = fullfile(pwd, 'figs');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
fprintf('Saving figures to: %s\n\n', outdir);

% ---------- Load literature parameters ----------
Pardef

global Z0 tt
global LambdaH LambdaL LambdaV muH muL muV muM
global gammaH deltaH gammaL epslonL epslonV RoV
global sigma0 Phi0 Phi1 BetaHV BetaLV BetaVL Kappa2 Kappa4 mu2

% ---------- Fixed parameters (Table 3) ----------
global theta BetaHL BetaR bmax Kappa1 Kappa3 mu1 mu3 mu4 gamma Alpha
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

% ---------- Fitted parameters (Table 5) ----------
global sigma1 BetaA phiF
sigma1 = 0.954;
BetaA  = 1.159;
phiF   = 1.365;

% ---------- Initial condition ----------
Astar = 0.617;
S0H0  = LambdaH / (muH + BetaA * Astar);
S1H0  = (BetaA * Astar / muH) * S0H0;
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

% ---------- Time grid ----------
T    = 0.5;
Nmax = 2000;
tt   = 0 : T/Nmax : T;

% =========================================================
%  BASELINE SIMULATION
% =========================================================
[TT_b, X_b] = run_model(sigma1, BetaA, phiF, theta, BetaHL, ...
                        BetaR, bmax, Kappa1, Kappa3, ...
                        mu1, mu3, mu4, gamma, Alpha);

% =========================================================
%  FIGURE 2 — Baseline outbreak trajectory (4 panels)
% =========================================================
fig = figure('Units','pixels','Position',[100 100 1200 800]);

ax1 = subplot(2,2,1); hold(ax1, 'on');
shade_climate_window(ax1);
plot(TT_b, X_b(:,1), 'b-', 'LineWidth', 1.3);
plot(TT_b, X_b(:,2), 'r-', 'LineWidth', 1.3);
plot(TT_b, X_b(:,3), 'k-', 'LineWidth', 1.3);
grid(ax1, 'on');
xlabel('t (yr)'); ylabel('Individuals');
legend('S_L','E_L','I_L','Location','best');
title('(a) Livestock compartments');

ax2 = subplot(2,2,2); hold(ax2, 'on');
shade_climate_window(ax2);
plot(TT_b, X_b(:,7) + X_b(:,8), 'm-', 'LineWidth', 1.5);
grid(ax2, 'on');
xlabel('t (yr)'); ylabel('I_H (aggregate)');
title('(b) Aggregate human infectives');

ax3 = subplot(2,2,3); hold(ax3, 'on');
shade_climate_window(ax3);
yyaxis(ax3, 'left');
plot(TT_b, X_b(:,9), 'b-', 'LineWidth', 1.5); ylabel(ax3, 'A');
yyaxis(ax3, 'right');
plot(TT_b, X_b(:,10), 'r-', 'LineWidth', 1.5); ylabel(ax3, 'R');
xlabel('t (yr)'); grid(ax3, 'on');
title('(c) Awareness A and rumour R');

ax4 = subplot(2,2,4); hold(ax4, 'on');
shade_climate_window(ax4);
plot(TT_b, X_b(:,13), 'r-', 'LineWidth', 1.3);
plot(TT_b, X_b(:,14), 'k-', 'LineWidth', 1.3);
grid(ax4, 'on');
xlabel('t (yr)'); ylabel('Individuals / eggs');
legend('I_V','M_V','Location','best');
title('(d) Vector compartments');

sgtitle('Baseline outbreak trajectory of the fitted model');
saveas(fig, fullfile(outdir, 'fig2_baseline_traj.png'));
fprintf('Saved fig2_baseline_traj.png\n');

% =========================================================
%  FIGURE 3 — Counterfactual with sigma1 = 0
% =========================================================
[TT_c, X_c] = run_model(0, BetaA, phiF, theta, BetaHL, ...
                        BetaR, bmax, Kappa1, Kappa3, ...
                        mu1, mu3, mu4, gamma, Alpha);

fig = figure('Units','pixels','Position',[100 100 1200 800]);

ax1 = subplot(2,2,1); hold(ax1, 'on');
shade_climate_window(ax1);
plot(TT_b, X_b(:,3), 'b-', 'LineWidth', 1.4);
plot(TT_c, X_c(:,3), 'r--', 'LineWidth', 1.4);
grid(ax1, 'on'); xlabel('t (yr)'); ylabel('I_L');
legend('Baseline (\sigma_1 = 0.954)','No rumour slaughter (\sigma_1 = 0)', ...
       'Location','best');
title('(a) Livestock infectious');

ax2 = subplot(2,2,2); hold(ax2, 'on');
shade_climate_window(ax2);
plot(TT_b, X_b(:,7) + X_b(:,8), 'b-', 'LineWidth', 1.4);
plot(TT_c, X_c(:,7) + X_c(:,8), 'r--', 'LineWidth', 1.4);
grid(ax2, 'on'); xlabel('t (yr)'); ylabel('I_H');
legend('Baseline','No rumour slaughter','Location','best');
title('(b) Aggregate human infectives');

ax3 = subplot(2,2,3); hold(ax3, 'on');
shade_climate_window(ax3);
cum_b = cumtrapz(TT_b, X_b(:,7) + X_b(:,8));
cum_c = cumtrapz(TT_c, X_c(:,7) + X_c(:,8));
plot(TT_b, cum_b, 'b-', 'LineWidth', 1.4);
plot(TT_c, cum_c, 'r--', 'LineWidth', 1.4);
grid(ax3, 'on'); xlabel('t (yr)'); ylabel('Cumulative I_H');
legend('Baseline','No rumour slaughter','Location','best');
title('(c) Cumulative human infectives');

ax4 = subplot(2,2,4); hold(ax4, 'on');
shade_climate_window(ax4);
plot(TT_b, X_b(:,10), 'b-', 'LineWidth', 1.4);
plot(TT_c, X_c(:,10), 'r--', 'LineWidth', 1.4);
grid(ax4, 'on'); xlabel('t (yr)'); ylabel('R');
legend('Baseline','No rumour slaughter','Location','best');
title('(d) Rumour pressure');

sgtitle('Counterfactual without rumour-driven slaughtering');
saveas(fig, fullfile(outdir, 'fig3_counterfactual_sigma1.png'));
fprintf('Saved fig3_counterfactual_sigma1.png\n');

% =========================================================
%  FIGURE 4 — Counterfactual with phiF = 3
% =========================================================
[TT_c2, X_c2] = run_model(sigma1, BetaA, 3, theta, BetaHL, ...
                          BetaR, bmax, Kappa1, Kappa3, ...
                          mu1, mu3, mu4, gamma, Alpha);

fig = figure('Units','pixels','Position',[100 100 1200 800]);

ax1 = subplot(2,2,1); hold(ax1, 'on');
shade_climate_window(ax1);
plot(TT_b, X_b(:,7) + X_b(:,8), 'b-', 'LineWidth', 1.4);
plot(TT_c2, X_c2(:,7) + X_c2(:,8), 'r:', 'LineWidth', 1.4);
grid(ax1, 'on'); xlabel('t (yr)'); ylabel('I_H');
legend('Fitted (\phi_F = 1.365)','Strong (\phi_F = 3)','Location','best');
title('(a) Aggregate human infectives');

ax2 = subplot(2,2,2); hold(ax2, 'on');
shade_climate_window(ax2);
cum_b  = cumtrapz(TT_b,  X_b(:,7)  + X_b(:,8));
cum_c2 = cumtrapz(TT_c2, X_c2(:,7) + X_c2(:,8));
plot(TT_b, cum_b, 'b-', 'LineWidth', 1.4);
plot(TT_c2, cum_c2, 'r:', 'LineWidth', 1.4);
grid(ax2, 'on'); xlabel('t (yr)'); ylabel('Cumulative I_H');
legend('Fitted','Strong','Location','best');
title('(b) Cumulative human infectives');

ax3 = subplot(2,2,3); hold(ax3, 'on');
shade_climate_window(ax3);
plot(TT_b, X_b(:,10), 'b-', 'LineWidth', 1.4);
plot(TT_c2, X_c2(:,10), 'r:', 'LineWidth', 1.4);
grid(ax3, 'on'); xlabel('t (yr)'); ylabel('R');
legend('Fitted','Strong','Location','best');
title('(c) Rumour pressure');

ax4 = subplot(2,2,4); hold(ax4, 'on');
shade_climate_window(ax4);
plot(TT_b, X_b(:,9), 'b-', 'LineWidth', 1.4);
plot(TT_c2, X_c2(:,9), 'r:', 'LineWidth', 1.4);
grid(ax4, 'on'); xlabel('t (yr)'); ylabel('A');
legend('Fitted','Strong','Location','best');
title('(d) Awareness');

sgtitle('Counterfactual with strong fact-checking');
saveas(fig, fullfile(outdir, 'fig4_counterfactual_phi.png'));
fprintf('Saved fig4_counterfactual_phi.png\n');

% =========================================================
%  FIGURE 5 — Counterfactual with alpha = 0
% =========================================================
[TT_c3, X_c3] = run_model(sigma1, BetaA, phiF, theta, BetaHL, ...
                          BetaR, bmax, Kappa1, Kappa3, ...
                          mu1, mu3, mu4, gamma, 0);

fig = figure('Units','pixels','Position',[100 100 1200 600]);

ax1 = subplot(1,2,1); hold(ax1, 'on');
shade_climate_window(ax1);
semilogy(TT_b, max(X_b(:,3), 1e-6), 'b-', 'LineWidth', 1.4);
semilogy(TT_c3, max(X_c3(:,3), 1e-6), 'r--', 'LineWidth', 1.4);
grid(ax1, 'on'); xlabel('t (yr)'); ylabel('I_L (log scale)');
legend('Fitted (\alpha = 0.63)','No climate (\alpha = 0)','Location','best');
title('(a) Infectious livestock');

ax2 = subplot(1,2,2); hold(ax2, 'on');
shade_climate_window(ax2);
semilogy(TT_b, max(X_b(:,7) + X_b(:,8), 1e-6), 'b-', 'LineWidth', 1.4);
semilogy(TT_c3, max(X_c3(:,7) + X_c3(:,8), 1e-6), 'r--', 'LineWidth', 1.4);
grid(ax2, 'on'); xlabel('t (yr)'); ylabel('I_H (log scale)');
legend('Fitted','No climate','Location','best');
title('(b) Aggregate human infectives');

sgtitle('Counterfactual without climate forcing');
saveas(fig, fullfile(outdir, 'fig5_counterfactual_alpha.png'));
fprintf('Saved fig5_counterfactual_alpha.png\n');

% =========================================================
%  FIGURE 6 — Numerical global stability of the epizootic core
% =========================================================
r_grid = linspace(0.5, 2.0, 30);
IL_dfe  = zeros(size(r_grid));
IL_ende = zeros(size(r_grid));
for k = 1:numel(r_grid)
    r = r_grid(k);
    bLV = BetaLV * r;
    bVL = BetaVL / r;
    IL_dfe(k)  = frozen_core_steady_state(bLV, bVL, 'dfe');
    IL_ende(k) = frozen_core_steady_state(bLV, bVL, 'endemic');
end

fig = figure('Units','pixels','Position',[100 100 900 500]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
plot(ax, r_grid, IL_dfe,  'o-', 'LineWidth', 1.3, 'MarkerSize', 6, ...
     'DisplayName', 'DFE start');
plot(ax, r_grid, IL_ende, 's-', 'LineWidth', 1.3, 'MarkerSize', 6, ...
     'DisplayName', 'Endemic start');
xlabel(ax, 'Asymmetry ratio r = \beta_{LV} N_L^*/(\beta_{VL} N_V^*)');
ylabel(ax, 'I_L^* at steady state (head)');
title(ax, 'Numerical global stability check for the epizootic core');
legend(ax, 'Location', 'best');
text(ax, 1.55, max(IL_dfe)*0.6, ...
     sprintf('max residual = %.1e head', max(abs(IL_dfe - IL_ende))));
saveas(fig, fullfile(outdir, 'fig6_global_stability.png'));
fprintf('Saved fig6_global_stability.png\n');

% =========================================================
%  FIGURE 7 — Local elasticity bar chart
% =========================================================
params = {'\sigma_1', '\beta_A', '\phi_F'};
E_CH = [-0.014, -0.171, -0.006];
E_RA = [-0.008,  0,      0    ];
E_RH = [ 0.046, -0.082, -0.003];

fig = figure('Units','pixels','Position',[100 100 900 500]);
bar([E_CH; E_RA; E_RH]', 'grouped');
set(gca, 'XTickLabel', params);
ylabel('Elasticity');
legend('C_H','R_A','R_H','Location','best');
grid on;
title('Local elasticities with respect to the estimated parameters');
saveas(fig, fullfile(outdir, 'fig7_elasticity.png'));
fprintf('Saved fig7_elasticity.png\n');

fprintf('\nAll figures saved to: %s\n', outdir);

end

% =============================================================
%  Helper functions
% =============================================================
function [TT, X] = run_model(s1, bA, pF, th, bHL, bR, bm, k1, k3, ...
                             m1, m3, m4, gm, al)
    global tt Z0
    global sigma1 BetaA phiF theta BetaHL BetaR bmax Kappa1 Kappa3 ...
           mu1 mu3 mu4 gamma Alpha
    sigma1 = s1; BetaA = bA; phiF = pF; theta = th; BetaHL = bHL;
    BetaR = bR; bmax = bm; Kappa1 = k1; Kappa3 = k3; mu1 = m1;
    mu3 = m3; mu4 = m4; gamma = gm; Alpha = al;
    opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9, 'NonNegative', 1:14);
    [TT, X] = ode15s(@hsecond, tt, Z0, opts);
end

function shade_climate_window(ax)
    % Draw the shading as the first child of the axes. Call this before
    % any plot() call on the same axes so the shading stays behind the
    % data. Sets the current axes to 'ax' as a side effect.
    axes(ax);                      %#ok<LAXES>
    yl = ylim(ax);
    patch(ax, [0 1/3 1/3 0], [yl(1) yl(1) yl(2) yl(2)], ...
          [0.9 0.9 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold(ax, 'on');
end

function IL_ss = frozen_core_steady_state(bLV, bVL, init_type)
    global LambdaL muL gammaL epslonL LambdaV muV epslonV RoV sigma0
    if strcmp(init_type, 'dfe')
        x0 = [LambdaL/muL; 0; 0; LambdaV/muV; 0; 0; 0];
    else
        x0 = [LambdaL/muL*0.5; 1; 50; LambdaV/muV*0.5; 1; 1; 10];
    end
    tspan = [0 200];
    o = odeset('RelTol', 1e-7, 'AbsTol', 1e-9, 'NonNegative', 1:7);
    [~, X] = ode15s(@(t,x) core_ode(t, x, bLV, bVL), tspan, x0, o);
    IL_ss = X(end, 3);
end

function dx = core_ode(~, x, bLV, bVL)
    global LambdaL muL gammaL epslonL LambdaV muV epslonV RoV sigma0
    SL = x(1); EL = x(2); IL = x(3);
    SV = x(4); EV = x(5); IV = x(6); MV = x(7);
    NL = SL + EL + IL;
    NV = SV + EV + IV;
    bL = 182 * NL / (NL + 1000);
    lamL = bL * bLV * IV / max(NV, 1e-6);
    lamV = bL * bVL * IL / max(NL, 1e-6);
    dx = zeros(7, 1);
    dx(1) = LambdaL - lamL*SL - muL*SL;
    dx(2) = lamL*SL - epslonL*EL - muL*EL;
    dx(3) = epslonL*EL - (gammaL + muL)*IL - sigma0*IL;
    dx(4) = LambdaV + 0.2*MV - lamV*SV - muV*SV;
    dx(5) = lamV*SV - epslonV*EV - muV*EV;
    dx(6) = epslonV*EV - muV*IV;
    dx(7) = RoV*IV - 0.2*MV - 0.5*MV;
end