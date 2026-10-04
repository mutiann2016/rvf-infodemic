function [qopt, FVAL] = main
% RVF 3-parameter fit (real data only)
% Estimates: sigma1, BetaA, phiF
% Fixed:     theta = 0.47, BetaHL = 0.05, BetaR = 0.33, bmax = 182,
%            Kappa1 = 0.29, Kappa3 = 0.21, mu1 = 1.7e-4, mu3 = 0.38,
%            mu4 = 0.26, gamma = 0.031, Alpha = 0.63

clc; close all;
set(0,'defaultaxesfontsize',13,'defaultaxeslinewidth',.7, ...
     'defaultlinelinewidth',.9,'defaultpatchlinewidth',.7);

diary_file = fullfile(pwd, 'RVF_3param_run.log');
diary(diary_file);
fprintf('\n=== RVF 3-parameter fit: %s ===\n', datestr(now));
fprintf('Working folder: %s\n', pwd);

% ---------- Literature parameters ----------
Pardef

global Z0 tt
global LambdaH LambdaL LambdaV muH muL muV muM
global gammaH deltaH gammaL epslonL epslonV RoV
global sigma0 Phi0 Phi1 BetaHV BetaLV BetaVL Kappa2 Kappa4 mu2
global sigma1 BetaA phiF
global theta BetaHL BetaR bmax Kappa1 Kappa3 mu1 mu3 mu4 gamma Alpha

% ---------- Initial condition: DFE + one exposed animal ----------
BetaA_init = 0.41;
Astar      = 0.421;
S0H0 = LambdaH / (muH + BetaA_init * Astar);
S1H0 = (BetaA_init * Astar / muH) * S0H0;
I0H0 = 0;
I1H0 = 0;

SL0  = LambdaL / muL;
EL0  = 1;
IL0  = 0;
RL0  = 0;
A0   = Astar;
R0v  = 0;
SV0  = LambdaV / muV;
EV0  = 0;
IV0  = 0;
MV0  = 0;

Z0 = [SL0; EL0; IL0; RL0; S0H0; S1H0; I0H0; I1H0; A0; R0v; SV0; EV0; IV0; MV0];

% ---------- Time grid ----------
T    = 1.0;
Nmax = 500;
tt   = 0 : T/Nmax : T;

% ---------- Estimated-parameter vector: q = [sigma1; BetaA; phiF] ----------
q0 = [0.84; 0.41; 0.12];
LB = [0.001; 0.01; 0.0 ];
UB = [5.0;   2.0;  3.0 ];

q0 = max(q0, LB);
q0 = min(q0, UB);

% ---------- fmincon options ----------
OPTIONS = optimset('fmincon');
OPTIONS = optimset(OPTIONS, 'Display','iter', ...
    'MaxFunEvals', 3000, 'MaxIter', 300, 'TolFun', 1e-5, 'TolX', 1e-5);

% ---------- SMOKE TEST ----------
fprintf('\n--- SMOKE TEST ---\n');
[TT0, X0] = solution(q0);
if isempty(TT0) || any(~isfinite(X0(:)))
    fprintf('SMOKE TEST FAILED.\n');
    qopt = q0; FVAL = 1e12; diary off; return;
end
fprintf('SMOKE TEST PASSED. I_L(end) = %g, I_H(end) = %g\n', ...
    X0(end,3), X0(end,7) + X0(end,8));

% ---------- fmincon ----------
fprintf('\n--- fmincon ---\n');
[qopt, FVAL, EXITFLAG] = fmincon(@objective, q0, [],[],[],[], ...
                                 LB, UB, [], OPTIONS);

% ---------- Estimated parameter names ----------
names = {'sigma1','BetaA','phiF'};

% ---------- Save results to file ----------
outfile = fullfile(pwd, 'RVF_fit_output.txt');
fid = fopen(outfile, 'w');
if fid < 0
    fprintf('ERROR: cannot open output file.\n');
else
    fprintf(fid, 'RVF 3-parameter fit (real data only)\n');
    fprintf(fid, 'Fit run: %s\n\n', datestr(now));
    fprintf(fid, 'Estimated parameters:\n');
    for k = 1:numel(names)
        fprintf(fid, '  %-8s = %.8g\n', names{k}, qopt(k));
    end
    fprintf(fid, '\nFVAL      = %.8g\n', FVAL);
    fprintf(fid, 'EXITFLAG  = %d\n',   EXITFLAG);
    fprintf(fid, '\nFixed parameters (Table 3, Table 5 of the manuscript):\n');
    fprintf(fid, '  theta    = %.4g\n', theta);
    fprintf(fid, '  BetaHL   = %.4g\n', BetaHL);
    fprintf(fid, '  BetaR    = %.4g\n', BetaR);
    fprintf(fid, '  bmax     = %.4g\n', bmax);
    fprintf(fid, '  Kappa1   = %.4g\n', Kappa1);
    fprintf(fid, '  Kappa3   = %.4g\n', Kappa3);
    fprintf(fid, '  mu1      = %.4g\n', mu1);
    fprintf(fid, '  mu3      = %.4g\n', mu3);
    fprintf(fid, '  mu4      = %.4g\n', mu4);
    fprintf(fid, '  gamma    = %.4g\n', gamma);
    fprintf(fid, '  Alpha    = %.4g\n', Alpha);
    fclose(fid);
    fprintf('\n>>> SAVED: %s\n', outfile);
end

% =========================================================
%  Derived quantities
% =========================================================
sigma1_ = qopt(1);   % fitted sigma1
BetaA_  = qopt(2);   % fitted BetaA
phiF_   = qopt(3);   % fitted phiF

Cbar = (0.53 + 0.73 + 1.03 + 0.70) / 4;

% Awareness-only equilibrium A*, Eq. (20) of manuscript
a_ = Kappa2;
b_ = -(Kappa2 - Kappa1 * Alpha * Cbar - Kappa4);
c_ = -Kappa1 * Alpha * Cbar;
Astar_fit = (-b_ + sqrt(b_^2 - 4*a_*c_)) / (2*a_);

% R_A, Eq. (22) of manuscript
RA_fit = sqrt( (bmax^2 * BetaLV * BetaVL * epslonL * epslonV) / ...
               ((gammaL + muL + sigma0) * (epslonL + muL) * ...
                (epslonV + muV) * muV) );

% R_H, Eq. (21) of manuscript
S0Hstar = LambdaH / (muH + BetaA_ * Astar_fit);
S1Hstar = BetaA_ * Astar_fit * S0Hstar / muH;
RH_fit  = (BetaHL * sigma0 * (S0Hstar + (1 - theta) * S1Hstar) / ...
           (LambdaL/muL)) / (gammaL + muL + sigma0) * ...
          (epslonL / (epslonL + muL));

% Fact-checking threshold, Eq. (41) of manuscript
phiFc = mu3 - mu2 - mu4 * Astar_fit;

% ---------- Report ----------
fprintf('\n========== ESTIMATED PARAMETERS ==========\n');
for k = 1:numel(names)
    fprintf('  %-8s = %g\n', names{k}, qopt(k));
end
fprintf('  FVAL      = %g\n', FVAL);
fprintf('  EXITFLAG  = %d\n', EXITFLAG);

fprintf('\n========== DERIVED QUANTITIES ==========\n');
fprintf('  A*        = %.6f\n', Astar_fit);
fprintf('  R_A       = %.6f\n', RA_fit);
fprintf('  R_H       = %.6f\n', RH_fit);
fprintf('  phi_F^c   = %.6f\n', phiFc);
fprintf('========================================\n');

% ---------- Append derived quantities to file ----------
fid = fopen(outfile, 'a');
if fid > 0
    fprintf(fid, '\nDerived quantities\n');
    fprintf(fid, '  A*        = %.6f\n', Astar_fit);
    fprintf(fid, '  R_A       = %.6f\n', RA_fit);
    fprintf(fid, '  R_H       = %.6f\n', RH_fit);
    fprintf(fid, '  phi_F^c   = %.6f\n', phiFc);
    fclose(fid);
    fprintf('>>> Appended derived quantities to %s\n', outfile);
end

% =========================================================
%  Plot model fit vs data
% =========================================================
[TT, X] = solution(qopt);

CHobs_months = [20; 340; 244; 75; 25];
CHtimes      = [2/12; 3/12; 4/12; 5/12; 6/12];

CHmodel = zeros(5,1);
for m = 1:5
    tm = CHtimes(m);
    lo = tm - 1/24;  hi = tm + 1/24;
    I1 = find(TT >= lo, 1, 'first');
    I2 = find(TT >= hi, 1, 'first');
    if isempty(I1) || isempty(I2) || I2 <= I1, continue; end
    NL_loc = X(I1:I2,1) + X(I1:I2,2) + X(I1:I2,3) + X(I1:I2,4);
    NL_loc = max(NL_loc, 1e-6);
    sigma_loc = sigma0 + sigma1 * X(I1:I2,10);
    lam0_loc  = BetaHL * sigma_loc .* X(I1:I2,3) ./ NL_loc;
    inc = lam0_loc .* X(I1:I2,5) + (1 - theta) * lam0_loc .* X(I1:I2,6);
    CHmodel(m) = trapz(TT(I1:I2), inc);
end

figure(1);
bar([CHobs_months, CHmodel]); grid on;
set(gca, 'XTickLabel', {'Nov06','Dec06','Jan07','Feb07','Mar07'});
ylabel('Confirmed human cases');
legend('Observed (WHO)','Model','Location','best');
title('Monthly human cases: model vs. WHO data');

figure(2);
plot(TT, X(:,7) + X(:,8), 'b-', 'LineWidth', 1.4); grid on;
xlabel('t (yr, t=0 is Sep 2006)'); ylabel('I_H (aggregate human infectives)');
title('Aggregate human infective trajectory (model)');

figure(3);
plot(TT, X(:,3), 'r-', 'LineWidth', 1.4); grid on;
xlabel('t (yr, t=0 is Sep 2006)'); ylabel('I_L (livestock infectious)');
title('Livestock infectious trajectory (model)');

diary off;
end