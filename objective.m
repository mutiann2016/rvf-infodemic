function [J] = objective(q)
% Reduced objective: 3 estimated parameters, 7 real scalar targets.
% q = [sigma1; BetaA; phiF]

global LambdaH LambdaL LambdaV muH muL muV muM gammaH deltaH gammaL epslonL
global epslonV RoV sigma0 Phi0 Phi1 BetaHV BetaLV BetaVL Kappa2 Kappa4 mu2

global sigma1 BetaA phiF
global theta BetaHL BetaR bmax Kappa1 Kappa3 mu1 mu3 mu4 gamma Alpha

sigma1 = q(1);
BetaA  = q(2);
phiF   = q(3);

global tt Z0

% ---------- Real published observations ----------
CHobs_months  = [20; 340; 244; 75; 25];              % Nov 06 - Mar 07
CHtimes       = [2/12; 3/12; 4/12; 5/12; 6/12];
CHsig         = [20; 40; 35; 15; 8];
CH_total_obs  = 684;   CH_total_sig = 50;
DH_total_obs  = 155;   DH_total_sig = 10;
A_L_obs       = 4.1e4; A_L_sig      = 6.0e3;

% ---------- Solve model ----------
[TT, X] = solution(q);
if isempty(TT) || any(~isfinite(X(:)))
    J = 1e10;
    return;
end

J = 0;

% Monthly confirmed cases
for m = 1:5
    tm = CHtimes(m);
    lo = tm - 1/24; hi = tm + 1/24;
    I1 = find(TT >= lo, 1, 'first');
    I2 = find(TT >= hi, 1, 'first');
    if isempty(I1) || isempty(I2) || I2 <= I1, continue; end
    NL_loc = X(I1:I2,1) + X(I1:I2,2) + X(I1:I2,3) + X(I1:I2,4);
    NL_loc = max(NL_loc, 1e-6);
    sigma_loc = sigma0 + sigma1 * X(I1:I2,10);
    lam0_loc  = BetaHL * sigma_loc .* X(I1:I2,3) ./ NL_loc;
    inc = lam0_loc .* X(I1:I2,5) + (1 - theta) * lam0_loc .* X(I1:I2,6);
    CHmodel = trapz(TT(I1:I2), inc);
    J = J + (CHmodel - CHobs_months(m))^2 / CHsig(m)^2;
end

% Cumulative livestock abortions
A_L_model = gammaL * trapz(TT, X(:,3));
J = J + (A_L_model - A_L_obs)^2 / A_L_sig^2;

% Cumulative human cases and deaths
IHinc = (sigma0 + sigma1 * X(:,10)) .* X(:,3) ./ max(X(:,1)+X(:,2)+X(:,3)+X(:,4), 1e-6);
inc_total = BetaHL * IHinc .* X(:,5) + BetaHL * (1 - theta) * IHinc .* X(:,6);
CH_total_model = trapz(TT, inc_total);
J = J + (CH_total_model - CH_total_obs)^2 / CH_total_sig^2;

% Deaths (approx 23% of cases)
DH_total_model = 0.23 * CH_total_model;
J = J + (DH_total_model - DH_total_obs)^2 / DH_total_sig^2;

% Tikhonov regularisation
L = diff(eye(3), 2);
mu_reg = 3.2e-3;
J = J + mu_reg * norm(L * q)^2;

if ~isfinite(J), J = 1e10; end
end