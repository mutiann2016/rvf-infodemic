function [Time, Zol] = solution(q)
% Integrate the model with the 3-vector of estimated parameters.
% q = [sigma1; BetaA; phiF]

global tt Z0

% Estimated parameters
global sigma1 BetaA phiF
sigma1 = q(1);
BetaA  = q(2);
phiF   = q(3);

% Fixed parameters — manuscript Table 3 and Table 5 values
global theta BetaHL BetaR bmax Kappa1 Kappa3 mu1 mu3 mu4 gamma Alpha
theta  = 0.47;         % awareness efficacy (manuscript)
BetaHL = 0.05;         % livestock-to-human transmission (manuscript, fixed on identifiability grounds)
BetaR  = 0.33;
bmax   = 182;
Kappa1 = 0.29;
Kappa3 = 0.21;
mu1    = 1.7e-4;
mu3    = 0.38;
mu4    = 0.26;
gamma  = 0.031;
Alpha  = 0.63;

opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9, 'NonNegative', 1:14);
try
    [TT, Z] = ode15s(@hsecond, tt, Z0, opts);
catch
    TT = []; Z = [];
end
Time = TT; Zol = Z;
end