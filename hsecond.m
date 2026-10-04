function Z = hsecond(t, x)
% RHS of the 14-dimensional coupled infodemic-epizootic RVF model.
%
% State ordering:
%   x(1) SL    x(2) EL    x(3) IL    x(4) RL
%   x(5) S0H   x(6) S1H   x(7) I0H   x(8) I1H
%   x(9) A     x(10) R    x(11) SV   x(12) EV
%   x(13) IV   x(14) MV

% Literature parameters
global LambdaH LambdaL LambdaV muH muL muV muM gammaH deltaH gammaL epslonL
global epslonV RoV sigma0 Phi0 Phi1 BetaHV BetaLV BetaVL Kappa2 Kappa4 mu2

% Parameters — 5 estimated, 8 fixed
global sigma1 BetaA theta BetaHL phiF              % estimated
global BetaR bmax Kappa1 Kappa3 mu1 mu3 mu4 gamma Alpha  % fixed

% ---------- Climate forcing C(t), Eq. (7) ----------
% t in years, t = 0 corresponds to 1 Sep 2006.
if     t < 0
    Ct = 0;
elseif t < 1/12
    Ct = 0.53;      % Sep 2006
elseif t < 2/12
    Ct = 0.73;      % Oct 2006
elseif t < 3/12
    Ct = 1.03;      % Nov 2006
elseif t < 4/12
    Ct = 0.70;      % Dec 2006
else
    Ct = 0;
end

sigmaLR = sigma0 + sigma1 * x(10);
phit    = Phi0   + Phi1 * Alpha * Ct;

NL = x(1) + x(2) + x(3) + x(4);
NH = x(5) + x(6) + x(7) + x(8);
NV = x(11) + x(12) + x(13);

if NL < 1e-6, NL = 1e-6; end
if NH < 1e-6, NH = 1e-6; end
if NV < 1e-6, NV = 1e-6; end

bLt = bmax * NL / (NL + NH);
bHt = bmax * NH / (NL + NH);

lambda0H = BetaHL * sigmaLR * x(3) / NL + bHt * BetaHV * x(13) / NV;
lambda1H = (1 - theta) * lambda0H;
lambdaL  = bLt * BetaLV * x(13) / NV;
lambdaV  = bLt * BetaVL * x(3)  / NL;

IH = x(7) + x(8);

Z1  = LambdaL - lambdaL * x(1) - muL * x(1);
Z2  = lambdaL * x(1) - epslonL * x(2) - muL * x(2);
Z3  = epslonL * x(2) - (gammaL + muL) * x(3) - sigmaLR * x(3);
Z4  = gammaL * x(3) - muL * x(4);
Z5  = LambdaH - lambda0H * x(5) - BetaA * x(9) * x(5) ...
      + BetaR * x(10) * x(6) - muH * x(5);
Z6  = BetaA * x(9) * x(5) - lambda1H * x(6) ...
      - BetaR * x(10) * x(6) - muH * x(6);
Z7  = lambda0H * x(5) - BetaA * x(9) * x(7) ...
      + BetaR * x(10) * x(8) - (gammaH + muH + deltaH) * x(7);
Z8  = lambda1H * x(6) + BetaA * x(9) * x(7) ...
      - BetaR * x(10) * x(8) - (gammaH + muH + deltaH) * x(8);
Z9  = Kappa1 * Alpha * Ct * (1 - x(9)) + Kappa2 * x(9) * (1 - x(9)) ...
      - Kappa3 * x(10) * x(9) - Kappa4 * x(9);
Z10 = mu1 * IH - mu2 * x(10) + mu3 * (x(10) / (1 + gamma * x(10))) ...
      - mu4 * x(10) * x(9) - phiF * x(10);
Z11 = LambdaV + phit * x(14) - lambdaV * x(11) - muV * x(11);
Z12 = lambdaV * x(11) - epslonV * x(12) - muV * x(12);
Z13 = epslonV * x(12) - muV * x(13);
Z14 = RoV * x(13) - phit * x(14) - muM * x(14);

Z = [Z1; Z2; Z3; Z4; Z5; Z6; Z7; Z8; Z9; Z10; Z11; Z12; Z13; Z14];
end