function Pardef
% Literature parameters (Table 3 of the manuscript).

global LambdaH LambdaL LambdaV muH muL muV muM gammaH deltaH gammaL epslonL
global epslonV RoV sigma0 Phi0 Phi1 BetaHV BetaLV BetaVL Kappa2 Kappa4 mu2

LambdaH = 1000;
LambdaL = 1000;
LambdaV = 10000;
muH     = 0.014;
muL     = 0.1;
muV     = 12;
muM     = 0.5;
gammaH  = 26;
deltaH  = 0.30 * gammaH;         % = 7.8 yr^-1
gammaL  = 50;
epslonL = 36;
epslonV = 18;
RoV     = 0.1;
sigma0  = 0.5;
Phi0    = 0.2;
Phi1    = 0.8;
BetaHV  = 0.01;
BetaLV  = 0.3;
BetaVL  = 0.3;
Kappa2  = 0.3;
Kappa4  = 0.2;
mu2     = 0.3;
end