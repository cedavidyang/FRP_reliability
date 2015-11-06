function shearTotalKN = shear_total_ACInew_MC(iDesignCase, hSmp, betaSmp, fFrpSmp, fcSmp, sqrtFcSmp, fsSmp)
% total shear resistance based on ACI 440.2R-08

load tmpdata.mat

FRP_STRN_LMT_ABS = 1e8; % 0.004;
FRP_STRN_LMT_REL = 0.8; % 0.75;

% read variables        
nCase = length(hSmp);
% geometrical properties
dBeamMM = hSmp - AS_DESIGN_ARRAY_MM(iDesignCase);
bBeamMM = hSmp/H2B_DESIGN_ARRAY(iDesignCase);
dFrpTopMM = DFRP_TOP_DESIGN_ARRAY_MM(iDesignCase)*ones(nCase, 1);
% Concrete properties
fcMPA = fcSmp;
sqrt_fc = sqrt(fcSmp); sqrt_fc(fcSmp<0) = 0;
% steel properties
fsMPA = fsSmp;
sdMM = SD_DESIGN_ARRAY_MM(iDesignCase)*ones(nCase, 1);
ssMM = SS_DESIGN_ARRAY_MM(iDesignCase)*ones(nCase, 1);
% FRP properties
frpForm = FRP_FORM_DESIGN_ARRAY(iDesignCase)*ones(nCase, 1);
fFrpMPA = fFrpSmp;
EFrpMPA = E_FRP_DESIGN_ARRAY_MPA(iDesignCase)*ones(nCase, 1);
betaFrpDEG = betaSmp;
tFrpMM = T_FRP_DESIGN_ARRAY_MM(iDesignCase)*ones(nCase, 1);
widthFrpMM = W_FRP_DESIGN_ARRAY_MM(iDesignCase)*ones(nCase, 1);
sFrpMM = S_FRP_DESIGN_ARRAY_MM(iDesignCase)*ones(nCase, 1);
% limit values
fsLimit = 420 * FS_BIAS;
sqrtFcLimit = 8.3 * FCT_BIAS;

%% ACI 440.2R: FRP contribution (Eq. 11-3)

dfvMM = dBeamMM - dFrpTopMM;
ep_fu = fFrpMPA ./ EFrpMPA;
betaRAD = betaFrpDEG /180*pi;

isSide = (frpForm == 1);
isU = (frpForm == 2);
isUorSide = (frpForm == 1) | (frpForm == 2);
isW = (frpForm == 3);

%% Side-bonding or U-jacketing
leMM = zeros(nCase, 1);
lmaxMM = zeros(nCase, 1);
k1 = ones(nCase, 1);
k2 = ones(nCase, 1);
k3 = ones(nCase, 1);
kv = FRP_STRN_LMT_REL*ones(nCase, 1);
ep_fe = zeros(nCase, 1);
lmbd = ones(nCase,1);
% active bond length Eq. (11-8)
leMM(isUorSide) = sqrt(EFrpMPA(isUorSide).*tFrpMM(isUorSide)./sqrt(fcMPA(isUorSide)));
lmaxMM(isU) = dfvMM(isU) ./ sin(betaRAD(isU));
lmaxMM(isSide) = dfvMM(isSide) ./ (2*sin(betaRAD(isSide)));
lmbd(isUorSide) = lmaxMM(isUorSide)./leMM(isUorSide);
% parameter new k1             
k1(isUorSide) = sqrt(fcMPA(isUorSide)/10);             
% parameter new k2
indx1 = isUorSide&(lmbd<=1);
k2(indx1) = 0.7909*lmbd(indx1) - 0.02818*lmbd(indx1).^2 - 0.1261*lmbd(indx1).^3;
indx2 = isUorSide&(lmbd>1);
k2(indx2) = 1-(pi-2)./(pi*lmbd(indx2));
% new parameter k3
k3 = sqrt( (2-widthFrpMM./(sFrpMM.*sin(betaRAD)))./(1+widthFrpMM./(sFrpMM.*sin(betaRAD))) );
k3 = ones(nCase, 1);
% parameter kv Eq. (11-7)
kv(isUorSide) = k1(isUorSide) .* k2(isUorSide) .* k3(isUorSide) .* leMM(isUorSide) ./...
    (EFrpMPA(isUorSide).*tFrpMM(isUorSide).*ep_fu(isUorSide));
kv(isUorSide & (kv<0) ) = 0;
kv(isUorSide & (kv>FRP_STRN_LMT_REL) ) = FRP_STRN_LMT_REL;

ep_fe(isUorSide) = kv(isUorSide) .* ep_fu(isUorSide);

covEp = std(ep_fu) / mean(ep_fu);
ep_fe( ep_fe>FRP_STRN_LMT_ABS ) = normrnd(FRP_STRN_LMT_ABS, FRP_STRN_LMT_ABS*covEp, sum(ep_fe>FRP_STRN_LMT_ABS), 1);

%% Wrapping

% covEpFe = F_FRP_COV;
% ep_fe(isW) = normrnd(FRP_STRN_LMT_ABS, FRP_STRN_LMT_ABS*covEpFe, sum(isW), 1);
% isWandOver = (isW) & (ep_fe > FRP_STRN_LMT_REL*ep_fu);
% ep_fe( isWandOver ) = FRP_STRN_LMT_REL*ep_fu(isWandOver);
% ep_fe( ep_fe<0 ) = 0;

% uncomment the folllwing block if you want to used HK model to predict W
% it's not necceaary since FRP_STRN_LMT_ABS will control nearly all W tests.
MAX_FRP = FRP_STRN_LMT_REL; TOP_DUMMY_HEIGHT=0.1;
zt = dFrpTopMM;
zb = (dBeamMM)-TOP_DUMMY_HEIGHT*dBeamMM;
distributeFrp = zeros(nCase, 1); sFrpDesignMax = zeros(nCase,1);
zeta = zt./zb;
distributeFrp(isW) = (1+zeta(isW))/2;
sFrpDesignMax(isW) = MAX_FRP * fFrpMPA(isW);
fFrpDesignEff = sFrpDesignMax .* distributeFrp;
ep_fe(isW) = fFrpDesignEff(isW) ./ EFrpMPA(isW);
isWandOverRel = (isW) & (ep_fe > FRP_STRN_LMT_REL*ep_fu);
ep_fe( isWandOverRel ) = FRP_STRN_LMT_REL*ep_fu(isWandOverRel);
iwWandOverAbs = (isW) & (ep_fe > FRP_STRN_LMT_ABS);
ep_fe( iwWandOverAbs ) = normrnd(FRP_STRN_LMT_ABS, FRP_STRN_LMT_ABS*F_FRP_COV, sum(iwWandOverAbs), 1);

%% FRP contribution to shear resistance Eq. (11-3)

f_feMPA = ep_fe .* EFrpMPA;
shearFrpKN = 2*tFrpMM .* widthFrpMM .* f_feMPA .*...
        (sin(betaRAD)+cos(betaRAD)) .* dfvMM ./ sFrpMM /1e3;
shearFrpKN( shearFrpKN<0 ) = 0;

%% ACI 318-11: Steel Contribution
nLarger = sum( fsMPA>fsLimit );
covLarger = FS_COV;
fsMPA( fsMPA>fsLimit ) = normrnd(fsLimit, fsLimit*covLarger, nLarger,1);
% fsMPA( fsMPA>fsLimit ) = fsLimit;

areaVMM2 = 2*pi/4*sdMM.^2;
shearSteelKN = areaVMM2 .* fsMPA .* dBeamMM ./ ssMM /1e3;
shearSteelKN( ssMM == 0 ) = 0;
shearSteelKN( shearSteelKN<0 ) = 0;

%% contribution of concrete

areaV_Min = 0.062*sqrt(fcMPA) .* bBeamMM .* ssMM ./ fsMPA;
tmp = (0.35*bBeamMM .* sdMM) ./ fsMPA;
areaV_Min( areaV_Min < tmp ) = tmp(areaV_Min < tmp); 

isAdjust = (areaVMM2<areaV_Min) & (sqrt_fc>sqrtFcLimit);
nAdjust = sum( isAdjust );
covAdjust = std(sqrt_fc)/mean(sqrt_mean);
sqrt_fc( isAdjust ) = normrnd( sqrtFcLimit, sqrtFcLimit*covAdjust, nAdjust, 1);
% sqrt_fc( isAdjust ) = sqrtFcLimit;

shearConcreteKN = 0.17*sqrt_fc .* bBeamMM .* dBeamMM / 1e3;
shearConcreteKN(shearConcreteKN<0) = 0;

%% FRP and steel interaction
isWithSteel = areaVMM2~=0;
mu = zeros(nCase,1);
mu(isSide&isWithSteel) = shearSteelKN(isSide&isWithSteel) ./ shearFrpKN(isSide&isWithSteel) ;
mu(isSide&(shearFrpKN == 0) ) = 0;
kfrp = -0.2*mu+1; kfrp(kfrp<0) = 0;
shearFrpKN(isSide) = kfrp(isSide) .* shearFrpKN(isSide);
shearFrpKN( isSide&shearFrpKN<0 ) = 0;

%% over reinforcement check

shearReinforceKN = shearSteelKN + shearFrpKN;

isOverReinforce = (shearReinforceKN > 0.66*sqrt_fc.*bBeamMM .* dBeamMM / 1e3);
shearReinforceKN( isOverReinforce ) = 0.66*sqrt_fc(isOverReinforce).*...
                       bBeamMM(isOverReinforce).*dBeamMM(isOverReinforce)/1e3;
                   
shearTotalKN = shearConcreteKN + shearReinforceKN;

return
end