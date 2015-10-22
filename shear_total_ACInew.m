function [shearTotalKN, isOverReinforce,...
         sidewarning, shearReinforceKN] = shear_total_ACInew(FLAG, factorFrp)
% function [shearTotalKN, isOverReinforce,...
%          ro, shearReinforceKN] = shear_total_ACI(FLAG, factorFrp)
% total shear resistance based on ACI 440.2R-08
PHI = 0.85; % resistance reduce factor (ACI 318-11 9.3.2.3)
FRP_STRN_LMT_ABS = 1e8;
FRP_STRN_LMT_REL = 0.8;

load tmpdata.mat

% read variables
switch FLAG
    case 'MODEL_ERROR'
        nCase = length(FC_TEST_ARRAY_MPA);
        % geometrical properties
        dBeamMM = D_TEST_ARRAY_MM;
        bBeamMM = B_TEST_ARRAY_MM;
        dFrpTopMM = DFRP_TOP_TEST_ARRAY_MM;
        % Concrete properties
        fcMPA = FC_TEST_ARRAY_MPA;
        sqrt_fc = sqrt( FC_TEST_ARRAY_MPA );
        % steel properties
        fsMPA = FS_TEST_ARRAY_MPA;
        sdMM = SD_TEST_ARRAY_MM;
        ssMM = SS_TEST_ARRAY_MM;
        % FRP properties
        frpForm = FRP_FORM_TEST_ARRAY;
        fFrpMPA = F_FRP_TEST_ARRAY_MPA;
        EFrpMPA = E_FRP_TEST_ARRAY_MPA;
        betaFrpDEG = BETA_TEST_ARRAY_DEG;
        tFrpMM = T_FRP_TEST_ARRAY_MM;
        widthFrpMM = W_FRP_TEST_ARRAY_MM;
        sFrpMM = S_FRP_TEST_ARRAY_MM;   
        % limit values
        fsLimit = 420 * FS_BIAS;
        sqrtFcLimit = 8.3 * FCT_BIAS;
    case 'DESIGN_VALUE'
        nCase = length(FC_DESIGN_ARRAY_MPA);
        % geometrical properties
        dBeamMM = D_DESIGN_ARRAY_MM;
        bBeamMM = B_DESIGN_ARRAY_MM;
        dFrpTopMM = DFRP_TOP_DESIGN_ARRAY_MM;
        % Concrete properties
        fcMPA = FC_DESIGN_ARRAY_MPA;
        sqrt_fc = sqrt( FC_DESIGN_ARRAY_MPA );
        % steel properties
        fsMPA = FS_DESIGN_ARRAY_MPA;
        sdMM = SD_DESIGN_ARRAY_MM;
        ssMM = SS_DESIGN_ARRAY_MM;
        % FRP properties
        frpForm = FRP_FORM_DESIGN_ARRAY;
        fFrpMPA = F_FRP_DESIGN_ARRAY_MPA;
        EFrpMPA = E_FRP_DESIGN_ARRAY_MPA;
        betaFrpDEG = BETA_DESIGN_ARRAY_DEG;
        tFrpMM = T_FRP_DESIGN_ARRAY_MM;
        widthFrpMM = W_FRP_DESIGN_ARRAY_MM;
        sFrpMM = S_FRP_DESIGN_ARRAY_MM;
        % limit values
        fsLimit = 420;
        sqrtFcLimit = 8.3;        
    case 'MONTE_CARLO'
    otherwise
end

% ACI 440.2R: FRP contribution (Eq. 11-3)

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
% leMM(isUorSide) = 23300 ./ (tFrpMM(isUorSide) .* EFrpMPA(isUorSide)).^0.58;
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
% parameter kv Eq. (11-7)
kv(isUorSide) = k1(isUorSide) .* k2(isUorSide) .* k3(isUorSide) .* leMM(isUorSide) ./...
    (EFrpMPA(isUorSide).*tFrpMM(isUorSide).*ep_fu(isUorSide));
kv(isUorSide & (kv<0) ) = 0;
kv(isUorSide & (kv>FRP_STRN_LMT_REL) ) = FRP_STRN_LMT_REL;

% failure strain Eq. (11-6)
ep_fe(isUorSide) = kv(isUorSide) .* ep_fu(isUorSide);
ep_fe( (isUorSide) & (ep_fe>FRP_STRN_LMT_ABS) ) = FRP_STRN_LMT_ABS;

% keep track of side warning for side bonding (2Le>hfrp):
indx = 2*leMM>dfvMM;
sidewarning = zeros(nCase, 1);
sidewarning(indx) = 1;

%% Wrapping

% ep_fe(isW) = FRP_STRN_LMT_ABS;
% isWandOver = (isW) & (ep_fe > FRP_STRN_LMT_REL*ep_fu);
% ep_fe( isWandOver ) = FRP_STRN_LMT_REL * ep_fu( isWandOver );

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
ep_fe( isW&ep_fe>FRP_STRN_LMT_ABS ) = FRP_STRN_LMT_ABS;

%% FRP contribution to shear resistance Eq. (11-3)

f_feMPA = ep_fe .* EFrpMPA;
shearFrpKN = 2*tFrpMM .* widthFrpMM .* f_feMPA .*...
        (sin(betaRAD)+cos(betaRAD)) .* dfvMM ./ sFrpMM /1e3;
shearFrpKN( shearFrpKN<0 ) = 0;

%% ACI 318-11: Steel Contribution

fsMPA( fsMPA>fsLimit ) = fsLimit;

areaVMM2 = 2*pi/4*sdMM.^2;
shearSteelKN = areaVMM2 .* fsMPA .* dBeamMM ./ ssMM /1e3;
shearSteelKN( ssMM == 0 ) = 0;
shearSteelKN( shearSteelKN<0 ) = 0;

%% contribution of concrete

areaV_Min = 0.062*sqrt(fcMPA) .* bBeamMM .* ssMM ./ fsMPA;
tmp = (0.35*bBeamMM .* sdMM) ./ fsMPA;
areaV_Min( areaV_Min < tmp ) = tmp(areaV_Min < tmp); 

sqrt_fc( areaVMM2<areaV_Min & sqrt_fc>sqrtFcLimit ) = sqrtFcLimit;

shearConcreteKN = 0.17*sqrt_fc .* bBeamMM .* dBeamMM / 1e3;
shearConcreteKN(shearConcreteKN<0) = 0;

%% FRP and steel interaction
mu = zeros(nCase,1);
isWithSteel = areaVMM2~=0;
mu(isSide&isWithSteel) = shearSteelKN(isSide&isWithSteel) ./ shearFrpKN(isSide&isWithSteel) ;
mu(isSide&(shearFrpKN==0) ) = 0;
kfrp = -0.2*mu+1; kfrp(kfrp<0) = 0;
shearFrpKN(isSide) = kfrp(isSide) .* shearFrpKN(isSide);
shearFrpKN( isSide&shearFrpKN<0 ) = 0;

%% over reinforcement check

shearReinforceKN = shearSteelKN + shearFrpKN;
ro = shearSteelKN ./ shearReinforceKN;
isOverReinforce = (shearReinforceKN > 0.66*sqrt_fc.*bBeamMM .* dBeamMM / 1e3);
shearReinforceKN( isOverReinforce ) = 0.66*sqrt_fc(isOverReinforce).*...
                       bBeamMM(isOverReinforce).*dBeamMM(isOverReinforce)/1e3;
switch FLAG
    case 'MODEL_ERROR'
        shearTotalKN = shearConcreteKN + shearReinforceKN;
    case 'DESIGN_VALUE'
        shearTotalKN = PHI * (shearConcreteKN + ro.*shearReinforceKN +...
                              factorFrp .* (1-ro).*shearReinforceKN);
    otherwise
end

return
end