function [shearTotalKN, isOverReinforce,...
         ro, shearReinforceKN] = shear_total_HK(FLAG, factorFrp)
% Determine the contribution of FRP, New HK guidelines
% Considering the steel-FRP interaction

nFactorFrp = length(factorFrp);
if nFactorFrp == 2
    gammaFrp = factorFrp(1);
    gammaBond = factorFrp(2);
elseif nFactorFrp == 1
    gammaFrp = factorFrp;
    gammaBond = factorFrp;
end

% important constants
RO_CYLINDE_2_CUBE = 0.8; % strength ratio of cylinders to cubes
TOP_DUMMY_HEIGHT = 0.1; % ratio of top dummy height to total depth
MAX_FRP = 0.8; % maximum tolerant FRP strength
N_STEEL_LEG = 2; % number of steel legs
N_FRP_LEG = 2; % number of FRP legs
RO_BEND_STEEL = 0.01; % ratio of flexure steel

load tmpdata.mat

switch FLAG
    case 'MODEL_ERROR'
        nCase = length(FC_TEST_ARRAY_MPA);
        % constant partial safety factors
        gammaConcrete = 1.00;
        gammaSteel = 1.00;
        % geometrical properties
        hBeamMM = H_TEST_ARRAY_MM;
        dBeamMM = D_TEST_ARRAY_MM;
        bBeamMM = B_TEST_ARRAY_MM;
        dFrpMM = DFRP_TEST_ARRAY_MM;
        dFrpTopMM = DFRP_TOP_TEST_ARRAY_MM;
        % Concrete mean properties
        fcMPA = FC_TEST_ARRAY_MPA;
        fcuMPA = fcMPA / RO_CYLINDE_2_CUBE;
        % steel mean properties
        barType = BAR_TYPE_TEST_ARRAY;
        fsMPA = FS_TEST_ARRAY_MPA;
        sdMM = SD_TEST_ARRAY_MM;
        ssMM = SS_TEST_ARRAY_MM;
        % FRP mean properties
        frpForm = FRP_FORM_TEST_ARRAY;
        fFrpMPA = F_FRP_TEST_ARRAY_MPA;
        EFrpMPA = E_FRP_TEST_ARRAY_MPA;
        betaFrpDEG = BETA_TEST_ARRAY_DEG;
        tFrpMM = T_FRP_TEST_ARRAY_MM;
        widthFrpMM = W_FRP_TEST_ARRAY_MM;
        sFrpMM = S_FRP_TEST_ARRAY_MM;
        % limit values
        fcuLowLimit4Steel = 40*FC_BIAS;
        fcuLowLimit = 25*FC_BIAS;
        fcuHighLimit = 80*FC_BIAS;
        vrLimit = 0.4;
        sqrtFcuLimit1 = sqrt(fcuMPA)*FC_BIAS;
        sqrtFcuLimit2 = (7*1.25)*FC_BIAS;
        % ratio of in-situ to cylinder strength
        roInsitu = 1.00;        
    case 'DESIGN_VALUE'
        nCase = length(FC_DESIGN_ARRAY_MPA);
        % constant partial safety factors
        gammaConcrete = 1.25;
        gammaSteel = 1.15;        
        % geometrical properties
        hBeamMM = H_DESIGN_ARRAY_MM;
        dBeamMM = D_DESIGN_ARRAY_MM;
        bBeamMM = B_DESIGN_ARRAY_MM;
        dFrpMM = DFRP_DESIGN_ARRAY_MM;
        dFrpTopMM = DFRP_TOP_DESIGN_ARRAY_MM;        
        % Concrete characteristic properties
        fcMPA = FC_DESIGN_ARRAY_MPA;
        fcuMPA = fcMPA / RO_CYLINDE_2_CUBE;                
        % steel characteristic properties
        barType = BAR_TYPE_DESIGN_ARRAY;
        fsMPA = FS_DESIGN_ARRAY_MPA;
        sdMM = SD_DESIGN_ARRAY_MM;
        ssMM = SS_DESIGN_ARRAY_MM;
        % FRP characteristic properties
        frpForm = FRP_FORM_DESIGN_ARRAY;
        fFrpMPA = F_FRP_DESIGN_ARRAY_MPA;
        EFrpMPA = E_FRP_DESIGN_ARRAY_MPA;
        betaFrpDEG = BETA_DESIGN_ARRAY_DEG;
        tFrpMM = T_FRP_DESIGN_ARRAY_MM;
        widthFrpMM = W_FRP_DESIGN_ARRAY_MM;
        sFrpMM = S_FRP_DESIGN_ARRAY_MM;
        % limit values
        fcuLowLimit4Steel = 40;
        fcuLowLimit = 25;
        fcuHighLimit = 80;
        vrLimit = 0.4;
        sqrtFcuLimit1 = sqrt(fcuMPA) / gammaConcrete;
        sqrtFcuLimit2 = (7*1.25) / gammaConcrete;
        % ratio of in-situ to cylinder strength
        roInsitu = 0.8375;          
    otherwise
end

zt = dFrpTopMM;
zb = (dBeamMM-(hBeamMM-dFrpMM))-TOP_DUMMY_HEIGHT*dBeamMM;
hFrpEff = zb-zt;
betaFrpRAD = betaFrpDEG/180*pi;

%% FRP contribution (without FRP-steel interaction)
isSide = (frpForm == 1);
isU = (frpForm == 2);
isUorSide = (frpForm == 1) | (frpForm == 2);
isW = (frpForm == 3);

roFrp = widthFrpMM ./ (sFrpMM.*sin(betaFrpRAD));
tmp = (2-roFrp)./(1+roFrp);
tmp( tmp<0 ) = 0;
betaW = sqrt(tmp);

lMax = zeros(nCase, 1);
lMax(isSide) = hFrpEff(isSide) ./ 2 ./ sin(betaFrpRAD(isSide));
lMax(isU) = hFrpEff(isU) ./ sin(betaFrpRAD(isU));

lEff = sqrt( EFrpMPA.*tFrpMM ./ sqrt(fcMPA) );
roLmax2Leff = lMax ./ lEff;

betaL = zeros(nCase, 1);
isLongEnough = roLmax2Leff >= 1;
betaL( isLongEnough ) = 1;
betaL( ~isLongEnough ) = sin(pi*roLmax2Leff(~isLongEnough)/2);

% Partial safety factor for debonding strength
sDebond = 0.315*betaW.*betaL.*sqrt(EFrpMPA.*sqrt(fcMPA)./tFrpMM);

failMode = zeros(nCase, 1);
sFrpDesignMax = zeros(nCase, 1);

sFrpDesignMax(isW) = MAX_FRP * fFrpMPA(isW)/gammaFrp;
failMode(isW) = 1;
sFrpDesignMax(isUorSide) = min( sDebond(isUorSide)/gammaBond, MAX_FRP * fFrpMPA(isUorSide)/gammaFrp);
isPossibleDebond = sDebond/gammaBond <= MAX_FRP * fFrpMPA/gammaFrp;
failMode(isUorSide & (~isPossibleDebond)) = 1;
failMode(isUorSide & isPossibleDebond) = 2;
isDebond = (failMode==2);
isRupture = (failMode==1);

distributeFrp = zeros(nCase, 1);
zeta = zt./zb;
distributeFrp(isRupture) = (1+zeta(isRupture))/2;
distributeFrp( (isDebond & isLongEnough) ) = 1 - (pi-2)./(pi*roLmax2Leff(isDebond & isLongEnough));
distributeFrp( isDebond & (~isLongEnough) ) = 2/pi./roLmax2Leff(isDebond & (~isLongEnough)) .*...
                                               (1-cos(pi/2*roLmax2Leff(isDebond & (~isLongEnough)))) ./....
                                               sin(pi/2*roLmax2Leff(isDebond & (~isLongEnough)));
fFrpDesignEff = sFrpDesignMax .* distributeFrp;
shearFrpKN = N_FRP_LEG*fFrpDesignEff.*tFrpMM.*widthFrpMM.*hFrpEff.*(sin(betaFrpRAD)+cos(betaFrpRAD))./sFrpMM/1000;
shearFrpKN( shearFrpKN < 0 ) = 0;

%% Steel contribution (Building Department HK 2013)
isWithSteel = barType~=0;
areaVMM2 = N_STEEL_LEG*pi/4*sdMM.^2;
shearSteelKN = zeros(nCase, 1);
fsDesignMPA = fsMPA/gammaSteel;
shearSteelKN(isWithSteel) = 1e-3*areaVMM2(isWithSteel) .* fsDesignMPA(isWithSteel) .*...
                            dBeamMM(isWithSteel) ./ ssMM(isWithSteel);

%% steel-FRP interaction (Chen's dissertation)
phiSteel = zeros(nCase, 1);
A = zeros(nCase, 1);
mu = zeros(nCase, 1);
kfrp = ones(nCase,1);

isDeformed = barType == 1;
isRound = barType == 2;

phiSteel(isSide & isDeformed) = 1e5 ./ (sdMM(isSide & isDeformed).^1.13 .* fsMPA(isSide & isDeformed).^ 1.71);
A(isSide & isDeformed) = phiSteel(isSide & isDeformed) .* (2.045*2*sin(betaFrpRAD(isSide & isDeformed)).*roLmax2Leff(isSide & isDeformed)+3.24);
phiSteel(isSide & isRound) = 1e5 ./ (sdMM(isSide & isRound).^0.834 .* fsMPA(isSide & isRound).^ 1.88);
A(isSide & isRound) = phiSteel(isSide & isRound) .* (1.01*2*sin(betaFrpRAD(isSide & isRound)).*roLmax2Leff(isSide & isRound)+2.13);
mu(isSide&isWithSteel) = shearSteelKN(isSide&isWithSteel) ./ shearFrpKN(isSide&isWithSteel) ;
mu(isSide&(shearFrpKN == 0) ) = 0;
kfrp( isSide&isWithSteel ) = A(isSide&isWithSteel)./ (A(isSide&isWithSteel)+mu(isSide&isWithSteel));

% Update FRP contribution (consider FRP-steel interaction)
shearFrpKN(isSide) = kfrp(isSide) .* shearFrpKN(isSide);
shearFrpKN( isSide&shearFrpKN<0 ) = 0;

%% concrete contribution (Building Department HK 2013)
tmp1 = (400./ dBeamMM).^(1/4);

vr = zeros(nCase,1);
isNormal = (fcuMPA>fcuLowLimit4Steel) & (fcuMPA<fcuHighLimit);
is2Low = fcuMPA <= fcuLowLimit4Steel;
is2High = fcuMPA >= fcuHighLimit;
vr( is2Low ) = vrLimit;
vr( isNormal ) = vrLimit*(fcuMPA(isNormal)/fcuLowLimit4Steel).^(2/3);
vr( is2High ) = vrLimit*(fcuHighLimit/fcuLowLimit4Steel).^(2/3);
areaSteelMin = vr.* bBeamMM .* ssMM ./ (fsMPA/gammaSteel);
tmp1( tmp1<0.8*roInsitu ) = 0.8*roInsitu; 
tmp1( areaVMM2 >= areaSteelMin & tmp1<1 ) = 1;

vcDesign = 0.79*(100*RO_BEND_STEEL).^(1/3).*tmp1 .* 1/gammaConcrete;
shearConcreteKN = vcDesign .* bBeamMM .* dBeamMM *1e-3;

isNormal = fcuMPA>fcuLowLimit & fcuMPA<fcuHighLimit;
is2High = fcuMPA>=fcuHighLimit;

shearConcreteKN( isNormal ) = shearConcreteKN( isNormal ).* ...
                              (fcuMPA(isNormal)/fcuLowLimit).^(1/3);
shearConcreteKN( is2High ) = shearConcreteKN( is2High ).* ...
                              (fcuHighLimit/fcuLowLimit)^(1/3);
shearConcreteKN( shearConcreteKN<0 ) = 0;

shearTotalKN = shearFrpKN + shearSteelKN + shearConcreteKN;

ro = shearSteelKN ./ (shearSteelKN + shearFrpKN);

sqrtFcuLimitDesign = min(sqrtFcuLimit1, sqrtFcuLimit2);
isOverReinforce = shearTotalKN*1e3./(bBeamMM.*dBeamMM) > sqrtFcuLimitDesign;
shearTotalKN( isOverReinforce ) = sqrtFcuLimitDesign(isOverReinforce) .* bBeamMM(isOverReinforce) .* dBeamMM(isOverReinforce) / 1000;

shearReinforceKN = shearTotalKN - shearConcreteKN;

return
end