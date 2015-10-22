function shearTotalKN = shear_total_HK_MC(iDesignCase, hSmp, betaSmp, fFrpSmp, fcSmp, sqrtFcSmp, fsSmp)
% Determine the contribution of FRP, New HK guidelines
% Considering the steel-FRP interaction

% important constants
RO_CYLINDE_2_CUBE = 0.8; % strength ratio of cylinders to cubes
TOP_DUMMY_HEIGHT = 0.1; % ratio of top dummy height to total depth
MAX_FRP = 0.8; % maximum tolerant FRP strength
N_STEEL_LEG = 2; % number of steel legs
N_FRP_LEG = 2; % number of FRP legs
RO_BEND_STEEL = 0.01; % ratio of flexure steel

load tmpdata.mat

% read variables
nCase = length(hSmp);
% constant partial safety factors
gammaConcrete = 1.00;
gammaSteel = 1.00;
gammaFrp = 1.00;
gammaBond = 1.00;
% geometrical properties
hBeamMM = hSmp;
dBeamMM = hBeamMM - AS_DESIGN_ARRAY_MM(iDesignCase);
bBeamMM = hSmp/H2B_DESIGN_ARRAY(iDesignCase);
dFrpMM = hSmp;
dFrpTopMM = DFRP_TOP_DESIGN_ARRAY_MM(iDesignCase)*ones(nCase, 1);
% Concrete properties
fcMPA = fcSmp;
fcuMPA = fcMPA / RO_CYLINDE_2_CUBE;
% sqrtFcuSmp = sqrtFcSmp./sqrt(0.8);
% steel properties
barType = BAR_TYPE_DESIGN_ARRAY(iDesignCase)*ones(nCase, 1);
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
fcuLowLimit4Steel = 40*FC_BIAS;
fcuLowLimit = 25*FC_BIAS;
fcuHighLimit = 80*FC_BIAS;
vrLimit = 0.4;
sqrtFcuLimit1Mean = sqrt( mean(fcuMPA)/FC_BIAS )*FC_BIAS;
sqrtFcuLimit2Mean = (7*1.25)*FC_BIAS;
sqrtFcuLimit1 = normrnd(sqrtFcuLimit1Mean, sqrtFcuLimit1Mean*FC_COV, nCase, 1);
sqrtFcuLimit2 = normrnd(sqrtFcuLimit2Mean, sqrtFcuLimit2Mean*FC_COV, nCase, 1);
% ratio of in-situ to cylinder strength
roInsitu = 0.8375;

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
roLmax2Leff( imag(lEff)~=0 ) = 0;

betaL = zeros(nCase, 1);
isLongEnough = roLmax2Leff >= 1;
betaL( isLongEnough ) = 1;
betaL( ~isLongEnough ) = sin(pi*roLmax2Leff(~isLongEnough)/2);

% Partial safety factor for debonding strength
sDebond = 0.427*betaW.*betaL.*sqrt(EFrpMPA.*sqrtFcSmp./tFrpMM);
sDebond( imag(sDebond)~=0 ) = 0;

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
% kfrp = ones(size(kfrp));
% kfrp = -0.2*mu+1; kfrp(kfrp<0) = 0;
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
% determine variation of tmp1
tmp1Mean = mean(tmp1);
tmp1Mean( tmp1Mean<0.8*roInsitu ) = 0.8*roInsitu; 
tmp1Mean( mean(areaVMM2) >= mean(areaSteelMin) & tmp1Mean<1 ) = 1;
tmp1Cov = FCT_COV;
tmp1 = normrnd(tmp1Mean, tmp1Mean*tmp1Cov, nCase, 1);
tmp1(tmp1<0) = 0;
% isNormal = (tmp1>=1) | (areaVMM2<areaSteelMin & tmp1>=0.8 & tmp1<1);
% smpNormal = tmp1(isNormal);
% covNormal = std(smpNormal)/mean(smpNormal);
% covNotNormal = min(FCT_COV, covNormal);
% isSmallType1 = (tmp1<0.8);
% nSmallType1 = sum(isSmallType1);
% isSmallType2 = (areaVMM2 >= areaSteelMin & tmp1<1);
% nSmallType2 = sum(isSmallType2);
% tmp1( isSmallType1 ) = normrnd(0.8, 0.8*covNotNormal, nSmallType1, 1); 
% tmp1( isSmallType2 ) = normrnd(1, covNotNormal, nSmallType2, 1);

vcDesign = 0.79*(100*RO_BEND_STEEL).^(1/3).*tmp1 .* 1/gammaConcrete;
shearConcreteKN = vcDesign .* bBeamMM .* dBeamMM *1e-3;

isNormal = fcuMPA>fcuLowLimit & fcuMPA<fcuHighLimit;
is2High = fcuMPA>=fcuHighLimit;

shearConcreteKN( isNormal ) = shearConcreteKN( isNormal ).* ...
                              ( mean(fcuMPA) /fcuLowLimit).^(1/3);
shearConcreteKN( is2High ) = shearConcreteKN( is2High ).* ...
                              (fcuHighLimit/fcuLowLimit)^(1/3);
shearConcreteKN( shearConcreteKN<0 ) = 0;

shearTotalKN = shearFrpKN + shearSteelKN + shearConcreteKN;

sqrtFcuLimitDesign = min(sqrtFcuLimit1, sqrtFcuLimit2);
isOverReinforce = shearTotalKN*1e3./(bBeamMM.*dBeamMM) > sqrtFcuLimitDesign;
shearTotalKN( isOverReinforce ) = sqrtFcuLimitDesign(isOverReinforce) .* bBeamMM(isOverReinforce) .* dBeamMM(isOverReinforce) / 1000;

return
end