function shearTotalKN = shear_total_GB_MC(iDesignCase, hSmp, betaSmp, fFrpSmp, fcSmp, fctSmpNoBrittleFactor, fsSmp)
% Determine the contribution of FRP, New HK guidelines
% Considering the steel-FRP interaction

% important constants
RO_CYLINDE_2_CUBE = 0.8; % strength ratio of cylinders to cubes
TOP_DUMMY_HEIGHT = 0.1; % ratio of top dummy height to total depth
N_STEEL_LEG = 2; % number of steel legs
N_FRP_LEG = 2; % number of FRP legs
PHI0_SIDE = 1.0;
PHI0_U = 1.3;
S2D_LOWER = 1.5;
S2D_UPPER = 3.0;

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
s2d = S2D_DESIGN_ARRAY(iDesignCase)*ones(nCase, 1);
% Concrete properties
fcMPA = fcSmp;
fcuMPA = fcMPA / RO_CYLINDE_2_CUBE;
fctMPA = fctSmpNoBrittleFactor;
% sqrtFcuSmp = sqrtFcSmp./sqrt(0.8);
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
fcuLimit1 = 50*FC_BIAS;
fcuLimit2 = 80*FC_BIAS;
fcuBrittleLimit1 = 40*FC_BIAS;
fcuBrittleLimit2 = 80*FC_BIAS;
% ratio of in-situ to cylinder strength
roInsitu = 0.88;

zt = dFrpTopMM;
zb = (dBeamMM-(hBeamMM-dFrpMM))-TOP_DUMMY_HEIGHT*dBeamMM;
hFrpEff = zb-zt;
betaFrpRAD = betaFrpDEG/180*pi;

% % the following block is likely to amplify the std a litte bit when fcu
% % is high, but the difference is little
% alphaBrittle = ones(nCase, 1);
% isHighLimit1 = ( fcuMPA > fcuBrittleLimit1) & (fcuMPA < fcuBrittleLimit2);
% isHighLimit2 = (fcuMPA >= fcuBrittleLimit2);
% alphaBrittle(isHighLimit1) = (1.0-0.87) / (fcuBrittleLimit1-fcuBrittleLimit2) *...
%                              (fcuMPA(isHighLimit1)-fcuBrittleLimit2) + 0.87;
% alphaBrittle(isHighLimit2) = 0.87;
% fctMPA = alphaBrittle .* fctMPA;
prism2CubeMean = 0.76;
isHighLimit1 = ( mean(fcuMPA) > fcuLimit1) & (mean(fcuMPA) < fcuLimit2);
isHighLimit2 = ( mean(fcuMPA) >= fcuLimit2);
prism2CubeMean(isHighLimit1) = (0.76-0.82) / (fcuLimit1-fcuLimit2) *...
                             (mean(fcuMPA)-fcuLimit2) + 0.82;
prism2CubeMean( isHighLimit2 ) = 0.82;   

alphaBrittleMean = 1;
isHighLimit1 = ( mean(fcuMPA) > fcuBrittleLimit1) & (mean(fcuMPA) < fcuBrittleLimit2);
isHighLimit2 = ( mean(fcuMPA) >= fcuBrittleLimit2);
alphaBrittleMean(isHighLimit1) = (1.0-0.87) / (fcuBrittleLimit1-fcuBrittleLimit2) *...
                             (mean(fcuMPA)-fcuBrittleLimit2) + 0.87;
alphaBrittleMean( isHighLimit2 ) = 0.87;

fctMPA = roInsitu * alphaBrittleMean .* fctMPA;
fcMPA = roInsitu * prism2CubeMean .* alphaBrittleMean .* fcuMPA;

%% concrete and steel contribution

s2d( s2d<S2D_LOWER ) = S2D_LOWER;
s2d( s2d>S2D_UPPER ) = S2D_UPPER;
shearConcreteKN = 1.75 ./ (s2d + 1) .* fctMPA .* bBeamMM.*dBeamMM * 1e-3 / gammaConcrete;
areaVMM2 = N_STEEL_LEG * pi* sdMM.^2/4;
shearSteelKN = fsMPA .* areaVMM2 ./ ssMM .* dBeamMM * 1e-3 / gammaSteel;
shearSteelKN( ssMM==0 ) = 0;

%% FRP contribution
shearFrpKN = zeros(nCase, 1);
isSide = (frpForm == 1);
isU = (frpForm == 2);
isUorSide = (frpForm == 1) | (frpForm == 2);
isW = (frpForm == 3);

lambdaEf = N_FRP_LEG * widthFrpMM(isW) .* tFrpMM(isW) ./...
           (bBeamMM(isW).*sFrpMM(isW).*sin(betaFrpRAD(isW))) .*...
           EFrpMPA(isW) ./ (fctMPA(isW)/gammaConcrete);
eFrpEff = 8 ./ ( sqrt(lambdaEf)+10 ) .* (fFrpMPA(isW)/gammaFrp) ./ EFrpMPA(isW);  
eFrpEff(lambdaEf<=0) = 0;
sFrpEff = min( (fFrpMPA(isW)/gammaFrp), EFrpMPA(isW).*eFrpEff);
shearFrpKN(isW) = 1e-3 * N_FRP_LEG .* widthFrpMM(isW) .* tFrpMM(isW) ./...
                  sFrpMM(isW) .* sFrpEff .* hFrpEff(isW) .* ...
                  ( sin(betaFrpRAD(isW)) + cos(betaFrpRAD(isW)) );              

betaW = sqrt( (2.25-widthFrpMM(isUorSide) ./ (sFrpMM(isUorSide).*sin(betaFrpRAD(isUorSide)))) ./...
              (1.25+widthFrpMM(isUorSide) ./ (sFrpMM(isUorSide).*sin(betaFrpRAD(isUorSide)))) );
bondStrength = 1.2 * betaW .* (fctMPA(isUorSide)/gammaConcrete) / gammaBond;
debondFactor = sin(betaFrpRAD(isUorSide)) .* sqrt(EFrpMPA(isUorSide).*tFrpMM(isUorSide)) ./...
              ( sin(betaFrpRAD(isUorSide)).*sqrt(EFrpMPA(isUorSide).*tFrpMM(isUorSide)) +...
              0.3*hFrpEff(isUorSide).*(fctMPA(isUorSide)/gammaConcrete) ) ;
shearFrpKN(isUorSide) = 1e-3 * debondFactor .* bondStrength .* widthFrpMM(isUorSide) .*...
             hFrpEff(isUorSide).^2 ./ sFrpMM(isUorSide) .*...
             (sin(betaFrpRAD(isUorSide)) + cos(betaFrpRAD(isUorSide)));

shearFrpKN(isSide) = PHI0_SIDE * shearFrpKN(isSide);
shearFrpKN(isU) = PHI0_U * shearFrpKN(isU);

shearTotalKN = shearConcreteKN+ shearSteelKN + shearFrpKN;

% elimination of diagonal-compression failure

% betaC = ones(nCase, 1);
% isHighLimit1 = (fcuMPA > fcuLimit1) & (fcuMPA < fcuLimit2);
% isHighLimit2 = (fcuMPA >= fcuLimit2);
% betaC(isHighLimit1) = (1.0-0.8) / (fcuLimit1-fcuLimit2) * (fcuMPA(isHighLimit1)-fcuLimit2) + 0.8;
% betaC(isHighLimit2) = 0.8;
% isOverReinforce = shearTotalKN > 0.25*betaC.*fcuMPA .* bBeamMM .* dBeamMM *1e-3;
% shearTotalKN(isOverReinforce) = 0.25 * betaC(isOverReinforce) .* fcMPA(isOverReinforce) .*...
%                                 bBeamMM(isOverReinforce) .* dBeamMM(isOverReinforce) * 1e-3;
betaCMean = 1;
isHighLimit1 = (mean(fcuMPA) > fcuLimit1) & (mean(fcuMPA) < fcuLimit2);
isHighLimit2 = (mean(fcuMPA) >= fcuLimit2);
betaCMean(isHighLimit1) = (1.0-0.8) / (fcuLimit1-fcuLimit2) * (mean(fcuMPA)-fcuLimit2) + 0.8;
betaCMean(isHighLimit2) = 0.8;
isOverReinforce = shearTotalKN > 0.25*betaCMean.*fcuMPA .* bBeamMM .* dBeamMM *1e-3;
shearTotalKN(isOverReinforce) = 0.25 * betaCMean .* fcMPA(isOverReinforce) .*...
                                bBeamMM(isOverReinforce) .* dBeamMM(isOverReinforce) * 1e-3;

return
end