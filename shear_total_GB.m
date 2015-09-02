function [shearTotalKN, isOverReinforce,...
         ro, shearReinforceKN] = shear_total_GB(FLAG, factorFrp)
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
N_STEEL_LEG = 2; % number of steel legs
N_FRP_LEG = 2; % number of FRP legs
PHI0_SIDE = 1.0;
PHI0_U = 1.3;
S2D_LOWER = 1.5;
S2D_UPPER = 3.0;

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
        s2d = S2D_TEST_ARRAY;
        % Concrete mean properties
        fcMPA = FC_TEST_ARRAY_MPA;
        fcuMPA = fcMPA / RO_CYLINDE_2_CUBE;
        fctMPA = 0.395 .* fcuMPA.^0.55;
        % steel mean properties
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
        fcuLimit1 = 50*FC_BIAS;
        fcuLimit2 = 80*FC_BIAS;
        fcuBrittleLimit1 = 40*FC_BIAS;
        fcuBrittleLimit2 = 80*FC_BIAS;
        % ratio of in-situ to cylinder strength
        roInsitu = 1.0;  
    case 'DESIGN_VALUE'
        nCase = length(FC_DESIGN_ARRAY_MPA);
        % constant partial safety factors
        gammaConcrete = 1.40;
        gammaSteel = 1.10;        
        % geometrical properties
        hBeamMM = H_DESIGN_ARRAY_MM;
        dBeamMM = D_DESIGN_ARRAY_MM;
        bBeamMM = B_DESIGN_ARRAY_MM;
        dFrpMM = DFRP_DESIGN_ARRAY_MM;
        dFrpTopMM = DFRP_TOP_DESIGN_ARRAY_MM; 
        s2d = S2D_DESIGN_ARRAY;
        % Concrete characteristic properties
        fcMPA = FC_DESIGN_ARRAY_MPA;
        fcuMPA = fcMPA / RO_CYLINDE_2_CUBE; 
        fctMPA = 0.395 .* fcuMPA.^0.55 .* (1-1.645*FCT_COV).^(0.45);
        % steel characteristic properties
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
        fcuLimit1 = 50;
        fcuLimit2 = 80;
        fcuBrittleLimit1 = 40;
        fcuBrittleLimit2 = 80;
        % ratio of in-situ to cylinder strength
        roInsitu = 0.88;         
    otherwise
end

zt = dFrpTopMM;
zb = (dBeamMM-(hBeamMM-dFrpMM))-TOP_DUMMY_HEIGHT*dBeamMM;
hFrpEff = zb-zt;
betaFrpRAD = betaFrpDEG/180*pi;

prism2Cube = 0.76*ones(nCase, 1);
isHighLimit1 = (fcuMPA > fcuLimit1) & (fcuMPA < fcuLimit2);
isHighLimit2 = (fcuMPA >= fcuLimit2);
prism2Cube(isHighLimit1) = (0.76-0.82) / (fcuLimit1-fcuLimit2) *...
                             (fcuMPA(isHighLimit1)-fcuLimit2) + 0.82;
prism2Cube(isHighLimit2) = 0.82;

alphaBrittle = ones(nCase, 1);
isHighLimit1 = (fcuMPA > fcuBrittleLimit1) & (fcuMPA < fcuBrittleLimit2);
isHighLimit2 = (fcuMPA >= fcuBrittleLimit2);
alphaBrittle(isHighLimit1) = (1.0-0.87) / (fcuBrittleLimit1-fcuBrittleLimit2) *...
                             (fcuMPA(isHighLimit1)-fcuBrittleLimit2) + 0.87;
alphaBrittle(isHighLimit2) = 0.87;

fctMPA = roInsitu * alphaBrittle .* fctMPA;
fcMPA = roInsitu * prism2Cube .* alphaBrittle .* fcuMPA;

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

switch FLAG
    case 'MODEL_ERROR'
        shearTotalKN = shearConcreteKN+ shearSteelKN + shearFrpKN;
    case 'DESIGN_VALUE'
        psi_f = 0.75*ones(nCase, 1);
        psi_f(isSide) = 0.5;
        shearTotalKN = shearConcreteKN+ shearSteelKN + psi_f.*shearFrpKN;
end
ro = shearSteelKN ./ (shearSteelKN + shearFrpKN);

% elimination of diagonal-compression failure
betaC = ones(nCase, 1);
isHighLimit1 = (fcuMPA > fcuLimit1) & (fcuMPA < fcuLimit2);
isHighLimit2 = (fcuMPA >= fcuLimit2);
betaC(isHighLimit1) = (1.0-0.8) / (fcuLimit1-fcuLimit2) * (fcuMPA(isHighLimit1)-fcuLimit2) + 0.8;
betaC(isHighLimit2) = 0.8;
isOverReinforce = shearTotalKN > 0.2*betaC.*fcMPA .* bBeamMM .* dBeamMM *1e-3;
shearTotalKN(isOverReinforce) = 0.2 * betaC(isOverReinforce) .* fcMPA(isOverReinforce) .*...
                                bBeamMM(isOverReinforce) .* dBeamMM(isOverReinforce) * 1e-3;

shearReinforceKN = shearTotalKN - shearConcreteKN;

return
end