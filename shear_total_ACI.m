function [shearTotalKN, isOverReinforce,...
         ro, shearReinforceKN] = shear_total_ACI(FLAG, factorFrp)
% total shear resistance based on ACI 440.2R-08
PHI = 0.75; % resistance reduce factor (ACI 318-11 9.3.2.3)

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
k1 = zeros(nCase, 1);
k2 = zeros(nCase, 1);
kv = zeros(nCase, 1);
ep_fe = zeros(nCase, 1);
% active bond length Eq. (11-8)
leMM(isUorSide) = 23300 ./ (tFrpMM(isUorSide) .* EFrpMPA(isUorSide)).^0.58;
% parameter k1 Eq. (11-9)              
k1(isUorSide) = (fcMPA(isUorSide) ./ 27) .^ (2/3);              
% parameter k2 Eq. (11-10) 
k2(isU) = ( dfvMM(isU)-leMM(isU) ) ./ dfvMM(isU);
k2(isSide) = ( dfvMM(isSide)-2*leMM(isSide) ) ./ dfvMM(isSide);
% parameter kv Eq. (11-7)
kv(isUorSide) = k1(isUorSide) .* k2(isUorSide) .* leMM(isUorSide) ./...
                ( 11900*ep_fu(isUorSide) );
kv(isUorSide & (kv<0) ) = 0;
kv(isUorSide & (kv>0.75) ) = 0.75;
% switch FLAG
%     case 'MODEL_ERROR'
%         kv(isUorSide & (kv>0.75) ) = 0.75;
%     case 'MONTE_CARLO'
%         nLargeKv = sum( isUorSide & (kv>0.75) );
%         meanKv = mean( kv(isUorSide & (kv<=0.75)) );
%         stdKv = std( kv(isUorSide & (kv<=0.75)) );
%         covKv = stdKv / meanKv;
%         kv( isUorSide & (kv>0.75) ) = normrnd(0.75, 0.75*covKv, nLargeKv, 1);        
%     otherwise
% end
% failure strain Eq. (11-6)
ep_fe(isUorSide) = kv(isUorSide) .* ep_fu(isUorSide);
ep_fe( (isUorSide) & (ep_fe>0.004) ) = 0.004;

%% Wrapping

ep_fe(isW) = 0.004;
isWandOver = (isW) & (ep_fe > 0.75*ep_fu);
ep_fe( isWandOver ) = 0.75 * ep_fu( isWandOver );

%% FRP contribution to shear resistance Eq. (11-3)

f_feMPA = ep_fe .* EFrpMPA;
shearFrpKN = 2*tFrpMM .* widthFrpMM .* f_feMPA .*...
        (sin(betaRAD)+cos(betaRAD)) .* dfvMM ./ sFrpMM /1e3;
shearFrpKN( shearFrpKN<0 ) = 0;

% ACI 318-11: Steel Contribution

fsMPA( fsMPA>fsLimit ) = fsLimit;

areaVMM2 = 2*pi/4*sdMM.^2;
shearSteelKN = areaVMM2 .* fsMPA .* dBeamMM ./ ssMM /1e3;
shearSteelKN( ssMM == 0 ) = 0;
shearSteelKN( shearSteelKN<0 ) = 0;

% contribution of concrete

areaV_Min = 0.062*sqrt(fcMPA) .* bBeamMM .* ssMM ./ fsMPA;
tmp = (0.35*bBeamMM .* sdMM) ./ fsMPA;
areaV_Min( areaV_Min < tmp ) = tmp(areaV_Min < tmp); 

sqrt_fc( areaVMM2<areaV_Min & sqrt_fc>sqrtFcLimit ) = sqrtFcLimit;

shearConcreteKN = 0.17*sqrt_fc .* bBeamMM .* dBeamMM / 1e3;
shearConcreteKN(shearConcreteKN<0) = 0;

% over reinforcement check

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