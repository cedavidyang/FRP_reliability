function shearTotalKN = shear_total_ACI_MC(iDesignCase, hSmp, betaSmp, fFrpSmp, fcSmp, sqrtFcSmp, fsSmp)
% total shear resistance based on ACI 440.2R-08

load tmpdata.mat

% read variables        
nCase = length(hSmp);
% geometrical properties
% dBeamMM = D_DESIGN_ARRAY_MM(iDesignCase)*ones(nCase, 1);
% bBeamMM = B_DESIGN_ARRAY_MM(iDesignCase)*ones(nCase, 1);
dBeamMM = hSmp - AS_DESIGN_ARRAY_MM(iDesignCase);
bBeamMM = hSmp/H2B_DESIGN_ARRAY(iDesignCase);
dFrpTopMM = DFRP_TOP_DESIGN_ARRAY_MM(iDesignCase)*ones(nCase, 1);
% Concrete properties
fcMPA = fcSmp;
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

covEp = std(ep_fu) / mean(ep_fu);
ep_fe( ep_fe>0.004 ) = normrnd(0.004, 0.004*covEp, sum(ep_fe>0.004), 1);

%% Wrapping

covEpFe = F_FRP_COV;
ep_fe(isW) = normrnd(0.004, 0.004*covEpFe, sum(isW), 1);
isWandOver = (isW) & (ep_fe > 0.75*ep_fu);
ep_fe( isWandOver ) = 0.75*ep_fu(isWandOver);
ep_fe( ep_fe<0 ) = 0;


%% FRP contribution to shear resistance Eq. (11-3)

f_feMPA = ep_fe .* EFrpMPA;
shearFrpKN = 2*tFrpMM .* widthFrpMM .* f_feMPA .*...
        (sin(betaRAD)+cos(betaRAD)) .* dfvMM ./ sFrpMM /1e3;
shearFrpKN( shearFrpKN<0 ) = 0;

% ACI 318-11: Steel Contribution
nLarger = sum( fsMPA>fsLimit );
covLarger = FS_COV;
fsMPA( fsMPA>fsLimit ) = normrnd(fsLimit, fsLimit*covLarger, nLarger,1);
% fsMPA( fsMPA>fsLimit ) = fsLimit;

areaVMM2 = 2*pi/4*sdMM.^2;
shearSteelKN = areaVMM2 .* fsMPA .* dBeamMM ./ ssMM /1e3;
shearSteelKN( ssMM == 0 ) = 0;
shearSteelKN( shearSteelKN<0 ) = 0;

% contribution of concrete

areaV_Min = 0.062*sqrt(fcMPA) .* bBeamMM .* ssMM ./ fsMPA;
tmp = (0.35*bBeamMM .* sdMM) ./ fsMPA;
areaV_Min( areaV_Min < tmp ) = tmp(areaV_Min < tmp); 

sqrt_fc = sqrt(fcMPA);
isAdjust = (areaVMM2<areaV_Min) & (sqrt_fc>sqrtFcLimit);
nAdjust = sum( isAdjust );
covAdjust = std(sqrt_fc)/mean(sqrt_fc);
sqrt_fc( isAdjust ) = normrnd( sqrtFcLimit, sqrtFcLimit*covAdjust, nAdjust, 1);

shearConcreteKN = 0.17*sqrt_fc .* bBeamMM .* dBeamMM / 1e3;
shearConcreteKN(shearConcreteKN<0) = 0;

% over reinforcement check

shearReinforceKN = shearSteelKN + shearFrpKN;

isOverReinforce = (shearReinforceKN > 0.66*sqrt_fc.*bBeamMM .* dBeamMM / 1e3);
shearReinforceKN( isOverReinforce ) = 0.66*sqrt_fc(isOverReinforce).*...
                       bBeamMM(isOverReinforce).*dBeamMM(isOverReinforce)/1e3;
                   
shearTotalKN = shearConcreteKN + shearReinforceKN;

return
end