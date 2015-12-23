function [shearTotalKN, isOverReinforce,...
         yield_warning, shearReinforceKN] = shear_total_fib(FLAG, factorFrp)
% Determine the contribution of FRP, New HK guidelines
% Considering the steel-FRP interaction

nFactorFrp = length(factorFrp);
if nFactorFrp == 3
    gammaFrp = factorFrp(1);
    gammaBond = factorFrp(2);
    psi_f = factorFrp(3);
elseif nFactorFrp == 2
    gammaFrp = factorFrp(1);
    gammaBond = factorFrp(2);
    psi_f = 1;
elseif nFactorFrp == 1
    gammaFrp = factorFrp;
    gammaBond = factorFrp;
    psi_f = 1;
end

% important constants

TOP_DUMMY_HEIGHT = 0.1; % ratio of top dummy height to total depth
N_STEEL_LEG = 2; % number of steel legs
N_FRP_LEG = 2; % number of FRP legs
REL_TOL = 1e-6;
ES_MPA = 206e3;

ext_data = load('tmpdata.mat');

switch FLAG
    case 'MODEL_ERROR'
        nCase = length(ext_data.FC_TEST_ARRAY_MPA);
        % constant partial safety factors
        gammaConcrete = 1.00;
        gammaSteel = 1.00;
        gammaEfrp = 1.00;
        % geometrical properties
        hBeamMM = ext_data.H_TEST_ARRAY_MM;
        dBeamMM = ext_data.D_TEST_ARRAY_MM;
        bBeamMM = ext_data.B_TEST_ARRAY_MM;
        dFrpMM = ext_data.DFRP_TEST_ARRAY_MM;
        dFrpTopMM = ext_data.DFRP_TOP_TEST_ARRAY_MM;
        s2d = ext_data.S2D_TEST_ARRAY;
        % Concrete mean properties
        fcMPA = ext_data.FC_TEST_ARRAY_MPA;
        fckMPA = fcMPA - 8;
        fctMPA = 0.3 .* fckMPA.^(2/3.0);
        fcm = fcMPA;
        fctm = fctMPA;
        % steel mean properties
        fsMPA = ext_data.FS_TEST_ARRAY_MPA;
        sdMM = ext_data.SD_TEST_ARRAY_MM;
        ssMM = ext_data.SS_TEST_ARRAY_MM;
        % FRP mean properties
        frpForm = ext_data.FRP_FORM_TEST_ARRAY;
        fFrpMPA = ext_data.F_FRP_TEST_ARRAY_MPA;
        EFrpMPA = ext_data.E_FRP_TEST_ARRAY_MPA;
        betaFrpDEG = ext_data.BETA_TEST_ARRAY_DEG;
        tFrpMM = ext_data.T_FRP_TEST_ARRAY_MM;
        widthFrpMM = ext_data.W_FRP_TEST_ARRAY_MM;
        sFrpMM = ext_data.S_FRP_TEST_ARRAY_MM;
        % ratio of in-situ to cylinder strength
        roInsitu = 1.0;  
    case 'DESIGN_VALUE'
        nCase = length(ext_data.FC_DESIGN_ARRAY_MPA);
        % constant partial safety factors
        gammaConcrete = 1.50;
        gammaSteel = 1.15;  
        gammaEfrp = 1.00;
        % geometrical properties
        hBeamMM = ext_data.H_DESIGN_ARRAY_MM;
        dBeamMM = ext_data.D_DESIGN_ARRAY_MM;
        bBeamMM = ext_data.B_DESIGN_ARRAY_MM;
        dFrpMM = ext_data.DFRP_DESIGN_ARRAY_MM;
        dFrpTopMM = ext_data.DFRP_TOP_DESIGN_ARRAY_MM; 
        s2d = ext_data.S2D_DESIGN_ARRAY;
        % Concrete characteristic properties
        fcMPA = ext_data.FC_DESIGN_ARRAY_MPA;
        fckMPA = fcMPA;
        fcm = fckMPA+8;
        fctm = 0.3 .* fckMPA.^(2/3.0);
        % steel characteristic properties
        fsMPA = ext_data.FS_DESIGN_ARRAY_MPA;
        sdMM = ext_data.SD_DESIGN_ARRAY_MM;
        ssMM = ext_data.SS_DESIGN_ARRAY_MM;
        % FRP characteristic properties
        frpForm = ext_data.FRP_FORM_DESIGN_ARRAY;
        fFrpMPA = ext_data.F_FRP_DESIGN_ARRAY_MPA;
        EFrpMPA = ext_data.E_FRP_DESIGN_ARRAY_MPA;
        betaFrpDEG = ext_data.BETA_DESIGN_ARRAY_DEG;
        tFrpMM = ext_data.T_FRP_DESIGN_ARRAY_MM;
        widthFrpMM = ext_data.W_FRP_DESIGN_ARRAY_MM;
        sFrpMM = ext_data.S_FRP_DESIGN_ARRAY_MM;
        % ratio of in-situ to cylinder strength
        roInsitu = 0.85;         
    otherwise
end

%% general parameters
yield_warning = zeros(nCase, 1);
shearTotalKN = zeros(nCase, 1);
mufc = (30./fckMPA).^(1/3); mufc(mufc>1) = 1;
zt = dFrpTopMM;
zb = (dBeamMM-(hBeamMM-dFrpMM))-TOP_DUMMY_HEIGHT*dBeamMM;
hFrpEff = zb-zt;
theta_final = zeros(nCase,1);
beta = betaFrpDEG/180*pi;

%% FRP contribution (independent from theta, 2010 Code, Lb is based on Chen and Teng)
isSide = (frpForm == 1);
isU = (frpForm == 2);
isW = (frpForm == 3);

roFrp = (2*tFrpMM./bBeamMM).*(widthFrpMM./sFrpMM);

switch FLAG
    case 'MODEL_ERROR'
        k = 1.0;
    case 'DESIGN_VALUE'
        k = 0.8;
    otherwise
end
edeb = psi_f*k*0.65*(fcm.^(2/3)./(1e-3*EFrpMPA.*roFrp)).^(0.65)*1e-3/gammaBond;
efu = (fFrpMPA./gammaFrp)./(EFrpMPA./gammaEfrp);
erup = psi_f*k*0.17.*(fcm.^(2/3)./(EFrpMPA.*1e-3.*roFrp)).^(0.30).*efu;
efe = min(edeb, erup);
efe(isW) = erup(isW);

Vf = efe.*(EFrpMPA./gammaEfrp).*roFrp.*bBeamMM.*hFrpEff.*(sin(beta)+cos(beta));
Vf(Vf<0) =0;

%% Without steel stirrups or less than minimum shear reinforcement
indx01 = ssMM==0 | sdMM==0;
sqrtfck = sqrt(fckMPA); sqrtfck(sqrtfck>8) = 8;
ro_s = ones(nCase, 1);
Asw = N_STEEL_LEG * pi* sdMM.^2/4;
ro_s(~indx01) = Asw(~indx01)./(bBeamMM(~indx01).*ssMM(~indx01));
indx02 = ro_s<0.08*sqrtfck./fsMPA;
indx0 = indx01 | indx02;
kv = 180./(1000+1.25*(0.9*dBeamMM));
Vcd = kv.*sqrt(fcm)./gammaConcrete.*(0.9*dBeamMM).*bBeamMM;
shearConc = Vcd;
shearTotalKN(indx0) = 1e-3*(shearConc(indx0)+Vf(indx0));
kep = 0.55;     % Eq. 7.3-37
Vrdmax0kN = kep.*mufc.* (fcMPA/gammaConcrete).* bBeamMM .* (0.9*dBeamMM) ./ 2*1e-3;
VrdkN = min(shearTotalKN, Vrdmax0kN);
shearTotalKN(indx0) = VrdkN(indx0);

%% FRP failure and steel yielding with theta=45\deg
mindeg = 45;
theta = mindeg/180*pi;

Vsd = Asw./ssMM.*(0.9*dBeamMM).*(fsMPA/gammaSteel).*cot(theta);
Vrdsf = Vcd+Vsd+Vf;
VrdkN = min(Vrdsf*1e-3, Vrdmax0kN);
shearTotalKN(~indx0) = VrdkN(~indx0);

isOverReinforce = shearTotalKN==Vrdmax0kN;
yield_warning = zeros(nCase,1);
shearReinforceKN = Vf*1e-3;

return
end