function shearTotalKN = shear_total_fib_MC(iDesignCase, hSmp, betaSmp, fFrpSmp, fcSmp, fctSmp, sqrtFctSmp, fsSmp)
% Determine the contribution of FRP, New HK guidelines
% Considering the steel-FRP interaction

% important constants
TOP_DUMMY_HEIGHT = 0.1; % ratio of top dummy height to total depth
N_STEEL_LEG = 2; % number of steel legs
N_FRP_LEG = 2; % number of FRP legs
REL_TOL = 1e-3;
ES_MPA = 206e3;

ext_data = load('tmpdata.mat');

% read variables
nCase = length(hSmp);
% constant partial safety factors
gammaConcrete = 1.00;
gammaSteel = 1.00;
gammaFrp = 1.00;
gammaEfrp = 1.00;
gammaBond = 1.00;
% geometrical properties
hBeamMM = hSmp;
dBeamMM = hBeamMM - ext_data.AS_DESIGN_ARRAY_MM(iDesignCase);
bBeamMM = hSmp/ext_data.H2B_DESIGN_ARRAY(iDesignCase);
dFrpMM = hSmp;
dFrpTopMM = ext_data.DFRP_TOP_DESIGN_ARRAY_MM(iDesignCase)*ones(nCase, 1);
s2d = ext_data.S2D_DESIGN_ARRAY(iDesignCase)*ones(nCase, 1);
% Concrete properties
fcMPA = fcSmp;
fckMPA = fcMPA-8; fckMPA(fckMPA<0)=0;
fctMPA = 0.3 .* fckMPA.^(2/3.0);
fcm = fcMPA;
fctm = fctMPA;
% steel properties
fsMPA = fsSmp;
sdMM = ext_data.SD_DESIGN_ARRAY_MM(iDesignCase)*ones(nCase, 1);
ssMM = ext_data.SS_DESIGN_ARRAY_MM(iDesignCase)*ones(nCase, 1);
% FRP properties
frpForm = ext_data.FRP_FORM_DESIGN_ARRAY(iDesignCase)*ones(nCase, 1);
fFrpMPA = fFrpSmp;
EFrpMPA = ext_data.E_FRP_DESIGN_ARRAY_MPA(iDesignCase)*ones(nCase, 1);
betaFrpDEG = betaSmp;
tFrpMM = ext_data.T_FRP_DESIGN_ARRAY_MM(iDesignCase)*ones(nCase, 1);
widthFrpMM = ext_data.W_FRP_DESIGN_ARRAY_MM(iDesignCase)*ones(nCase, 1);
sFrpMM = ext_data.S_FRP_DESIGN_ARRAY_MM(iDesignCase)*ones(nCase, 1);

%% general parameters
yield_warning = zeros(nCase, 1);
shearTotalKN = zeros(nCase, 1);
mufc = (30./fckMPA).^(1/3); mufc(mufc>1) = 1;
zt = dFrpTopMM;
zb = (dBeamMM-(hBeamMM-dFrpMM))-TOP_DUMMY_HEIGHT*dBeamMM;
hFrpEff = zb-zt;
beta = betaFrpDEG/180*pi;

%% FRP contribution (independent from theta)
isSide = (frpForm == 1);
isU = (frpForm == 2);
isW = (frpForm == 3);

roFrp = (2*tFrpMM./bBeamMM).*(widthFrpMM./sFrpMM);
k=1.0;
edeb = k*0.65*(fcm.^(2/3)./(1e-3*EFrpMPA.*roFrp)).^(0.65)*1e-3/gammaBond;
efu = (fFrpMPA./gammaFrp)./(EFrpMPA./gammaEfrp);
erup = k*0.17.*(fcm.^(2/3)./(EFrpMPA.*1e-3.*roFrp)).^(0.30).*efu;
efe = min(edeb, erup);
efe(isW) = erup(isW);

Vf = efe.*(EFrpMPA./gammaEfrp).*roFrp.*bBeamMM.*hFrpEff.*(sin(beta)+cos(beta));
Vf(Vf<0) =0;

%% Without steel stirrups
indx01 = ssMM==0 | sdMM==0;
sqrtfck = sqrt(fckMPA); sqrtfck(sqrtfck>8) = 8;
ro_s = ones(nCase, 1);
Asw = N_STEEL_LEG * pi* sdMM.^2/4;
ro_s(~indx01) = Asw(~indx01)./(bBeamMM(~indx01).*ssMM(~indx01));
indx02 = ro_s<0.08*sqrtfck./fsMPA;
indx0 = indx01 | indx02;
kv = 180./(1000+1.25*(0.9*dBeamMM));
Vcd = kv.*sqrt(fcMPA)./gammaConcrete.*(0.9*dBeamMM).*bBeamMM;
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

return
end