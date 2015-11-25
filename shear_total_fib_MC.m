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

kbl = 2.0;
lEff = sqrt( EFrpMPA.*tFrpMM ./ (kbl.*fctm) ); lEff(fctm==0) = hFrpEff(fctm==0)./sin(beta(fctm==0));
lb = lEff;
lb(isSide) = hFrpEff(isSide) ./ 2 ./ sin(beta(isSide));
lb(isU) = hFrpEff(isU) ./ sin(beta(isU));

betaL = lb./lEff .* (2-lb./lEff); betaL(lEff==0) = 1.0;
betaL( lb>lEff ) = 1.0;

tmp = (2-widthFrpMM./(sFrpMM.*sin(beta)))./(1+widthFrpMM./(sFrpMM.*sin(beta)));
tmp( tmp<0 ) = 0;
kb = sqrt(tmp);
kb(kb<1) = 1.0;

kk = 0.17; kc = 1.0;
edeb = kc*(kk./gammaBond).*kb.*betaL.*sqrt(2*fcm.^(2/3.0)./(EFrpMPA.*tFrpMM));
roFrp = (2*tFrpMM./bBeamMM).*(widthFrpMM./sFrpMM);
efu = (fFrpMPA./gammaFrp)./(EFrpMPA./gammaEfrp);
erup = 0.8*0.17.*(fcm.^(2/3.0)./(EFrpMPA.*1e-3.*roFrp)).^(0.30).*efu;
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
shearConc = kv.*sqrt(fcMPA)./gammaConcrete.*(0.9*dBeamMM).*bBeamMM;
shearTotalMinKN = 1e-3*(shearConc+Vf);
shearTotalKN(indx0) = 1e-3*(shearConc(indx0)+Vf(indx0));
kep = 0.55;     % Eq. 7.3-37
Vrdmax0kN = kep.*mufc.* (fcMPA/gammaConcrete).* bBeamMM .* (0.9*dBeamMM) ./ 2*1e-3;
indx0c = (indx0 & shearTotalKN>Vrdmax0kN);
shearTotalKN(indx0c) = Vrdmax0kN(indx0c);
shearTotalMinKN(shearTotalMinKN>Vrdmax0kN) = shearTotalMinKN(shearTotalMinKN>Vrdmax0kN);

%% FRP failure and steel yielding with theta=22\deg
mindeg = 30;
theta = mindeg/180*pi;
% local function
function Vrdsf = shear_resistance(theta)
    Vrd = Asw./ssMM.*(0.9*dBeamMM).*(fsMPA/gammaSteel).*cot(theta);
    Vrdsf = Vrd+Vf;
end
Vrdsfmindeg = shear_resistance(theta);
Vrdmaxmindeg = kep.*mufc .* bBeamMM .* (0.9*dBeamMM) .* (fcMPA/gammaConcrete) .* (sin(theta)+cos(theta));

indx1 = (Vrdmaxmindeg>Vrdsfmindeg) & (~indx0);
shearTotalKN(indx1) = Vrdsfmindeg(indx1) *1e-3;

yield_warning(indx1&(fsMPA./ES_MPA>efe)) = 1;
theta_final(indx1) = theta;

%% concrete crushing with 45\deg strut
theta0 = 45/180*pi;
Vrdsf45 = shear_resistance(theta0);
Vrdmax45 = kep.*mufc.* (fcMPA/gammaConcrete).* bBeamMM .* (0.9*dBeamMM) ./ 2;
indx2 = (Vrdmax45<Vrdsf45) & (~indx0);
shearTotalKN(indx2) = Vrdmax45(indx2) *1e-3;
theta_final(indx2) = theta0;

%% concrete crushing with other deg
indx3 = (~indx1)&(~indx2)&(~indx0);
theta1 = mindeg/180*pi*ones(nCase,1); theta1 = theta1(indx3);
theta2 = 45/180*pi*ones(nCase,1); theta2 = theta2(indx3);
delta = 1; delta1 = 0.5; delta2 = 0.5;
conv_tol = sum(delta) / (sum(abs(delta1))+sum(abs(delta2)));
if ~isempty(indx3)
    max_shear = @(theta) kep*mufc(indx3) .* bBeamMM(indx3) .* (0.9*dBeamMM(indx3)) .* (fcMPA(indx3)/gammaConcrete) .* (sin(theta)+cos(theta));

%     while norm(theta1-theta2)>THETA_TOL*sqrt(sum(indx3))
    while conv_tol>REL_TOL
        thetatmp = (mindeg/180*pi)*ones(nCase,1); thetatmp(indx3) = theta1;
        Vtmp1 = shear_resistance(thetatmp); Vtmp1 = Vtmp1(indx3);
        thetatmp = (mindeg/180*pi)*ones(nCase,1); thetatmp(indx3) = theta2;
        Vtmp2 = shear_resistance(thetatmp); Vtmp2 = Vtmp2(indx3);    
        delta1 = max_shear(theta1)-Vtmp1;
        delta2 = max_shear(theta2)-Vtmp2;
        theta = theta1 - (theta2-theta1).*delta1./(delta2-delta1);

        thetatmp = (mindeg/180*pi)*ones(nCase,1); thetatmp(indx3) = theta;
        Vtmp = shear_resistance(thetatmp); Vtmp = Vtmp(indx3);  
        delta = max_shear(theta)-Vtmp;

        theta1(delta>0) = theta1(delta>0); theta2(delta>0)=theta(delta>0);
        theta1(delta<0) = theta(delta<0); theta2(delta<0)=theta2(delta<0);
        theta1(delta==0) = theta(delta==0); theta2(delta==0)=theta(delta==0);
        theta1(isnan(delta)) = theta1(isnan(delta)); theta2(isnan(delta))=theta2(isnan(delta));
        
        conv_tol = sum(delta) ./ (sum(abs(delta1))+sum(abs(delta2)));
    end
    thetatmp = (mindeg/180*pi)*ones(nCase,1); thetatmp(indx3) = (theta1+theta2)/2;
    Vtmp = shear_resistance(thetatmp); Vtmp = Vtmp(indx3); 
    shearTotalKN(indx3) = Vtmp*1e-3;
    theta_final(indx3) = thetatmp(indx3);
end

shearTotalKN(shearTotalKN<shearTotalMinKN) = shearTotalMinKN(shearTotalKN<shearTotalMinKN);
isOverReinforce = indx0c | indx2 | indx3;
shearReinforceKN = Vf*1e-3;

return
end