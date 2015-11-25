function shearTotalKN = shear_total_TR_MC(iDesignCase, hSmp, betaSmp, fFrpSmp, fcSmp, fctSmp, sqrtFctSmp, fsSmp)
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

% % limit values
% fcuLimit1 = 50*FC_BIAS;
% fcuLimit2 = 80*FC_BIAS;
% fcuBrittleLimit1 = 40*FC_BIAS;
% fcuBrittleLimit2 = 80*FC_BIAS;
% % ratio of in-situ to cylinder strength
% roInsitu = 0.88;

%% general parameters
yield_warning = zeros(nCase, 1);
shearTotalKN = zeros(nCase, 1);
alpha_cw = ones(nCase,1); % no axial force
v1 = 0.6*(1-fckMPA/250);
zt = dFrpTopMM;
zb = (dBeamMM-(hBeamMM-dFrpMM))-TOP_DUMMY_HEIGHT*dBeamMM;
hFrpEff = zb-zt;
beta = betaFrpDEG/180*pi;

%% FRP contribution (independent from theta)
isSide = (frpForm == 1);
isU = (frpForm == 2);
isW = (frpForm == 3);

ns = zeros(nCase, 1); ns(isU) = 1; ns(isSide) = 2;
Efd = EFrpMPA/gammaEfrp;
ltmax = 0.7*sqrt( Efd.*tFrpMM./fctm); ltmax(fctm==0) = hFrpEff(fctm==0)./sin(beta(fctm==0));
epf1 = 0.5*(fFrpMPA/gammaFrp)./Efd;
epf2 = 0.5*sqrt(fctm./(Efd.*tFrpMM)) / (gammaBond);
epf3 = epf1;
covEp = std(epf1) / mean(epf1);
epf3( epf3>0.004 ) = normrnd(0.004, 0.004*covEp, sum(epf3>0.004), 1);
epf = min([epf1'; epf2'; epf3'])';
Afw = N_FRP_LEG .* widthFrpMM .* tFrpMM;  
Vf = Afw./sFrpMM.*(hFrpEff-ns/3.*ltmax.*cos(beta)).*Efd.*epf.*(sin(beta)+cos(beta));

%% Without steel stirrups
indx0 = ssMM==0 | sdMM==0;
%Vrd = (Crd*k*(100*rhol*fck).^(1/3)+k1*scp)*bw*d;
k = 1+sqrt(200./dBeamMM);
shearConc = (0.18/gammaConcrete).*k.*(fckMPA).^(1/3).*bBeamMM.*dBeamMM ;
shearConcMin = 0.035*k.^(3/2).*sqrt(fckMPA).*bBeamMM.*dBeamMM;
shearConc(shearConc<shearConcMin) = shearConcMin(shearConc<shearConcMin);
shearTotalMinKN = 1e-3*(shearConc+Vf);
shearTotalKN(indx0) = 1e-3*(shearConc(indx0)+Vf(indx0));
Vrdmax0kN = alpha_cw .* bBeamMM .* (0.9*dBeamMM) .* v1 .* (fcMPA/gammaConcrete) ./ 2*1e-3;
indx0c = (indx0 & shearTotalKN>Vrdmax0kN);
shearTotalKN(indx0c) = Vrdmax0kN(indx0c);

shearTotalMinKN(shearTotalMinKN>Vrdmax0kN) = shearTotalMinKN(shearTotalMinKN>Vrdmax0kN);

%% FRP failure and steel yielding with theta=22\deg
theta = 22/180*pi;
% local function
function Vrdsf = shear_resistance(theta)
    Asw = N_STEEL_LEG * pi* sdMM.^2/4;
    Vrd = Asw./ssMM.*(0.9*dBeamMM).*(fsMPA/gammaSteel).*cot(theta);
    Vrdsf = Vrd+Vf;
end
Vrdsf22 = shear_resistance(theta);
Vrdmax22 = alpha_cw .* bBeamMM .* (0.9*dBeamMM) .* v1 .* (fcMPA/gammaConcrete) ./ (cot(theta)+tan(theta));

indx1 = (Vrdmax22>Vrdsf22) & (~indx0);
shearTotalKN(indx1) = Vrdsf22(indx1) *1e-3;

yield_warning(indx1&(fsMPA./ES_MPA>epf)) = 1;

%% concrete crushing with 45\deg strut
theta0 = 45/180*pi;
Vrdsf45 = shear_resistance(theta0);
Vrdmax45 = alpha_cw .* bBeamMM .* (0.9*dBeamMM) .* v1 .* (fcMPA/gammaConcrete) ./ (cot(theta0)+tan(theta0));
indx2 = (Vrdmax45<Vrdsf45) & (~indx0);
shearTotalKN(indx2) = Vrdmax45(indx2) *1e-3;

%% concrete crushing with other deg
indx3 = (~indx1)&(~indx2)&(~indx0);
theta1 = 22/180*pi*ones(nCase,1); theta1 = theta1(indx3);
theta2 = 45/180*pi*ones(nCase,1); theta2 = theta2(indx3);
delta = 1; delta1 = 0.5; delta2 = 0.5;
conv_tol = sum(delta) / (sum(abs(delta1))+sum(abs(delta2)));
if ~isempty(indx3)
    max_shear = @(theta) alpha_cw(indx3) .* bBeamMM(indx3) .* (0.9*dBeamMM(indx3)) .* v1(indx3) .* (fcMPA(indx3)/gammaConcrete) ./ (cot(theta)+tan(theta));

%     while norm(theta1-theta2)>THETA_TOL*sqrt(sum(indx3))
    while conv_tol>REL_TOL
        thetatmp = (22/180*pi)*ones(nCase,1); thetatmp(indx3) = theta1;
        Vtmp1 = shear_resistance(thetatmp); Vtmp1 = Vtmp1(indx3);
        thetatmp = (22/180*pi)*ones(nCase,1); thetatmp(indx3) = theta2;
        Vtmp2 = shear_resistance(thetatmp); Vtmp2 = Vtmp2(indx3);    
        delta1 = max_shear(theta1)-Vtmp1;
        delta2 = max_shear(theta2)-Vtmp2;
        theta = theta1 - (theta2-theta1).*delta1./(delta2-delta1);

        thetatmp = (22/180*pi)*ones(nCase,1); thetatmp(indx3) = theta;
        Vtmp = shear_resistance(thetatmp); Vtmp = Vtmp(indx3);  
        delta = max_shear(theta)-Vtmp;

        theta1(delta>0) = theta1(delta>0); theta2(delta>0)=theta(delta>0);
        theta1(delta<0) = theta(delta<0); theta2(delta<0)=theta2(delta<0);
        theta1(delta==0) = theta(delta==0); theta2(delta==0)=theta(delta==0);
        theta1(isnan(delta)) = theta1(isnan(delta)); theta2(isnan(delta))=theta2(isnan(delta));
        
        conv_tol = sum(delta) ./ (sum(abs(delta1))+sum(abs(delta2)));
    end
    thetatmp = (22/180*pi)*ones(nCase,1); thetatmp(indx3) = (theta1+theta2)/2;
    Vtmp = shear_resistance(thetatmp); Vtmp = Vtmp(indx3); 
    shearTotalKN(indx3) = Vtmp*1e-3;
end

shearTotalKN(shearTotalKN<shearTotalMinKN) = shearTotalMinKN(shearTotalKN<shearTotalMinKN);
isOverReinforce = indx0c | indx2 | indx3;
shearReinforceKN = Vf*1e-3;

return
end