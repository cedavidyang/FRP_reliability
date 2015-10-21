function [momentTotalKNM, failMode, isSteelYielding]= flexure_total_GB(FLAG, factorFrp)
% Compute maximum moment using ACI 440.2R-08
% Input: FLAG --- type of the analysis
%        factorFrp -- partial safety factor or resistance reduction factor
% Output: momentTotalKNM --- predicted maximum moment
%         failMode --- fail mode: 1 for IC debonding;
%                               2 for FRP rupture
%                               3 for Concrete crush;
%         isSteelYielding ---  1 (tensile steel yielding before
%                                 concrete crush or FRP failure)
%                              0 (tensile steel does not yield)

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

STRAIN_ULTI_CONCRETE = 0.0033;
STRAIN_INITIAL_FRP = 0;
MAX_INITIAL_GUESS = 10; % maximum number of initial guess for the depth of stress block
RO_CYLINDE_2_CUBE = 0.8; % strength ratio of cylinders to cubes
CONVERGE_CRITERION = 1E-5;

load tmpdata

%% data extracted from database
switch FLAG
    case{'MODEL_ERROR'}
        % constant partial safte factors
        gammaConcrete = 1.00;
        gammaSteel = 1.00;
        % geometric data
        bBeamMM = B_TEST_ARRAY_MM;
        hBeamMM = H_TEST_ARRAY_MM;
        bFlangeMM = BF_TEST_ARRAY_MM;
        tFlangeMM = TF_TEST_ARRAY_MM;
        dBeamMM = D_TEST_ARRAY_MM;
        dCmpMM = D_CMP_TEST_ARRAY_MM;
        aFrpMM = FRP_END_TEST_ARRAY_MM;
        shearMM = SHEAR_TEST_ARRAY_MM;
        % concrete properties
        fcMPA = FC_TEST_ARRAY_MPA;
        fcuMPA = fcMPA / RO_CYLINDE_2_CUBE;
        fctMPA = 0.5*sqrt(fcuMPA);
        % steel reinforcement
        EsMPA = ES_TEST_ARRAY_MPA;
        fsMPA = FS_TEST_ARRAY_MPA;
        areaSteelMM2 = AREA_STEEL_TEST_ARRAY_MM2;
        EsCmpMPA = ES_CMP_TEST_ARRAY_MPA;
        fsCmpMPA = FS_CMP_TEST_ARRAY_MPA;
        areaSteelCmpMM2 = AREA_STEEL_CMP_TEST_ARRAY_MM2;
        % FRP strengthening
        EfrpMPA = E_FRP_TEST_ARRAY_MPA;
        fFrpMPA = F_FRP_TEST_ARRAY_MPA;
        tFrpMM = T_FRP_TEST_ARRAY_MM;
        bFrpMM = B_FRP_TEST_ARRAY_MM;
        areaFrpMM2 = bFrpMM .* tFrpMM;
        isAnchor = ANCHOR_TEST_ARRAY;
        % limit values
        fcuLimit1 = 50*FC_BIAS;
        fcuLimit2 = 80*FC_BIAS;
        fcuBrittleLimit1 = 40*FC_BIAS;
        fcuBrittleLimit2 = 80*FC_BIAS;        
        % ratio of in-situ to cylinder strength
        roInsitu = 1.0;          
    case{'DESIGN_VALUE'}
        % constant partial safte factors
        gammaConcrete = 1.40;
        gammaSteel = 1.10;
        % geometric data
        bBeamMM = B_DESIGN_ARRAY_MM;
        hBeamMM = H_DESIGN_ARRAY_MM;
        bFlangeMM = BF_DESIGN_ARRAY_MM;
        tFlangeMM = TF_DESIGN_ARRAY_MM;
        dBeamMM = D_DESIGN_ARRAY_MM;
        dCmpMM = D_CMP_DESIGN_ARRAY_MM;
        aFrpMM = FRP_END_DESIGN_ARRAY_MM;
        shearMM = SHEAR_DESIGN_ARRAY_MM;
        % concrete properties
        fcMPA = FC_DESIGN_ARRAY_MPA;
        fcuMPA = fcMPA / RO_CYLINDE_2_CUBE;
        fctMPA = 0.395*fcuMPA.^0.55;
        % steel reinforcement
        EsMPA = ES_DESIGN_ARRAY_MPA;
        fsMPA = FS_DESIGN_ARRAY_MPA;
        areaSteelMM2 = AREA_STEEL_DESIGN_ARRAY_MM2;
        EsCmpMPA = ES_CMP_DESIGN_ARRAY_MPA;
        fsCmpMPA = FS_CMP_DESIGN_ARRAY_MPA;
        areaSteelCmpMM2 = AREA_STEEL_CMP_DESIGN_ARRAY_MM2;
        % FRP strengthening
        EfrpMPA = E_FRP_DESIGN_ARRAY_MPA;
        fFrpMPA = F_FRP_DESIGN_ARRAY_MPA;
        tFrpMM = T_FRP_DESIGN_ARRAY_MM;
        bFrpMM = B_FRP_DESIGN_ARRAY_MM;
        areaFrpMM2 = bFrpMM .* tFrpMM;
        isAnchor = ANCHOR_DESIGN_ARRAY;
        % limit values
        fcuLimit1 = 50;
        fcuLimit2 = 80;  
        fcuBrittleLimit1 = 40;
        fcuBrittleLimit2 = 80;        
        % ratio of in-situ to cylinder strength
        roInsitu = 0.88;          
    otherwise
end

%% GB 50608 model
%% determine debonding strain of FRP
nCase = length(bBeamMM);
% about concrete tensile strengths
% alphaBrittle = ones(nCase, 1);
% isHighLimit1 = (fcuMPA > fcuBrittleLimit1) & (fcuMPA < fcuBrittleLimit2);
% isHighLimit2 = (fcuMPA >= fcuBrittleLimit2);
% alphaBrittle(isHighLimit1) = (1.0-0.87) / (fcuBrittleLimit1-fcuBrittleLimit2) *...
%                              (fcuMPA(isHighLimit1)-fcuBrittleLimit2) + 0.87;
% alphaBrittle(isHighLimit2) = 0.87;
% fctMPA = alphaBrittle .* fctMPA;
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
% from eqution: ftk = 0.88*0.395*fcuk.^0.55.*(1-1.645*0.2).^0.45.*0.8;
% 0.88, as ratio of in-situ concrete strength to cylinder strength, should
% be eliminated

% fcMPA = roInsitu * fcMPA;
% fctMPA = roInsitu * fctMPA;
fctMPA = roInsitu .* alphaBrittle .* fctMPA;
fcMPA = roInsitu .* prism2Cube .* alphaBrittle .* fcuMPA;

failMode = ones(nCase , 1 );
betaW = sqrt( ( 2.25 - bFrpMM./bBeamMM ) ./ (1.25 + bFrpMM./bBeamMM) );
Ld = shearMM - aFrpMM;
eDebond = (1.1./sqrt(EfrpMPA.*tFrpMM) - 0.2./Ld) .* betaW .* (fctMPA/gammaConcrete) ./ gammaBond;
eFrpUlti = (fFrpMPA/gammaFrp) ./ EfrpMPA;
isRupture = eDebond > eFrpUlti;
failMode( isRupture|isAnchor ) = 2;
eFrpCritical = zeros(nCase, 1);
eFrpCritical(isRupture|isAnchor) = eFrpUlti( isRupture|isAnchor );
eFrpCritical((~isRupture)&(~isAnchor)) = eDebond( (~isRupture)&(~isAnchor) );
% if strcmp(FLAG, 'DESIGN_VALUE')
%     eFrpCritical( (isRupture|isAnchor)&(eFrpCritical>0.01) ) = 0.01;
% end
%% determine c
beta1 = 0.8*ones(nCase, 1);
beta1( fcuMPA>fcuLimit2 ) = 0.74;
beta1( fcuMPA<=fcuLimit2 & fcuMPA>=fcuLimit1 ) = ...
    (0.74 - 0.8)./(fcuLimit2-fcuLimit1).*(fcuMPA( fcuMPA<=fcuLimit2 & fcuMPA>=fcuLimit1 ) - fcuLimit1) + 0.8;

eFrpIni = STRAIN_INITIAL_FRP;
eCo = 0.002 + 0.5*(fcuMPA-fcuLimit1)*1e-5;
eConcUlti = STRAIN_ULTI_CONCRETE - (fcuMPA-fcuLimit1)*1e-5;
eConcUlti( eConcUlti>STRAIN_ULTI_CONCRETE ) = STRAIN_ULTI_CONCRETE;
isSteelYielding = ones( nCase, 1 ); 

c = zeros(nCase, 1);
fval = zeros(nCase, 1);
exitflag = zeros(nCase, 1);

for iCase = 1:nCase
    cIni = linspace(dCmpMM(iCase), dBeamMM(iCase), MAX_INITIAL_GUESS);
    for iC = 1:MAX_INITIAL_GUESS
        [c(iCase), fval(iCase), exitflag(iCase)] = fsolve( @(c) force_balance(beta1(iCase), eCo(iCase), EsMPA(iCase), EsCmpMPA(iCase), EfrpMPA(iCase),...
                                              eFrpIni, eConcUlti(iCase), eFrpCritical(iCase), bBeamMM(iCase), bFlangeMM(iCase), tFlangeMM(iCase), dBeamMM(iCase), dCmpMM(iCase),...
                                              hBeamMM(iCase), c, (fsMPA(iCase)/gammaSteel), (fsCmpMPA(iCase)/gammaSteel), (fcMPA(iCase)/gammaConcrete),...
                                              areaSteelMM2(iCase), areaFrpMM2(iCase), areaSteelCmpMM2(iCase)),...
                                   cIni(iC), optimset('Display','off'));
        if (exitflag(iCase) == 1 || abs(fval(iCase))<=CONVERGE_CRITERION) && (c(iCase)<=hBeamMM(iCase))
            break
        end
    end
    if iC > MAX_INITIAL_GUESS
        fprintf('Failed convergence at %d, exitflag = %d, c=%.2f\n', iCase, exitflag(iCase), c(iCase))
    end
%% other approaches
%     if dCmpMM(iCase) > 0.5 * beta1(iCase) * c(iCase)
%         % Approach 1: neglect compressive steel and assume c = 2d_c
% %        areaSteelCmpMM2(iCase) = 0;
% %        c = 2*dCmpMM(iCase)./beta1(iCase);
%         % Approach 2: neglect compressive steel and determine c again
% %         for iC = 1:10
% %             % assuming no compressive steel
% %             [c, fval, exitflag] = fsolve( @(c) force_balance(alpha1(i), beta1(i), EsMPA(i), EsCmpMPA(i), EfrpMPA(i),...
% %                 e_bi, eConcUlti(i), eFrpCritical(i), bBeamMM(i), bFlangeMM(i), tFlangeMM(i), dBeamMM(i), dCmpMM(i),...
% %                 hBeamMM(i), c, fsMPA(i), fsCmpMPA(i), fcMPA(i), areaSteelMM2(i), areaFrpMM2(i), 0),...
% %                 cIni(iC), optimset('Display','off'));
% %             x(i, 1) = c; x(i, 2) = exitflag;
% %             if exitflag == 1 || abs(fval)<=1e-5
% %                 break
% %             end
% %         end
%        % Approach 3: consider the effects of compressive steel --- Do nothing
% %       disp('Compressive rebar is too deep')
%     end
%     % uncomment if approach 2 is used
% %     if iC == 11
% %         fprintf('Failed convergence at %d, exitflag = %d, c=%.2f\n', i, exitflag, c)
% %     end 
%     x(iCase, 1) = c; x(iCase, 2) = exitflag;
end
[ffe, eFrpEff, failMode] = FRP_stress(EfrpMPA, eConcUlti, hBeamMM, c, eFrpIni, eFrpCritical, failMode);
[fs, es, isSteelYielding] = tension_steel(EsMPA, eFrpEff, eFrpIni, dBeamMM, hBeamMM, c, (fsMPA/gammaSteel), isSteelYielding);
fsc = comp_steel(eCo, eConcUlti, EsCmpMPA, eFrpEff, eFrpIni, dCmpMM, hBeamMM, c, (fsCmpMPA/gammaSteel));
switch FLAG
    case {'MODEL_ERROR'}
        phi = 1;
        momentTotalKNM = 1e-6 * phi .* ( areaSteelMM2.*fs.*(dBeamMM-0.5*beta1.*c) + areaFrpMM2.*ffe.*(hBeamMM-0.5*beta1.*c) + (0.5*beta1.*c-dCmpMM).*areaSteelCmpMM2.*fsc );    
    case {'DESIGN_VALUE'}
        phi = 1;
%         psi_f = psi_f*ones(nCase,1);
%         psi_f((~isRupture)&(~isAnchor)) = 0.90;
%         psi_f(isRupture|isAnchor) = 0.50;
        momentTotalKNM = 1e-6 * phi .* ( areaSteelMM2.*fs.*(dBeamMM-0.5*beta1.*c) + psi_f.*areaFrpMM2.*ffe.*(hBeamMM-0.5*beta1.*c) + (0.5*beta1.*c-dCmpMM).*areaSteelCmpMM2.*fsc );  
    otherwise
end

return
end

function [r,f] = force_balance(beta1, e_co, Es, Es_c, Efrp, e_ini, e_cu, e_fd, b, bf, tf, d, d_c, df, c, fyd, fy_cd, fcd, As, Afrp, As_c)
% balance residual
nCase = length(Es);
tmp = zeros(nCase,1);

cr = c.*beta1; %relative depth of neutral axis

[ffe, e_fe, ~] = FRP_stress(Efrp, e_cu, df, c, e_ini, e_fd, tmp);
[fs, ~] = tension_steel(Es, e_fe, e_ini, d, df, c, fyd, tmp);

% e_cf = (e_fe + e_ini ) .* c ./ (hBeamMM-c); 
% w = ones(nCase, 1);
% w( e_cf<e_cu ) = 0.5 + 0.5 * e_fe(e_cf<e_cu) ./ e_fd(e_cf<e_cu);
w = 0.5 + 0.5 * e_fe ./ e_fd;

fsc = comp_steel(e_co, e_cu, Es_c, e_fe, e_ini, d_c, df, c, fy_cd);
Ac = zeros(nCase, 1);
Ac(bf == b) = b(bf == b) .*cr(bf == b);
Ac( cr<=tf & bf~=b ) = bf(cr<=tf & bf~=b) .* cr( cr<=tf & bf~=b ); % Type I T section
Ac( cr>tf & bf~=b ) = b( cr>tf & bf~=b ) .* cr( cr>tf & bf~=b ) + ( bf( cr>tf & bf~=b )-b( cr>tf & bf~=b )) .*tf( cr>tf & bf~=b ); % Type II T section
r = As.*fs + Afrp.*ffe - As_c.*fsc - (w.*fcd).*Ac;
f = 0.5 * ( (As.*fs + Afrp.*ffe) + (As_c.*fsc + (w.*fcd).*Ac) );

return
end

function [ss e_s isSteelYielding] = tension_steel(Es, e_fe, e_ini, d, df, c, fyd, isSteelYielding)
% determine tensile steel stress
e_s = (e_fe+e_ini).*(d-c)./(df-c);
ss = Es .* e_s;
isSteelYielding( ss<fyd ) = 0;
ss( ss>fyd ) = fyd( ss>fyd );

return
end

function ss = comp_steel(e_co, e_cu, Es_c, e_fe, e_ini, d_c, df, c, fy_cd)
% determine compressive steel stress

e_sc = (c-d_c) ./ (df-c) .* (e_fe+e_ini);
ss = Es_c .* e_sc;
ss( ss>fy_cd ) = fy_cd( ss>fy_cd );
ss( ss<-fy_cd ) = -fy_cd( ss<-fy_cd  );

return
end

function [sf e_fe failMode] = FRP_stress(Efrp, e_cu, df, c, e_ini, e_fd, failMode)
% determine FRP stress
e_fe = e_cu.*(df-c)./c - e_ini;
failMode( e_fe<e_fd ) = 3;
e_fe( e_fe>e_fd ) = e_fd( e_fe>e_fd );
sf = Efrp .* e_fe;

return
end