function [momentTotalKNM, failMode, isSteelYielding]= flexure_total_HK(FLAG, factorFrp)
% Compute maximum moment using Hong Kong guideline
% Input: FLAG --- type of the analysis
%        factorFrp -- partial safety factor or resistance reduction factor
% Output: momentTotalKNM --- predicted maximum moment
%         failMode --- fail mode: 1 for IC debonding;
%                               2 for FRP rupture
%                               3 for Concrete crush;
%         isSteelYielding ---  1 (tensile steel yielding before
%                                 concrete crush or FRP failure, full ductility)
%                              2 (tensile steel yielding before
%                                 concrete crush or FRP failure, partial ductility)
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

PHI1 = 1.0;
PHI2 = 0.72;
ES_LIMIT = 0.005;
STRAIN_ULTI_CONCRETE = 0.0035;
STRAIN_INITIAL_FRP = 0;
MAX_INITIAL_GUESS = 10; % maximum number of initial guess for the depth of stress block
RO_CYLINDE_2_CUBE = 0.8; % strength ratio of cylinders to cubes
ANCHOR_STRENGTH = 0;
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
        fcuLimit = 60*FC_BIAS;
        % ratio of in-situ to cylinder strength
        roInsitu = 1.00;
    case{'DESIGN_VALUE'}
        % constant partial safte factors
        gammaConcrete = 1.50;
        gammaSteel = 1.15;
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
        fctMPA = 0.5*sqrt(fcuMPA);
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
        fcuLimit = 60; 
        % ratio of in-situ to cylinder strength
        roInsitu = 0.8375;        
    otherwise
end

%% determine debonding strain of FRP
nCase = length(bBeamMM);
failMode = ones( nCase , 1 );
isSteelYielding = ones(nCase, 1); 

eFrpIni = STRAIN_INITIAL_FRP;
fAnchor = ANCHOR_STRENGTH;
eConcUlti = zeros(nCase, 1);
eConcUlti( fcuMPA<=fcuLimit ) = STRAIN_ULTI_CONCRETE;
eConcUlti( fcuMPA>fcuLimit ) = STRAIN_ULTI_CONCRETE - 0.00006*sqrt( fcuMPA( fcuMPA>fcuLimit )-fcuLimit );

betaW = sqrt( (2-bFrpMM./bBeamMM) ./ (1+bFrpMM./bBeamMM) );
% betaW = sqrt( (2.25-bFrpMM./bBeamMM) ./ (1.25+bFrpMM./bBeamMM) );
% tMax = 1.5*betaW.*(fctMPA ./ gammaConcrete);
tMax = 1.5*betaW.*(fctMPA);
Lee = 0.228*sqrt(EfrpMPA.*tFrpMM);
Ld = shearMM - aFrpMM;
alpha = 3.41*Lee./Ld;
fFrpDebond = 0.114*(4.41-alpha).*tMax.*sqrt(EfrpMPA./tFrpMM) + fAnchor;
fFrpCritical = min(fFrpMPA ./ gammaFrp, fFrpDebond./gammaBond);
isRupture = (fFrpMPA/gammaFrp) < (fFrpDebond/gammaBond);
failMode(isRupture|isAnchor) = 2; 
fFrpCritical(isRupture|isAnchor) = fFrpMPA(isRupture|isAnchor)/gammaFrp;
Ec = 3460*sqrt(fcuMPA/gammaConcrete) + 3210;
% 0.67 = 0.8 * 0.85 = cylinder/cube * in-situ/cylinder, therefore, 0.67 is
% replaced by 0.8 in the following equation
switch FLAG
    case {'MODEL_ERROR'}
        eCo = 2*0.8*(fcuMPA/gammaConcrete) ./ Ec; 
    case {'DESIGN_VALUE'}
        eCo = roInsitu*2*0.8*(fcuMPA/gammaConcrete) ./ Ec; 
    otherwise
end
%% determine c
eFrpCritical = fFrpCritical./EfrpMPA;
% cIni = dCmpMM + (dBeamMM - dCmpMM)/MAX_INITIAL_GUESS;
% cIni = dCmpMM;

% [c, fval, exitflag] = fsolve( @(c) force_balance(eCo, Ec, EsMPA, EsCmpMPA, EfrpMPA,...
%     eFrpIni, eConcUlti, eFrpCritical, bBeamMM, bFlangeMM, tFlangeMM, dBeamMM, dCmpMM,...
%     hBeamMM, c, (fsMPA/gammaSteel), (fsCmpMPA/gammaSteel), (fcuMPA/gammaConcrete),...
%     areaSteelMM2, areaFrpMM2, areaSteelCmpMM2),...
%     cIni, optimset('Display','off'));
c = zeros(nCase, 1);
fval = zeros(nCase, 1);
exitflag = zeros(nCase, 1);
for iCase = 1:nCase
% for iCase = 88
    cIni = linspace(0.1*dCmpMM(iCase), dBeamMM(iCase), MAX_INITIAL_GUESS);
    for iC = 1:MAX_INITIAL_GUESS
        [c(iCase), fval(iCase), exitflag(iCase)] = fsolve( @(c) force_balance(eCo(iCase), Ec(iCase), EsMPA(iCase), EsCmpMPA(iCase), EfrpMPA(iCase),...
            eFrpIni, eConcUlti(iCase), eFrpCritical(iCase), bBeamMM(iCase), bFlangeMM(iCase), tFlangeMM(iCase), dBeamMM(iCase), dCmpMM(iCase),...
            hBeamMM(iCase), c, (fsMPA(iCase)/gammaSteel), (fsCmpMPA(iCase)/gammaSteel), (fcuMPA(iCase)/gammaConcrete),...
            areaSteelMM2(iCase), areaFrpMM2(iCase), areaSteelCmpMM2(iCase), roInsitu),...
            cIni(iC), optimset('Display','off'));
        if (exitflag(iCase) > 0 || abs(fval(iCase))<=CONVERGE_CRITERION) && c(iCase)<=hBeamMM(iCase)
            break
        end        
    end
    if iC > MAX_INITIAL_GUESS
        fprintf('Failed convergence at %d, exitflag = %d, c=%.2f\n', i_test, exitflag, c)
    end
end
%% other approaches
%     if dCmpMM(i_test) > k2_i * c
%         % Approach 1: neglect compressive steel and assume c = 2d_c
% %        areaSteelCmpMM2(i_test) = 0;
% %        c = dCmpMM(i_test) ./ k2_i;
%         % Approach 2: neglect compressive steel and determine c again
% %         for i_c_ini = 1:10
% %             % assuming no compressive steel
% %             [c, fval, exitflag] = fsolve( @(c) force_balance(alpha1(i), beta1(i), EsMPA(i), EsCmpMPA(i), EfrpMPA(i),...
% %                 e_bi, eConcUlti(i), eFrpCritical(i), bBeamMM(i), bFlangeMM(i), tFlangeMM(i), dBeamMM(i), dCmpMM(i),...
% %                 hBeamMM(i), c, fsMPA(i), fsCmpMPA(i), fcMPA(i), areaSteelMM2(i), areaFrpMM2(i), 0),...
% %                 cIni(i_c_ini), optimset('Display','off'));
% %             x(i, 1) = c; x(i, 2) = exitflag;
% %             if exitflag == 1 || abs(fval)<=1e-5
% %                 break
% %             end
% %         end
%        % Approach 3: consider the effects of compressive steel --- Do
%        % nothing
%     end
%     % uncomment if approach 2 is used
% %     if i_c_ini == 11
% %         fprintf('Failed convergence at %d, exitflag = %d, c=%.2f\n', i, exitflag, c)
% %     end  
%%
[ffe, eFrpEff, failMode] = frp_stress(EfrpMPA, eConcUlti, hBeamMM, c, eFrpIni, eFrpCritical, failMode);
[fs, es, isSteelYielding] = tension_steel(EsMPA, eFrpEff, eFrpIni, dBeamMM, hBeamMM, c, (fsMPA/gammaSteel), isSteelYielding);
eConcFail = (eFrpEff + eFrpIni ) .* c ./ (hBeamMM-c);
eConcFail( eConcFail>eConcUlti ) = eConcUlti( eConcFail>eConcUlti );
k2 = 0.33+0.045*eConcFail./eCo;
fsc = comp_steel(eCo, eConcUlti, EsCmpMPA, eFrpEff, eFrpIni, dCmpMM, hBeamMM, c, (fsCmpMPA/gammaSteel));
switch FLAG
    case {'MODEL_ERROR'}
        phi = 1;
        momentTotalKNM = 1e-6 * phi .* ( areaSteelMM2.*fs.*(dBeamMM-k2.*c) + areaFrpMM2.*ffe.*(hBeamMM-k2.*c) + (k2.*c-dCmpMM).*areaSteelCmpMM2.*fsc );    
    case {'DESIGN_VALUE'}
        phi = ones(nCase, 1);
%         esy = fsMPA./(gammaSteel*EsMPA);
%         isDuctile = es>ES_LIMIT;
%         isLessDuctile = (es>=esy) & es<=ES_LIMIT;
%         isBrittle = es<esy;
%         phi(isDuctile)  = PHI1;
%         phi(isLessDuctile) = PHI2 + 0.28*(es(isLessDuctile)-esy(isLessDuctile))./...
%             (ES_LIMIT-esy(isLessDuctile));
%         phi(isBrittle) = PHI2;
%         isSteelYielding(isLessDuctile) = 2;
%         psi_f = ones(nCase,1);
%         psi_f(isRupture|isAnchor) = 0.55;
        momentTotalKNM = 1e-6 * phi .* ( areaSteelMM2.*fs.*(dBeamMM-k2.*c) + psi_f.*areaFrpMM2.*ffe.*(hBeamMM-k2.*c) + (k2.*c-dCmpMM).*areaSteelCmpMM2.*fsc );
    otherwise
end

return
end

% function [r,f,k1] = force_balance(e_co, Ec, Es, Es_c, Efrp, e_ini, e_cu, e_fd, b, bf, tf, d, d_c, df, c, fyd, fy_cd, fcud, As, Afrp, As_c)
function r = force_balance(e_co, Ec, Es, Es_c, Efrp, e_ini, e_cu, e_fd, b, bf, tf, d, d_c, df, c, fyd, fy_cd, fcud, As, Afrp, As_c, roInsitu)
[nRow, nCol] = size(c);
cSave = c;
r = zeros(nRow, nCol);
for iCol = 1:nCol
    c = cSave(:, iCol);
    % balance residual
    h = df;
    nCase = length(Es);
    k1 = zeros(nCase, 1);
    tmp = zeros(nCase,1);
    Ac = zeros(nCase, 1);
    
    [ffe, eFrpEff, ~] = frp_stress(Efrp, e_cu, df, c, e_ini, e_fd, tmp);
    [fs, ~] = tension_steel(Es, eFrpEff, e_ini, d, df, c, fyd, tmp);
    
    % calculate concrete strain for FRP control cases. From subfunction
    % frp_stress, the following equation equals to e_cu for concrete control
    % cases.
    e_cf = (eFrpEff + e_ini ) .* c ./ (h-c);
    k1( e_cf>0 & e_cf<=e_co ) = (Ec( e_cf>0 & e_cf<=e_co ).*e_cf( e_cf>0 & e_cf<=e_co )./...
        (2*(0.8*roInsitu)*fcud( e_cf>0 & e_cf<=e_co ))).*(1-Ec( e_cf>0 & e_cf<=e_co ).*e_cf( e_cf>0 & e_cf<=e_co )./...
        (6*(0.8*roInsitu)*fcud( e_cf>0 & e_cf<=e_co )));
    k1(e_cf>e_co&e_cf<=e_cu) = 1 - e_co(e_cf>e_co&e_cf<=e_cu)/3./e_cf(e_cf>e_co&e_cf<=e_cu);
    
    fsc = comp_steel(e_co, e_cu, Es_c, eFrpEff, e_ini, d_c, df, c, fy_cd);
    Ac(bf == b) = b(bf == b) .*c(bf == b);
    Ac( c<=tf & bf~=b ) = bf(c<=tf & bf~=b) .* c( c<=tf & bf~=b ); % Type I T section
    Ac( c>tf & bf~=b ) = b( c>tf & bf~=b ) .* c( c>tf & bf~=b ) + ( bf( c>tf & bf~=b )-b( c>tf & bf~=b )) .*tf( c>tf & bf~=b ); % Type II T section
    r = As.*fs + Afrp.*ffe - As_c.*fsc - k1.*((0.8*roInsitu)*fcud).*Ac;
end

return
end

function [ss, e_s, isSteelYielding] = tension_steel(Es, e_fe, e_ini, d, df, c, fyd, isSteelYielding)
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

function [sf, e_fe, failMode] = frp_stress(Efrp, e_cu, df, c, e_ini, e_fd, failMode)
% determine FRP stress
e_fe = e_cu.*(df-c)./c - e_ini;
failMode( e_fe<e_fd ) = 3;
e_fe( e_fe>e_fd ) = e_fd( e_fe>e_fd );
sf = Efrp .* e_fe;

return
end