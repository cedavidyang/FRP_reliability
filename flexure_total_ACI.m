function [momentTotalKNM, failMode, isSteelYielding]= flexure_total_ACI(FLAG, factorFrp)
% Compute maximum moment using ACI 440.2R-08
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

PHI1 = 0.9;
PHI2 = 0.65;
ES_LIMIT = 0.005;
STRAIN_ULTI_CONCRETE = 0.003;
STRAIN_INITIAL_FRP = 0;
MAX_INITIAL_GUESS = 10; % maximum number of initial guess for the depth of stress block
CONVERGE_CRITERION = 1E-5;

load tmpdata
%% data extracted from database
switch FLAG
    case{'MODEL_ERROR'}
        % geometric data
        bBeamMM = B_TEST_ARRAY_MM;
        hBeamMM = H_TEST_ARRAY_MM;
        bFlangeMM = BF_TEST_ARRAY_MM;
        tFlangeMM = TF_TEST_ARRAY_MM;
        dBeamMM = D_TEST_ARRAY_MM;
        dCmpMM = D_CMP_TEST_ARRAY_MM;
        % concrete properties
        fcMPA = FC_TEST_ARRAY_MPA;
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
        fcLimit = 28*FC_BIAS;
        % ratio of in-situ to cylinder strength
        roInsitu = 1.00;    
    case{'DESIGN_VALUE'}
        % geometric data
        bBeamMM = B_DESIGN_ARRAY_MM;
        hBeamMM = H_DESIGN_ARRAY_MM;
        bFlangeMM = BF_DESIGN_ARRAY_MM;
        tFlangeMM = TF_DESIGN_ARRAY_MM;
        dBeamMM = D_DESIGN_ARRAY_MM;
        dCmpMM = D_CMP_DESIGN_ARRAY_MM;
        % concrete properties
        fcMPA = FC_DESIGN_ARRAY_MPA;
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
        fcLimit = 28;
        % ratio of in-situ to cylinder strength
        roInsitu = 0.85;        
    otherwise
end

%% determine debonding strain of FRP
nCase = length(bBeamMM);
failMode = ones(nCase , 1 );
eFrpUlti = fFrpMPA ./ EfrpMPA;
eDebond = 0.41 .* sqrt( fcMPA ./ (EfrpMPA.*tFrpMM) );
isRupture = eDebond > 0.9*eFrpUlti;
failMode( isRupture|isAnchor ) = 2;
eFrpCritical = zeros(nCase, 1);
eFrpCritical(isRupture|isAnchor) = 0.9*eFrpUlti( isRupture|isAnchor );
eFrpCritical((~isRupture)&(~isAnchor)) = eDebond( (~isRupture)&(~isAnchor) );
%% determine c
% ACI 318-11 10.2.7.3
% 0.85 shouldn't be used in Model Error, but should be used in Design values
alpha1 = roInsitu;
% depth of the stress block
beta1 = 0.85*ones(nCase, 1);
beta1( fcMPA>fcLimit ) = 0.85 - ( fcMPA(fcMPA>fcLimit) - fcLimit )/7 *0.05;
beta1( beta1<0.65 ) = 0.65;

eConcUlti =STRAIN_ULTI_CONCRETE;
eFrpIni = STRAIN_INITIAL_FRP;
isSteelYielding = ones(nCase, 1); 

% cIni = dCmpMM + (dBeamMM - dCmpMM)/MAX_INITIAL_GUESS;
% 
% index = 1:nCase;
% [c, fval, exitflag] =...
%     fsolve( @(c) force_balance(alpha1, beta1(index), EsMPA(index), EsCmpMPA(index),...
%     EfrpMPA(index),eFrpIni, eConcUlti, eFrpCritical(index),...
%     bBeamMM(index), bFlangeMM(index), tFlangeMM(index),...
%     dBeamMM(index), dCmpMM(index), hBeamMM(index), c,...
%     fsMPA(index), fsCmpMPA(index), fcMPA(index),...
%     areaSteelMM2(index), areaFrpMM2(index),areaSteelCmpMM2(index)),...
%     cIni(index), optimset('Display','off'));
c = zeros(nCase, 1);
fval = zeros(nCase, 1);
exitflag = zeros(nCase, 1);
for iCase = 1:nCase
    cIni = linspace(0.1*dCmpMM(iCase), dBeamMM(iCase), MAX_INITIAL_GUESS);
    for iC = 1:MAX_INITIAL_GUESS
        [c(iCase), fval(iCase), exitflag(iCase)] = fsolve( @(c) force_balance(alpha1, beta1(iCase), EsMPA(iCase), EsCmpMPA(iCase),...
            EfrpMPA(iCase),eFrpIni, eConcUlti, eFrpCritical(iCase),...
            bBeamMM(iCase), bFlangeMM(iCase), tFlangeMM(iCase),...
            dBeamMM(iCase), dCmpMM(iCase), hBeamMM(iCase), c,...
            fsMPA(iCase), fsCmpMPA(iCase), fcMPA(iCase),...
            areaSteelMM2(iCase), areaFrpMM2(iCase),areaSteelCmpMM2(iCase)),...
            cIni(iC), optimset('Display','off'));
        if exitflag(iCase) > 0 || abs(fval(iCase))<=CONVERGE_CRITERION && c(iCase)<=hBeamMM(iCase)
            break
        end        
    end
    if iC > MAX_INITIAL_GUESS
        fprintf('Failed convergence at %d, exitflag = %d, c=%.2f\n', i_test, exitflag, c)
    end
%% other approaches
% if dCmpMM(iCase) > beta1(iCase) .* c/2
%     % Approach 1: neglect compressive steel and assume c = 2d_c
%     %         As_c(i) = 0;
%     %         c = 2*d_c(i) ./ beta1(i);
%     % Approach 2: neglect compressive steel and determine c again
%     %         for i_c_ini = 1:10
%     %             % assuming no compressive steel
%     %             [c, fval, exitflag] = fsolve( @(c) force_balance(alpha1, beta1(i), Es(i), Es_c(i), Efrp(i),...
%     %                 e_bi, e_cu, e_fd(i), b(i), bf(i), tf(i), d(i), d_c(i),...
%     %                 h(i), c, fyk(i), fy_ck(i), fck(i), As(i), Afrp(i), 0),...
%     %                 c_ini(i_c_ini), optimset('Display','off'));
%     %             if exitflag == 1 || abs(fval)<=1e-5
%     %                 break
%     %             end
%     %         end
%     % Approach 3: consider the effects of compressive steel --- Do nothing
%     %        disp('Compressive rebar is too deep')
% end
% % uncomment if approach 2 is used
% %     if i_c_ini == 11
% %         fprintf('Failed convergence at %d, exitflag = %d, c=%.2f\n', i, exitflag, c)
% %     end    
end

[fFrpEff, eFrpEff, failMode] = frp_stress(EfrpMPA, eConcUlti, hBeamMM, c, eFrpIni, eFrpCritical, failMode);
[fs, es, isSteelYielding] = tension_steel(EsMPA, eFrpEff, eFrpIni, dBeamMM, hBeamMM, c, fsMPA, isSteelYielding);
fsc = comp_steel(beta1, EsCmpMPA, eFrpEff, eFrpIni, dCmpMM, hBeamMM, c, fsCmpMPA);

% cr = beta1.*c;
% isCmpSteel2Deep = dCmpMM > beta1 .* c /2;
% cr( isCmpSteel2Deep ) = 2 * dCmpMM( isCmpSteel2Deep );
% areaConcBlock = zeros(nCase, 1);
% isRec = bFlangeMM == bBeamMM;
% areaConcBlock(isRec) = bBeamMM(isRec) .*cr(isRec);
% isT1 = cr<=tFlangeMM & bFlangeMM~=bBeamMM;
% areaConcBlock(isT1) = bFlangeMM(isT1) .* cr(isT1);
% isT2 = cr>tFlangeMM & bFlangeMM~=bBeamMM;
% areaConcBlock(isT2) = bBeamMM(isT2) .* cr(isT2) + ( bFlangeMM(isT2)-bBeamMM(isT2)).*tFlangeMM(isT2);

% dFrp = hBeamMM;
% r = force_balance(alpha1, beta1, EsMPA, EsCmpMPA, EfrpMPA, eFrpIni, eConcUlti, eFrpCritical, bBeamMM, bFlangeMM, tFlangeMM, dBeamMM, dCmpMM, dFrp, c, fsMPA, fsCmpMPA, fcMPA, areaSteelMM2, areaFrpMM2, areaSteelCmpMM2);
switch FLAG
    case {'MODEL_ERROR'}
        psif = 1; phi = 1;
        momentTotalKNM = 1e-6 * phi .* ( areaSteelMM2.*fs.*(dBeamMM-beta1.*c/2) + psif*areaFrpMM2.*fFrpEff.*(hBeamMM-beta1.*c/2) + (beta1.*c/2-dCmpMM).*areaSteelCmpMM2.*fsc );  
%         Mpre = 1e-6 * phi .* ( alpha1*fck.*Ac.*(d-cr/2) + psif*Afrp.*ffe.*(df-d) + fsc.*As_c.*(d-d_c) );
    case {'DESIGN_VALUE'}
        psif = factorFrp;
        phi = zeros(nCase, 1);
        esy = (fsMPA./EsMPA);
        isDuctile = es>ES_LIMIT;
        isLessDuctile = (es>=esy) & es<=ES_LIMIT;
        isBrittle = es<esy;
        phi(isDuctile)  = PHI1;
        phi(isLessDuctile) = PHI2 + 0.25*(es(isLessDuctile)-esy(isLessDuctile))./...
                            (ES_LIMIT-esy(isLessDuctile));
        phi(isBrittle) = PHI2;
        isSteelYielding(isLessDuctile) = 2;
        momentTotalKNM = 1e-6 * phi .* ( areaSteelMM2.*fs.*(dBeamMM-beta1.*c/2) + psif*areaFrpMM2.*fFrpEff.*(hBeamMM-beta1.*c/2) + (beta1.*c/2-dCmpMM).*areaSteelCmpMM2.*fsc );
%         Mpre = 1e-6 * phi .* ( alpha1*fck.*Ac.*(d-cr/2) +
%         psif*Afrp.*ffe.*(df-d) + fsc.*As_c.*(d-d_c) );
    otherwise
end

return
end

function r = force_balance(alpha1, beta1, Es, Es_c, Efrp, e_bi, e_cu, e_fd, b, bf, tf, d, d_c, df, c, fyk, fy_ck, fck, As, Afrp, As_c)
[nRow, nCol] = size(c);
cSave = c;
r = zeros(nRow, nCol);
for iCol = 1:nCol
    c = cSave(:, iCol);
    % balance residual
    nCase = length(Es);
    tmp = zeros(nCase,1);
    [ffe, eFrpEff, ~] = frp_stress(Efrp, e_cu, df, c, e_bi, e_fd, tmp);
    [fs, ~] = tension_steel(Es, eFrpEff, e_bi, d, df, c, fyk, tmp);
    fsc = comp_steel(beta1, Es_c, eFrpEff, e_bi, d_c, df, c, fy_ck);
    cr = c.*beta1;
    Ac = zeros(nCase, 1);
    Ac(bf == b) = b(bf == b) .*cr(bf == b);
    Ac( cr<=tf & bf~=b ) = bf(cr<=tf & bf~=b) .* cr( cr<=tf & bf~=b ); % Type I T section
    Ac( cr>tf & bf~=b ) = b( cr>tf & bf~=b ) .* cr( cr>tf & bf~=b ) + ( bf( cr>tf & bf~=b )-b( cr>tf & bf~=b )) .*tf( cr>tf & bf~=b ); % Type II T section
    r(:, iCol) = As.*fs + Afrp.*ffe - As_c.*fsc - alpha1.*fck.*Ac;
end

return
end

function [ss, e_s, isSteelYielding] = tension_steel(Es, e_fe, e_bi, d, df, c, fyk, isSteelYielding)
% determine tensile steel stress
e_s = (e_fe+e_bi).*(d-c)./(df-c);
ss = Es .* e_s;
isSteelYielding( ss<fyk ) = 0;
ss( ss>fyk ) = fyk( ss>fyk );

return
end

function ss = comp_steel(beta1, Es_c, e_fe, e_bi, d_c, df, c, fy_ck)
% determine compressive steel stress
e_sc = (c-d_c) ./ (df-c) .* (e_fe+e_bi);
ss = Es_c .* e_sc;
ss( ss>fy_ck ) = fy_ck( ss>fy_ck );
ss( ss<-fy_ck ) = -fy_ck( ss<-fy_ck  );
% ss( d_c>beta1.*c/2 ) = 0;
return
end

function [sf, e_fe, f_mode] = frp_stress(Efrp, e_cu, df, c, e_bi, e_fd, f_mode)
% determine FRP stress
e_fe = e_cu.*(df-c)./c - e_bi;
f_mode( e_fe<e_fd ) = 3;
e_fe( e_fe>e_fd ) = e_fd( e_fe>e_fd );
sf = Efrp .* e_fe;

return
end
