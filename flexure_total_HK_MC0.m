function resistSmp = flexure_total_HK_MC(iDesignCase, hSmp, fcSmp, fctSmp, EFrpSmp, fFrpSmp, fsSmp, areaSteelSmp)
% Compute maximum moment using Hong Kong guideline
% Input: FLAG --- type of the analysis
%        factorFrp -- partial safety factor or resistance reduction factor
% Output: momentTotalKNM --- predicted maximum moment
%         failMode --- fail mode: 1 for IC debonding;
%                               2 for FRP rupture
%                               3 for Concrete crush;
%         isSteelYielding ---  1 (tensile steel yielding before
%                                 concrete crush or FRP failure)
%                              0 (tensile steel does not yield)

PHI1 = 1.0;
PHI2 = 0.72;
ES_LIMIT = 0.005;
STRAIN_ULTI_CONCRETE = 0.0035;
STRAIN_INITIAL_FRP = 0;
MAX_INITIAL_GUESS = 10; % maximum number of initial guess for the depth of stress block
RO_CYLINDE_2_CUBE = 0.8; % strength ratio of cylinders to cubes
ANCHOR_STRENGTH = 0;
CONVERGE_CRITERION = 1E-5;
MAX_INITIAL_GUESS = 10;

load tmpdata
%% data extracted from database
nCase = length(hSmp);

% constant partial safte factors
gammaConcrete = 1.00;
gammaSteel = 1.00;
gammaFrp = 1.00;
gammaBond = 1.00;
% geometric properties
hBeamMM = hSmp;
as = H_DESIGN_ARRAY_MM(iDesignCase) - D_DESIGN_ARRAY_MM(iDesignCase);
dBeamMM = hSmp - as;
h2b = H_DESIGN_ARRAY_MM(iDesignCase)/B_DESIGN_ARRAY_MM(iDesignCase);
bBeamMM = hSmp/h2b;
bFlangeMM = BF_DESIGN_ARRAY_MM(iDesignCase) * ones(nCase, 1);
tFlangeMM = TF_DESIGN_ARRAY_MM(iDesignCase) * ones(nCase, 1);
dCmpMM = D_CMP_DESIGN_ARRAY_MM(iDesignCase) * ones(nCase, 1);
aFrpMM = FRP_END_DESIGN_ARRAY_MM(iDesignCase) * ones(nCase, 1);
shearMM = SHEAR_DESIGN_ARRAY_MM(iDesignCase) * ones(nCase, 1);
% concrete properties
fcMPA = fcSmp;
fcuMPA = fcMPA / RO_CYLINDE_2_CUBE;
fctMPA = fctSmp;
% steel reinforcement
EsMPA = ES_DESIGN_ARRAY_MPA(iDesignCase) * ones(nCase, 1);
fsMPA = fsSmp;
areaSteelMM2 = areaSteelSmp;
EsCmpMPA = ES_CMP_DESIGN_ARRAY_MPA(iDesignCase) * ones(nCase, 1);
fsCmpMPA = FS_CMP_DESIGN_ARRAY_MPA(iDesignCase) * ones(nCase, 1);
areaSteelCmpMM2 = AREA_STEEL_CMP_DESIGN_ARRAY_MM2(iDesignCase) * ones(nCase, 1);
% FRP strengthening
EfrpMPA = EFrpSmp;
fFrpMPA = fFrpSmp;
tFrpMM = T_FRP_DESIGN_ARRAY_MM(iDesignCase) * ones(nCase, 1);
bFrpMM = B_FRP_DESIGN_ARRAY_MM(iDesignCase) * ones(nCase, 1);
areaFrpMM2 = bFrpMM .* tFrpMM;
isAnchor = ANCHOR_DESIGN_ARRAY(iDesignCase) * ones(nCase, 1);
% limit values
fcuLimit = 60*FC_BIAS;   
% ratio of in-situ to cylinder strength
roInsitu = 0.8375;

%% determine debonding strain of FRP
nCase = length(bBeamMM);
failMode = ones( nCase , 1 );
isSteelYielding = ones(nCase, 1); 

eFrpIni = STRAIN_INITIAL_FRP;
fAnchor = ANCHOR_STRENGTH;
eConcUlti = STRAIN_ULTI_CONCRETE;
isHigh = mean(fcuMPA)>fcuLimit;
eConcUlti( isHigh ) = STRAIN_ULTI_CONCRETE - 0.00006*sqrt( mean(fcuMPA)-fcuLimit );
% eConcUlti = zeros(nCase, 1);
% eConcUlti( fcuMPA<=fcuLimit ) = 0.0035;
% eConcUlti( fcuMPA>fcuLimit ) = 0.0035 - 0.00006*sqrt( fcuMPA( fcuMPA>fcuLimit )-fcuLimit );

betaW = sqrt( (2-bFrpMM./bBeamMM) ./ (1+bFrpMM./bBeamMM) );
% betaW = sqrt( (2.25-bFrpMM./bBeamMM) ./ (1.25+bFrpMM./bBeamMM) );
tMax = 1.5*betaW.*(fctMPA ./ gammaConcrete);
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
% replaced by 0.8 in model error analysis
eCo = 2*(0.8*roInsitu)*(fcuMPA/gammaConcrete) ./ Ec; 
%% determine c
eFrpCritical = fFrpCritical./EfrpMPA;
% cIni = dCmpMM + (dBeamMM - dCmpMM)/MAX_INITIAL_GUESS;
% cIni = dCmpMM;
% cIni = 2*dCmpMM;
% 
% [c, fval, exitflag] = fsolve( @(c) force_balance(eCo, Ec, EsMPA, EsCmpMPA, EfrpMPA,...
%     eFrpIni, eConcUlti, eFrpCritical, bBeamMM, bFlangeMM, tFlangeMM, dBeamMM, dCmpMM,...
%     hBeamMM, c, (fsMPA/gammaSteel), (fsCmpMPA/gammaSteel), (fcuMPA/gammaConcrete),...
%     areaSteelMM2, areaFrpMM2, areaSteelCmpMM2, roInsitu),...
%     cIni, optimset('Display','off'));

cIniArray = linspace(0.1*dCmpMM, dBeamMM, MAX_INITIAL_GUESS);
for iC = 1:MAX_INITIAL_GUESS
    cIni = cIniArray(iC);
    [c, fval, exitflag] = fsolve( @(c) force_balance(eCo, Ec, EsMPA, EsCmpMPA, EfrpMPA,...
        eFrpIni, eConcUlti, eFrpCritical, bBeamMM, bFlangeMM, tFlangeMM, dBeamMM, dCmpMM,...
        hBeamMM, c, (fsMPA/gammaSteel), (fsCmpMPA/gammaSteel), (fcuMPA/gammaConcrete),...
        areaSteelMM2, areaFrpMM2, areaSteelCmpMM2, roInsitu),...
        cIni, optimset('Display','off'));
    if exitflag > 0 && abs(fval)<=CONVERGE_CRITERION && c<=hBeamMM
        break
    end
end
if iC > MAX_INITIAL_GUESS
    fprintf('Failed convergence at %d, exitflag = %d, c=%.2f\n', i_test, exitflag, c)
end
    
%% useless older versions
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
eConcFail( eConcFail>eConcUlti ) = eConcUlti;
k2 = 0.33+0.045*eConcFail./eCo;
fsc = comp_steel(eCo, eConcUlti, EsCmpMPA, eFrpEff, eFrpIni, dCmpMM, hBeamMM, c, (fsCmpMPA/gammaSteel));

phi = 1;
resistSmp = 1e-6 * phi .* ( areaSteelMM2.*fs.*(dBeamMM-k2.*c) + areaFrpMM2.*ffe.*(hBeamMM-k2.*c) + (k2.*c-dCmpMM).*areaSteelCmpMM2.*fsc );

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
%     f = 0.5 * ( (As.*fs + Afrp.*ffe) + (As_c.*fsc + k1.*(0.8*fcud).*Ac) );
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