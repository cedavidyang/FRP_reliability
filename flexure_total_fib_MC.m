function [momentTotalKNM, failMode, isSteelYielding]= flexure_total_fib_MC(iDesignCase, hSmp, fcSmp, fctSmp, EFrpSmp, fFrpSmp, fsSmp, areaSteelSmp)
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

STRAIN_INITIAL_FRP = 0;
MAX_INITIAL_GUESS = 10; % maximum number of initial guess for the depth of stress block
CONVERGE_CRITERION = 1E-5;

load tmpdata
%% data extracted from database

nCase = length(hSmp);
failMode = ones( nCase , 1 );
isSteelYielding = ones(nCase, 1); 
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
fckMPA = fcMPA-8; fckMPA(fckMPA<0)=0;
fctMPA = 0.3 .* fckMPA.^(2/3.0);
fctkMPA =0.7.*fctMPA;
fcm = fcMPA;
fctm = fctMPA;
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
isAnchor = logical(ANCHOR_DESIGN_ARRAY(iDesignCase) * ones(nCase, 1));
% ratio of in-situ to cylinder strength
roInsitu = 0.85;

%% determine effective strain of FRP
kbl = 2.0;
lEff = sqrt( EfrpMPA.*tFrpMM ./ (kbl.*fctm) );
lb = shearMM - aFrpMM;
betaL = lb./lEff .* (2-lb./lEff);
betaL( lb>lEff ) = 1.0;
tmp = (2-bFrpMM./bBeamMM)./(1+bFrpMM./bBeamMM);
tmp( tmp<0 ) = 0;
kb = sqrt(tmp);
kb(kb<1) = 1.0;
kk = 0.25; kc = 2.0;
edeb = kc*(kk./gammaBond).*kb.*betaL.*sqrt(2*fcm.^(2/3.0)./(EfrpMPA.*tFrpMM));
efu = (fFrpMPA./gammaFrp)./EfrpMPA;
efe = min(edeb, efu);
efe(isAnchor) = efu(isAnchor);

%% determine c
eFrpIni = STRAIN_INITIAL_FRP;
eConcUlti = 0.0035;
eCo = 0.002;
eFrpCritical = efe;
c = zeros(nCase, 1);
fval = zeros(nCase, 1);
exitflag = zeros(nCase, 1);
for iCase = 1:nCase
    cIni = linspace(0.1*dCmpMM(iCase), dBeamMM(iCase), MAX_INITIAL_GUESS);
    for iC = 1:MAX_INITIAL_GUESS
        [c(iCase), fval(iCase), exitflag(iCase)] = fsolve( @(c) force_balance(eCo, EsMPA(iCase), EsCmpMPA(iCase), EfrpMPA(iCase),...
            eFrpIni, eConcUlti, eFrpCritical(iCase), bBeamMM(iCase), bFlangeMM(iCase), tFlangeMM(iCase), dBeamMM(iCase), dCmpMM(iCase),...
            hBeamMM(iCase), c, (fsMPA(iCase)/gammaSteel), (fsCmpMPA(iCase)/gammaSteel), (fcMPA(iCase)/gammaConcrete),...
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

%% determine capacity

[ffe, eFrpEff, failMode] = frp_stress(EfrpMPA, eConcUlti, hBeamMM, c, eFrpIni, eFrpCritical, failMode);
[fs, es, isSteelYielding] = tension_steel(EsMPA, eFrpEff, eFrpIni, dBeamMM, hBeamMM, c, (fsMPA/gammaSteel), isSteelYielding);

deltaG = zeros(nCase,1);
eConcFail = (eFrpEff + eFrpIni ) .* c ./ (hBeamMM-c);
indx1 = eConcFail>0&eConcFail<=eCo;
deltaG(indx1) = (8-1000*eConcFail(indx1))./(4.*(6-1000*eConcFail(indx1)));
indx2 = eConcFail>eCo&eConcFail<=eConcUlti;
deltaG(indx2) = (1000*eConcFail(indx2).*(2000*eConcFail(indx2)-4)+2)./...
    (2000*eConcFail(indx2).*(3000*eConcFail(indx2)-2));
deltaG( eConcFail>eConcUlti ) = (1000*eConcUlti.*(2000*eConcUlti-4)+2)./(2000*eConcUlti.*(3000*eConcUlti-2));

fsc = comp_steel(eCo, eConcUlti, EsCmpMPA, eFrpEff, eFrpIni, dCmpMM, hBeamMM, c, (fsCmpMPA/gammaSteel));

phi=1;
momentTotalKNM = 1e-6 * phi .* ( areaSteelMM2.*fs.*(dBeamMM-deltaG.*c) + areaFrpMM2.*ffe.*(hBeamMM-deltaG.*c) + (deltaG.*c-dCmpMM).*areaSteelCmpMM2.*fsc );

return
end

function r = force_balance(e_co, Es, Es_c, Efrp, e_ini, e_cu, e_fd, b, bf, tf, d, d_c, df, c, fyd, fy_cd, fcd, As, Afrp, As_c, roInsitu)
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
psi = zeros(nCase,1);
e_cf = (eFrpEff + e_ini ) .* c ./ (h-c);
indx1 = e_cf>0&e_cf<=e_co;
psi(indx1) = 1000*e_cf(indx1).*(0.5-1000/12*e_cf(indx1));
indx2 = e_cf>e_co&e_cf<=e_cu;
psi(indx2) = 1-2./3000./e_cf(indx2);
psi(e_cf>e_cu) = 1-2./3000./e_cu;

fsc = comp_steel(e_co, e_cu, Es_c, eFrpEff, e_ini, d_c, df, c, fy_cd);
Ac(bf == b) = b(bf == b) .*c(bf == b);
Ac( c<=tf & bf~=b ) = bf(c<=tf & bf~=b) .* c( c<=tf & bf~=b ); % Type I T section
Ac( c>tf & bf~=b ) = b( c>tf & bf~=b ) .* c( c>tf & bf~=b ) + ( bf( c>tf & bf~=b )-b( c>tf & bf~=b )) .*tf( c>tf & bf~=b ); % Type II T section
r = As.*fs + Afrp.*ffe - As_c.*fsc - roInsitu.*psi.*fcd.*Ac;

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