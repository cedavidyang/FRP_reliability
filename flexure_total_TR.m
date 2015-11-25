function [momentTotalKNM, failMode, isSteelYielding]= flexure_total_TR(FLAG, factorFrp)
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

load tmpdata

%% data extracted from database
switch FLAG
    case{'MODEL_ERROR'}
        % constant partial safte factors
        gammaConcrete = 1.00;
        gammaSteel = 1.00;
        gammaEfrp = (gammaBond/psi_f).^2;
        gammaEPfrp = (gammaFrp/psi_f)./gammaEfrp;
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
        fcm = fcMPA;
        fck = fcm - 8;
        fctm = 0.30*fck.^(2/3);
        fctk = 0.7*fctm;
        fc = fcm; fct = fctm;
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
        isAnchor = logical(ANCHOR_TEST_ARRAY);       
        % ratio of in-situ to cylinder strength
        roInsitu = 1.0;          
    case{'DESIGN_VALUE'}
        % constant partial safte factors
        gammaConcrete = 1.50;
        gammaSteel = 1.15;
        gammaEfrp = (gammaBond/psi_f).^2;
        gammaEPfrp = (gammaFrp/psi_f)./gammaEfrp;    
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
        fck = fcMPA;
        fcm = fck+8;
        fctm = 0.30*fck.^(2/3);
        fctk = 0.7*fctm;
        fc = fck; fct = fctk;
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
        isAnchor = logical(ANCHOR_DESIGN_ARRAY);     
        % ratio of in-situ to cylinder strength
        roInsitu = 0.85;         % BS EN 1992-1-NA-2004 (3.1.6) 
    otherwise
end

%% beam model
nCase = length(bBeamMM);
beam_models(nCase).name = [];
for i=1:nCase
    % concrete
    beam_models(i).gammaConcrete = gammaConcrete;
%     beam_models(i).fcm =  fcm(i);
%     beam_models(i).fck =  fck(i);
    beam_models(i).fc = fc(i);
    beam_models(i).fcd =  roInsitu*fc(i)/gammaConcrete;
%     beam_models(i).fctm =  fctm(i);
%     beam_models(i).fctk =  fctk(i);
    beam_models(i).fct = fct(i);
    beam_models(i).fctd =  fct(i)/gammaConcrete;
    beam_models(i).epc0 = 1.75/1000; beam_models(i).epcu = 3.5/1000;
    % steel
    beam_models(i).fyd = fsMPA(i)/gammaSteel; beam_models(i).fycd = fsCmpMPA(i)/gammaSteel;
    beam_models(i).Es = EsMPA(i); beam_models(i).Esc = EsCmpMPA(i);
    beam_models(i).epy = fsMPA(i)/gammaSteel/EsMPA(i);
    beam_models(i).epyc = fsCmpMPA(i)/gammaSteel/EsCmpMPA(i);
    % frp
    beam_models(i).efd = fFrpMPA(i)/EfrpMPA(i)/gammaEPfrp;
    beam_models(i).Efd =  EfrpMPA(i)/gammaEfrp;
    % geometric
    beam_models(i).isT = bFlangeMM(i)~=bBeamMM(i);
    beam_models(i).d = dBeamMM(i); beam_models(i).a = dCmpMM(i); beam_models(i).h = hBeamMM(i);
    beam_models(i).hf = tFlangeMM(i); beam_models(i).bf = bFlangeMM(i); beam_models(i).b = bBeamMM(i);
    beam_models(i).Asc = areaSteelCmpMM2(i); beam_models(i).As = areaSteelMM2(i);
    beam_models(i).Af = areaFrpMM2(i);
    beam_models(i).tfrp = tFrpMM(i);
    beam_models(i).bfrp = bFrpMM(i);
    beam_models(i).ss = shearMM(i);
end
d = [beam_models.d]'; h = [beam_models.h]'; a = [beam_models.a]';
epy = [beam_models.epy]'; epyc = [beam_models.epyc]';
Efd = [beam_models.Efd]'; efd = [beam_models.efd]';
epc0 = [beam_models.epc0]'; epcu = [beam_models.epcu]';
ss = [beam_models.ss]';
% Efd = EfrpMPA/gammaEfrp;
% fcd = fck/gammaConcrete;
% epc0 = 1.75/1000;
% epcu = 3.5/1000;
% fyd = fsMPA/gammaSteel;
% fycd = fsCmpMPA/gammaSteel;
% epy = fyd./EsMPA; epyc = fycd./EsCmpMPA;
% isT = bFlangeMM~=bBeamMM;
% d = dBeamMM; a = dCmpMM; h = hBeamMM; hf = tFlangeMM; bf = bFlangeMM; b=bBeamMM;
% Esc = EsCmpMPA; Asc = areaSteelCmpMM2; Es = EsMPA; As= areaSteelMM2; Af = areaFrpMM2;

%% determine yield momemt
% symbolic computation
syms x
xnum = zeros(nCase,1);

epc = x./(d-x).*epy;
ef = (h-x)./(d-x).*epy;
esc = (x-a)./(d-x+a).*epy;
es = epy;

My = zeros(nCase,1);
xy = zeros(nCase,1);
for ipoly = 1:nCase
    [My(ipoly), xy(ipoly)] = compute_moment(beam_models(ipoly), epc(ipoly), ef(ipoly),...
        esc(ipoly), es(ipoly));
end

%% determine FRP failure moment
% determine debonding strain of FRP
edeblimit = 0.008;
efe = min(edeblimit, efd);
efe(isAnchor) = efd(isAnchor);

epc = x./(h-x).*efe;
ef = efe;

Mf = zeros(nCase,1);
xf = zeros(nCase,1);
for i = 1:nCase 
    exit_flag = 0;
    for icase=1:4
        switch icase
            case 1    % all yield
                esc = epyc(i); es = epy(i);
                [Mf(i), xf(i)] = compute_moment(beam_models(i), epc(i), ef(i),esc, es);
                if (xf(i)-a(i))/(h(i)-xf(i)+a(i))*efe(i) >= epyc(i) &&...
                        (d(i)-xf(i))/(h(i)-xf(i))*efe(i) >= epy(i)
                    exit_flag = 1;
                end
            case 2    % tension yield and compression not yield
                esc = (x-a(i))/(h(i)-x+a(i))*efe(i); es = epy(i);
                [Mf(i), xf(i)] = compute_moment(beam_models(i), epc(i), ef(i),esc, es);
                if (xf(i)-a(i))/(h(i)-xf(i)+a(i))*efe(i) < epyc(i) &&...
                        (d(i)-xf(i))/(h(i)-xf(i))*efe(i) >= epy(i)
                    exit_flag = 1;
                end
            case 3    % compression yield and tension not yield
                esc = epyc(i); es = (d(i)-x)/(h(i)-x)*efe(i);
                [Mf(i), xf(i)] = compute_moment(beam_models(i), epc(i), ef(i),esc, es);
                if (xf(i)-a(i))/(h(i)-xf(i)+a(i))*efe(i) >= epyc(i) &&...
                        (d(i)-xf(i))/(h(i)-xf(i))*efe(i) < epy(i)
                    exit_flag = 1;
                end            
            case 4
                esc = (x-a(i))./(h(i)-x+a(i)).*efe(i); es = (d(i)-x)./(h(i)-x).*efe(i);
                [Mf(i), xf(i)] = compute_moment(beam_models(i), epc(i), ef(i),esc, es);
                if (xf(i)-a(i))/(h(i)-xf(i)+a(i))*efe(i) < epyc(i) &&...
                        (d(i)-xf(i))/(h(i)-xf(i))*efe(i) < epy(i)
                    exit_flag = 1;
                end 
            otherwise
        end
        if exit_flag == 1
            break
        end
    end
end

%% determine concrete crushing moment
Mc = zeros(nCase,1);
xc = zeros(nCase,1);

epc = epcu;
ef = (h-x)./x.*epc;
for i = 1:nCase
    exit_flag = 0;
    for icase=1:4
        switch icase
            case 1    % all yield
                esc = epyc(i); es = epy(i);
                [Mc(i), xc(i)] = compute_moment(beam_models(i), epc(i), ef(i),esc, es);
                if (xf(i)-a(i))/xf(i)*epc(i) >= epyc(i) && (d(i)-xf(i))/xf(i)*epc(i) >= epy(i)
                    exit_flag = 1;
                end
            case 2    % tension yield and compression not yield
                esc = (x-a(i))/x*epc(i); es = epy(i);
                [Mc(i), xc(i)] = compute_moment(beam_models(i), epc(i), ef(i),esc, es);
                if (xf(i)-a(i))/xf(i)*epc(i) < epyc(i) && (d(i)-xf(i))/xf(i)*epc(i) >= epy(i)
                    exit_flag = 1;
                end
            case 3    % compression yield and tension not yield
                esc = epyc(i); es = (d(i)-x)/x*epc(i);
                [Mc(i), xc(i)] = compute_moment(beam_models(i), epc(i), ef(i),esc, es);
                if (xf(i)-a(i))/xf(i)*epc(i) >= epyc(i) && (d(i)-xf(i))/xf(i)*epc(i) < epy(i)
                    exit_flag = 1;
                end            
            case 4
                esc = (x-a(i))./x.*epc(i); es = (d(i)-x)./x.*epc(i);
                [Mc(i), xc(i)] = compute_moment(beam_models(i), epc(i), ef(i),esc, es);
                if (xf(i)-a(i))/xf(i)*epc(i) < epyc(i) && (d(i)-xf(i))/xf(i)*epc(i) < epy(i)
                    exit_flag = 1;
                end 
            otherwise
        end
        if exit_flag == 1
            break
        end
    end
end

%% determine failure mode and revise FRP rupture cases
M = zeros(nCase,1);
% My<Mf<Mc
nlayer = 100;
indx1 = min([My, Mc, Mf]')'==My;
Med = min([Mc, Mf]')';
dx = ss./(1+My./Med);
efy = (h-xy)./(d-xy).*epy;
sfy = Efd.*efy;
xmin = xf;
xmin(Med==Mc) = xc(Med==Mc);
efe(Med==Mf) = efe(Med==Mf);
efe(Med==Mc) = xc(Med==Mc)./(h(Med==Mc)-xc(Med==Mc)).*epcu(Med==Mc);
sfe = Efd.*efe;
offlimit = isofflimit(dx, beam_models, sfy, sfe, My, Med);
offlimit = offlimit&(~isAnchor);
M(indx1&(~offlimit)) = Med(indx1&(~offlimit));
for icase=find((indx1&offlimit)==1)'
    [M(icase), ~] = adjustM(nlayer, xmin(icase), xy(icase), efe(icase), efy(icase), My(icase), beam_models(icase));
end
% M(indx1) = Med(indx1);
% Mf<Mf<Mc and Mf<My<Mc
M( min([My, Mc, Mf]')==Mf' ) = Mf( min([My, Mc, Mf]')==Mf' );
% Mc<Mf<My and Mc<My<Mf
M( min([My, Mc, Mf]')==Mc' ) = Mc( min([My, Mc, Mf]')==Mc' );

momentTotalKNM = M;
failMode = indx1&offlimit;
isSteelYielding = zeros(nCase,1);

return
end

function [M, xnum] = compute_moment(beam_model, epc, ef, esc, es)
    syms x
    Efd = beam_model.Efd;
    fcd = beam_model.fcd;
    epc0 = beam_model.epc0;
    epcu = beam_model.epcu;
    fyd = beam_model.fyd;
    fycd = beam_model.fycd;
    epy = beam_model.epy; epyc = beam_model.epyc;
    isT = beam_model.isT;
    d = beam_model.d; a = beam_model.a; h = beam_model.h; hf = beam_model.hf;
    bf = beam_model.bf; b=beam_model.b;
    Esc = beam_model.Esc; Asc = beam_model.Asc;
    Es = beam_model.Es; As= beam_model.As; Af = beam_model.Af;
    
    Fcc = bf .* x .* 0.5 .* fcd ./ epc0 .* epc;    % Type I T-beam and epc<ep0 
    Fsc = Esc .* esc .* Asc ;
    Fst = As .* Es .* es ;
    Fft = Af .* Efd .* ef;
    xnum = getx(Fcc, Fsc, Fst, Fft, d);
    efnum = double(subs(ef, x, xnum));
    beyondEpc0 = xnum./(h-xnum).*efnum>epc0;

    if beyondEpc0
        Fcc = Fcc .* (1-((epc-epc0)./epc).^2);    % Type I T-beam and epc>ep0 
        xnum = getx(Fcc, Fsc, Fst, Fft, d);
        M = moment_with_x(2, beam_model, xnum, epc, ef, esc, es);
    else
        M = moment_with_x(1, beam_model, xnum, epc, ef, esc, es);
    end
    isTypeI = ~(isT & xnum>hf);
    if ~isTypeI
        r = epc0./epc;
        Fcc = 0.5*bf.*fcd.*(2-r).*x - 0.5*(x-hf).*fcd.*(x-hf)./(r.*x).*(bf-b);
        xnum = getx(Fcc, Fsc, Fst, Fft, d);
        if double(subs(x.*(1-r), x, xnum))>hf
            Fcc = fcd.*hf.*bf+fcd.*(x-hf).*b-0.5*fcd.*epc0.*x./epc.*b;
            xnum = getx(Fcc, Fsc, Fst, Fft, d);
            M = moment_with_x(4, beam_model, xnum, epc, ef, esc, es);
        else
            M = moment_with_x(3, beam_model, xnum, epc, ef, esc, es);
        end 
    end
return
end

function xnum = getx(Fcc, Fsc, Fst, Fft, xsup)
    [N, D] = numden(Fcc+Fsc-(Fst+Fft));
    p = coeffs(N);
    p = fliplr(double(p));
    x_candid = roots(p);
    indx = imag(x_candid)==0 & real(x_candid)<xsup & real(x_candid)>0;
    if sum(indx)==0
        xnum=inf;
    else
        xnum = min(x_candid(indx));
    end
return
end

function M = moment_with_x(conc_blk, beam_model, xnum, epc, ef, esc, es)
syms x
Efd = beam_model.Efd;
fcd = beam_model.fcd;
epc0 = beam_model.epc0;
epcu = beam_model.epcu;
fyd = beam_model.fyd;
fycd = beam_model.fycd;
epy = beam_model.epy; epyc = beam_model.epyc;
isT = beam_model.isT;
d = beam_model.d; a = beam_model.a; h = beam_model.h; hf = beam_model.hf;
bf = beam_model.bf; b=beam_model.b;
Esc = beam_model.Esc; Asc = beam_model.Asc;
Es = beam_model.Es; As= beam_model.As; Af = beam_model.Af;
    
epc = double(subs(epc, x, xnum));
ef = double(subs(ef, x, xnum));
esc = double(subs(esc, x, xnum));
es = double(subs(es, x, xnum));
x = xnum;

switch conc_blk
    case 1
        dn = 1.0/3*x;
        M = 1e-6*(As.*(es.*Es).*(d-dn)+Af.*ef.*Efd.*(h-dn)-Asc.*(Esc.*esc).*(a-dn));
    case 2
        Atrapzoid = fcd.*x-0.5*fcd.*epc0./epc.*x;
        wtriangle = epc0./epc.*x;
        dn = (fcd.*x.*0.5.*x-0.5*fcd.*wtriangle.*(x-1/3*wtriangle))./Atrapzoid;
        M = 1e-6*(As.*(es.*Es).*(d-dn)+Af.*ef.*Efd.*(h-dn)-Asc.*(Esc.*esc).*(a-dn));
    case 3
        wtriangle1 = epc0./epc.*x+hf-x;
        Atriangle1 = 0.5*wtriangle1.^2.*(fcd./(epc0./epc.*x));
        wtriangle2 = (x-hf);
        Atriangle2 = 0.5*wtriangle2.^2.*(fcd./(epc0./epc.*x));
        dn = (fcd.*hf.*0.5.*hf+Atriangle2.*(b./bf).*(hf+1/3.*wtriangle2)-...
            Atriangle1.*(hf-1/3*wtriangle1))./ (fcd.*hf+Atriangle2.*(b./bf)-Atriangle1);
        M = 1e-6*(As.*(es.*Es).*(d-dn)+Af.*ef.*Efd.*(h-dn)-Asc.*(Esc.*esc).*(a-dn));     
    case 4
        wtriangle = epc0./epc.*x;
        Atriangle = 0.5*wtriangle.*fcd;
        wrec = x-hf-wtriangle;
        Arec = wrec.*fcd;
        dn = (fcd.*hf.*0.5.*hf+Arec.*(b./bf).*(hf+0.5.*wrec)+...
            Atriangle.*(b./bf).*(x-2/3*wtriangle))./ ...
            (fcd.*hf+Atriangle.*(b./bf)+Arec.*(b./bf));
        M = 1e-6*(As.*(es.*Es).*(d-dn)+Af.*ef.*Efd.*(h-dn)-Asc.*(Esc.*esc).*(a-dn));         
end
return
end

function bl = isofflimit(dx, beams, sfy, sfe, My, Med)
tf = [beams.tfrp]';
efd = [beams.efd]'; Efd = [beams.Efd]';
fct = [beams.fct]';
gamma = [beams.gammaConcrete]';
tm = tf.*(sfe-sfy)./dx; tsc = 7.8.*(1.1-My./Med).*fct;
cr1 = tm+tsc>(4.5*fct./gamma);
efmax = sfe./Efd;
efsc = 0.114*tsc./sqrt(Efd.*tf);
cr2 = efmax+efsc>efd;
bl = cr1 | cr2;
return
end

function [M, rf] = adjustM(nlayer, xmin, xy, efmax, efy, My, beam)
% efy = (beam.h-xy)/(beam.d-xy)*beam.epy;
sfy = beam.Efd*efy;
% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
options = optimoptions('fmincon','TolCon',0.1,'Display','off');
[xvector, fval] = fmincon(@myobj,[xy,efy],[],[],[],[],...
    [min(xy,xmin),min(efy,efmax)], [max(xy,xmin),max(efy,efmax)], @mycon, options);
    function Med = myobj(xvector)
        x = xvector(1);
        efe = xvector(2);
        es = (beam.d-x)/(beam.h-x)*efe;
        esc = (x-beam.a)/(beam.h-x)*efe;
        carray = linspace(0, x, nlayer);
        dc = x/nlayer;
        Fc = 0; Mc =0;
        for i=1:nlayer
            c = carray(i);
            epc = c/(beam.h-x)*efe;
            dFc = blk_width(c, x, beam.hf, beam.b, beam.bf)*dc*...
                conc(epc, beam.epc0, beam.epcu, beam.fcd);
            dMc = dFc*(x-c);
            Fc = Fc+dFc;
            Mc = Mc+dMc;
        end
        Fs = beam.As*steel(es, beam.epy, beam.fyd);
        Ms = Fs*beam.d;
        Fsc = beam.Asc*steel(esc, beam.epyc, beam.fycd);
        Msc = Fsc*beam.a;
        Ff = beam.Af*beam.Efd*efe;
        Mf = Ff*beam.h;
        Med = -(Ms+Mf-Mc-Msc)*1e-6;
        Fr = (Fs+Ff-Fc-Fsc)*1e-3;
%         dx = beam.ss/(1+My/Med);
%         offlimit = isofflimit(dx, beam, sfy, beam.Efd*efe, My, -Med);
%         if offlimit
%             Fr = Fr-Ff;
%             Med = -(Ms-Mc-Msc);
%         end
    end
    function [lb, Fr] = mycon(xvector)
        x = xvector(1);
        efe = xvector(2);
        es = (beam.d-x)/(beam.h-x)*efe;
        esc = (x-beam.a)/(beam.h-x)*efe;
        carray = linspace(0, x, nlayer);
        dc = x/nlayer;
        Fc = 0; Mc =0;
        for i=1:nlayer
            c = carray(i);
            epc = c/(beam.h-x)*efe;
            dFc = blk_width(c, x, beam.hf, beam.b, beam.bf)*dc*...
                conc(epc, beam.epc0, beam.epcu, beam.fcd);
            dMc = dFc*(x-c);
            Fc = Fc+dFc;
            Mc = Mc+dMc;
        end
        Fs = beam.As*steel(es, beam.epy, beam.fyd);
        Ms = Fs*beam.d;
        Fsc = beam.Asc*steel(esc, beam.epyc, beam.fycd);
        Msc = Fsc*beam.a;
        Ff = beam.Af*beam.Efd*efe;
        Mf = Ff*beam.h;
        Med = (Ms+Mf-Mc-Msc)*1e-6;
        Fr = (Fs+Ff-Fc-Fsc)*1e-3;
        dx = beam.ss/(1+My/Med);
        offlimit = isofflimit(dx, beam, sfy, beam.Efd*efe, My, Med);
        if offlimit
            Fr = Fr-Ff;
            Med = Ms-Mc-Msc;
        end
        lb=My-Med;
    end
    function w = blk_width(c, x, hf, b, bf)
        if c>(x-hf)
            w=bf;
        else
            w=b;
        end
    end
    function sc = conc(ec, ec0, ecu, fcd)
        if ec<ec0
            sc = fcd./ec0.*ec;
        elseif ec<=ecu
            sc = fcd;
        else
            sc = 0;
        end
    end
    function sst = steel(es, epy, fyd)
        if abs(es)<epy
            sst = es.*fyd./epy;
        else
            sst = sign(es)*fyd;
        end
    end
    function sf = frp(ef, Efd, isofflimit)
        if isofflimit
            sf = 0;
        else
            sf = ef.*Efd;
        end
    end
M = -fval;
Mmax = -myobj([xmin,efmax]);
[~,rf] = mycon(xvector);
return
end