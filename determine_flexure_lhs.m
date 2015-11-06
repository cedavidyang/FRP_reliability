function [resistMean, resistStd, resistSmp] = determine_flexure_lhs(iDesignCase, nSim)
% Monte Carlo simulation to determine statistical properties of resistance
RO_CYLINDE_2_CUBE = 0.8;

% load data
load tmpdata

nRv = 9;
cdfs = lhsdesign(nSim, nRv);

%% generate random properties with Rosenblueth's 2K+1 point estimate method
hBeamMean = H_DESIGN_ARRAY_MM(iDesignCase) - H_NORM_MEAN;
hBeamStd = H_STD_MM;
rv1 = makedist('norm', 'mu', hBeamMean, 'sigma', hBeamStd);
hSmp = rv1.icdf(cdfs(:,1));

EFrpMean = E_FRP_DESIGN_ARRAY_MPA(iDesignCase) * E_FRP_BIAS;
EFrpStd = EFrpMean * E_FRP_COV;
EFrpLogStd = sqrt( log( E_FRP_COV^2 + 1 ) );
EFrpLogMean = log( EFrpMean ) - .5*EFrpLogStd.^2;
rv2 = makedist('lognormal', 'mu', EFrpLogMean, 'sigma', EFrpLogStd);
EFrpSmp = rv2.icdf(cdfs(:,2));

fFrpMean = F_FRP_DESIGN_ARRAY_MPA(iDesignCase) * F_FRP_BIAS;
fFrpStd = fFrpMean * F_FRP_COV;
wblparam = fsolve(@(x) [x(1)*gamma(1+1./x(2)) - fFrpMean ; x(1).^2 * (gamma(1+2./x(2)) - ( gamma(1+1./x(2)).^2)) - fFrpStd^2],[fFrpMean;1.2/(fFrpStd/fFrpMean)], optimset('Display','off'));
rv3 = makedist('weibull', 'a', wblparam(1), 'b', wblparam(2));
fFrpSmp = rv3.icdf(cdfs(:,3));

fcMean = FC_DESIGN_ARRAY_MPA(iDesignCase) * FC_BIAS;
fcStd = fcMean * FC_COV;
rv4 = makedist('norm', 'mu', fcMean, 'sigma', fcStd);
fcSmp = rv4.icdf(cdfs(:,4)); fcSmp(fcSmp<0) = 0;

switch DESIGN_CODE
    case {'ACI', 'aci'}        
        sqrtFcMean = sqrt(fcMean);
        sqrtFcStd = sqrtFcMean * FCT_COV;
        rv5 = makedist('norm', 'mu', sqrtFcMean, 'sigma', sqrtFcStd);
        sqrtFcSmp = rv5.icdf(cdfs(:,5));
    case {'HK', 'hk'}
        fctMean = 0.5*sqrt(fcMean);
        fctStd = fctMean * FCT_COV;
        rv6 = makedist('norm', 'mu', fctMean, 'sigma', fctStd);
        fctSmp = rv6.icdf(cdfs(:,6));
    case {'GB', 'gb'}
        fcuMean = fcMean / RO_CYLINDE_2_CUBE;
        fctMean = 0.395 .* fcuMean.^0.55;
        fctStd = fctMean * FCT_COV;
        fctMeanNoBrittleFactor = fctMean;
        fctStdNoBrittleFactor = fctStd;
        rv7 = makedist('norm', 'mu', fctMeanNoBrittleFactor, 'sigma', fctStdNoBrittleFactor);
        fctSmpNoBrittleFactor = rv7.icdf(cdfs(:,7));   
    otherwise
end

fsMean = FS_DESIGN_ARRAY_MPA(iDesignCase) * FS_BIAS;
fsStd = fsMean * FS_COV;
fsLogStd = sqrt( log( FS_COV^2 + 1 ) );
fsLogMean = log( fsMean ) - .5*fsLogStd.^2;
rv8 = makedist('lognormal', 'mu', fsLogMean, 'sigma', fsLogStd);
fsSmp = rv8.icdf(cdfs(:,8));

areaSteelMean = AREA_STEEL_DESIGN_ARRAY_MM2(iDesignCase) * AS_BIAS;
areaSteelStd = areaSteelMean * AS_COV;
rv9 = makedist('norm', 'mu', areaSteelMean, 'sigma', areaSteelStd);
areaSteelSmp = rv9.icdf(cdfs(:,9));

switch DESIGN_CODE
    case {'aci' 'ACI'}
        resistSmp = flexure_total_ACI_MC(iDesignCase, hSmp, fcSmp, sqrtFcSmp, EFrpSmp, fFrpSmp, fsSmp, areaSteelSmp);
    case {'hk' 'HK'}
        resistSmp = flexure_total_HK_MC(iDesignCase, hSmp, fcSmp, fctSmp, EFrpSmp, fFrpSmp, fsSmp, areaSteelSmp);      
    case{'GB' 'gb'}
        resistSmp = flexure_total_GB_MC(iDesignCase, hSmp, fcSmp, fctSmpNoBrittleFactor, EFrpSmp, fFrpSmp, fsSmp, areaSteelSmp);   
    otherwise
end

resistMean = mean(resistSmp);
resistStd = std(resistSmp);

return
end