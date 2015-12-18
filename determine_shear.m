function [resistMean, resistStd, resistSmp] = determine_shear(iDesignCase, nSim)
% Monte Carlo simulation to determine statistical properties of resistance
RO_CYLINDE_2_CUBE = 0.8;
% load data
load tmpdata

%% generate random properties
hBeamMean = H_DESIGN_ARRAY_MM(iDesignCase) - H_NORM_MEAN;
hBeamStd = H_STD_MM;
hSmp = normrnd(hBeamMean, hBeamStd, nSim, 1);

betaMean = BETA_DESIGN_ARRAY_DEG(iDesignCase)*BETA_BIAS;
betaStd = BETA_STD_DEG;
betaSmp = normrnd(betaMean, betaStd, nSim, 1);

fFrpMean = F_FRP_DESIGN_ARRAY_MPA(iDesignCase) * F_FRP_BIAS;
fFrpStd = fFrpMean * F_FRP_COV;
wblparam = fsolve(@(x) [x(1)*gamma(1+1./x(2)) - fFrpMean ; x(1).^2 * (gamma(1+2./x(2)) - ( gamma(1+1./x(2)).^2)) - fFrpStd^2],[fFrpMean;1.2/(fFrpStd/fFrpMean)], optimset('Display','off'));
fFrpSmp = wblrnd(wblparam(1), wblparam(2), nSim, 1);

fcMean = FC_DESIGN_ARRAY_MPA(iDesignCase) * FC_BIAS;
fcStd = fcMean * FC_COV;
fcSmp = normrnd(fcMean, fcStd, nSim, 1);
fcSmp(fcSmp<0) = 0;

switch DESIGN_CODE
    case {'ACI', 'aci'}        
        sqrtFcMean = sqrt(fcMean);
        sqrrFcStd = sqrtFcMean * FCT_COV;
        sqrtFcSmp = normrnd(sqrtFcMean, sqrrFcStd, nSim, 1);
    case {'HK', 'hk'}
        sqrtFcMean = sqrt(fcMean);
        sqrrFcStd = sqrtFcMean * FCT_COV;
        sqrtFcSmp = normrnd(sqrtFcMean, sqrrFcStd, nSim, 1);
    case {'GB', 'gb', 'GBnew', 'gbnew'}
        fcuMean = fcMean / RO_CYLINDE_2_CUBE;
        fctMean = 0.395 .* fcuMean.^0.55;
        fctStd = fctMean * FCT_COV;
        fctSmpNoBrittleFactor = normrnd(fctMean, fctStd, nSim, 1);      
    case {'TR', 'tr'}
        fckMean = fcMean-8; fckMean(fckMean<0) = 0;
        fctMean = 0.3 .* fckMean.^(2/3.0);
        fctStd = fctMean * FCT_COV;
        fctSmp = normrnd(fctMean, fctStd, nSim, 1); 
        sqrtFctMean = sqrt(fctMean);
        sqrtFctStd = sqrtFctMean * FCT_COV;
        sqrtFctSmp = normrnd(sqrtFctMean, sqrtFctStd, nSim, 1);        
%         sqrtFcMean = sqrt(fcMean);
%         sqrrFcStd = sqrtFcMean * FCT_COV;
%         sqrtFcSmp = normrnd(sqrtFcMean, sqrrFcStd, nSim, 1);
%         cbrtFcMean = fcMean.^(1/3);
%         cbrtFcStd = cbrtFcMean * FCT_COV;
%         cbrtFcSmp = normrnd(cbrtFcMean, cbrtFcStd, nSim, 1);
    case {'fib', 'FIB'}
        fckMean = fcMean-8; fckMean(fckMean<0) = 0;
        fctMean = 0.3 .* fckMean.^(2/3.0);
        fctStd = fctMean * FCT_COV;
        fctSmp = normrnd(fctMean, fctStd, nSim, 1); 
        sqrtFctMean = sqrt(fctMean);
        sqrtFctStd = sqrtFctMean * FCT_COV;
        sqrtFctSmp = normrnd(sqrtFctMean, sqrtFctStd, nSim, 1);        
    case {'ACInew', 'acinew'}        
        sqrtFcMean = sqrt(fcMean);
        sqrtFcStd = sqrtFcMean * FCT_COV;
        sqrtFcSmp = normrnd(sqrtFcMean, sqrtFcStd, nSim, 1);
    otherwise
end

fsMean = FS_DESIGN_ARRAY_MPA(iDesignCase) * FS_BIAS;
% fsStd = fsMean * FS_COV;
fsLogStd = sqrt( log( FS_COV^2 + 1 ) );
fsLogMean = log( fsMean ) - .5*fsLogStd.^2;
fsSmp = lognrnd( fsLogMean, fsLogStd, nSim, 1);

switch DESIGN_CODE
    case {'aci' 'ACI'}
        resistSmp = shear_total_ACI_MC(iDesignCase, hSmp, betaSmp, fFrpSmp, fcSmp, sqrtFcSmp, fsSmp);
    case {'hk' 'HK'}
        resistSmp = shear_total_HK_MC(iDesignCase, hSmp, betaSmp, fFrpSmp, fcSmp, sqrtFcSmp, fsSmp);
    case{'GB' 'gb'}
        resistSmp = shear_total_GB_MC(iDesignCase, hSmp, betaSmp, fFrpSmp, fcSmp, fctSmpNoBrittleFactor, fsSmp);
    case{'GBnew' 'gbnew'}
        resistSmp = shear_total_GBnew_MC(iDesignCase, hSmp, betaSmp, fFrpSmp, fcSmp, fctSmpNoBrittleFactor, fsSmp);    
    case{'TR' 'tr'}
        resistSmp = shear_total_TR_MC(iDesignCase, hSmp, betaSmp, fFrpSmp, fcSmp, fctSmp, sqrtFctSmp, fsSmp);
    case{'FIB' 'fib'}
        resistSmp = shear_total_fib_MC(iDesignCase, hSmp, betaSmp, fFrpSmp, fcSmp, fctSmp, sqrtFctSmp, fsSmp);         
    case {'acinew' 'ACInew'}
        resistSmp = shear_total_ACInew_MC(iDesignCase, hSmp, betaSmp, fFrpSmp, fcSmp, sqrtFcSmp, fsSmp);       
    otherwise
end

resistMean = mean(resistSmp);
resistStd = std(resistSmp);

return
end