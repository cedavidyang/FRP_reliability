function [resistMean, resistStd, resistSmp] = determine_flexure_rosenblueth(iDesignCase, nSim)
% Monte Carlo simulation to determine statistical properties of resistance
RO_CYLINDE_2_CUBE = 0.8;

% load data
load tmpdata

%% generate random properties with Rosenblueth's 2K+1 point estimate method
hBeamMean = H_DESIGN_ARRAY_MM(iDesignCase) - H_NORM_MEAN;
hBeamStd = H_STD_MM;
hSmp = normrnd(hBeamMean, hBeamStd, nSim, 1);

EFrpMean = E_FRP_DESIGN_ARRAY_MPA(iDesignCase) * E_FRP_BIAS;
EFrpStd = EFrpMean * E_FRP_COV;
EFrpLogStd = sqrt( log( E_FRP_COV^2 + 1 ) );
EFrpLogMean = log( EFrpMean ) - .5*EFrpLogStd.^2;
EFrpSmp = lognrnd( EFrpLogMean, EFrpLogStd, nSim, 1);

fFrpMean = F_FRP_DESIGN_ARRAY_MPA(iDesignCase) * F_FRP_BIAS;
fFrpStd = fFrpMean * F_FRP_COV;
wblparam = fsolve(@(x) [x(1)*gamma(1+1./x(2)) - fFrpMean ; x(1).^2 * (gamma(1+2./x(2)) - ( gamma(1+1./x(2)).^2)) - fFrpStd^2],[fFrpMean;1.2/(fFrpStd/fFrpMean)], optimset('Display','off'));
fFrpSmp = wblrnd(wblparam(1), wblparam(2), nSim, 1);

fcMean = FC_DESIGN_ARRAY_MPA(iDesignCase) * FC_BIAS;
fcStd = fcMean * FC_COV;
fcSmp = normrnd(fcMean, fcStd, nSim, 1);

switch DESIGN_CODE
    case {'ACI', 'aci'}        
        sqrtFcMean = sqrt(fcMean);
        sqrtFcStd = sqrtFcMean * FCT_COV;
        sqrtFcSmp = normrnd(sqrtFcMean, sqrtFcStd, nSim, 1);
    case {'HK', 'hk'}
        fctMean = 0.5*sqrt(fcMean);
        fctStd = fctMean * FCT_COV;
        fctSmp = normrnd(fctMean, fctStd, nSim, 1);
    case {'GB', 'gb'}
        fcuMean = fcMean / RO_CYLINDE_2_CUBE;
        fctMean = 0.395 .* fcuMean.^0.55;
        fctStd = fctMean * FCT_COV;
        fctMeanNoBrittleFactor = fctMean;
        fctStdNoBrittleFactor = fctStd;
        fctSmpNoBrittleFactor = normrnd(fctMean, fctStd, nSim, 1);        
    otherwise
end

fsMean = FS_DESIGN_ARRAY_MPA(iDesignCase) * FS_BIAS;
fsStd = fsMean * FS_COV;
fsLogStd = sqrt( log( FS_COV^2 + 1 ) );
fsLogMean = log( fsMean ) - .5*fsLogStd.^2;
fsSmp = lognrnd( fsLogMean, fsLogStd, nSim, 1);

areaSteelMean = AREA_STEEL_DESIGN_ARRAY_MM2(iDesignCase) * AS_BIAS;
areaSteelStd = areaSteelMean * AS_COV;
areaSteelSmp = normrnd(areaSteelMean, areaSteelStd, nSim, 1);

switch DESIGN_CODE
    case {'aci' 'ACI'}
        xMeanArray = [hBeamMean; fcMean; sqrtFcMean; EFrpMean; fFrpMean; fsMean; areaSteelMean];
        xStdArray = [hBeamStd; fcStd; sqrtFcStd; EFrpStd; fFrpStd; fsStd; areaSteelStd];
        resist0 = flexure_total_ACI_MC(iDesignCase, hBeamMean, fcMean, sqrtFcMean, EFrpMean, fFrpMean, fsMean, areaSteelMean);
        nV = length( xMeanArray );
        resistArray = zeros( nV, 1);
        covResistArray = zeros( nV, 1 );
        for iV = 1:nV
            xArrayPos = xMeanArray;
            xArrayNeg = xMeanArray;
            xArrayPos(iV) = xMeanArray(iV) +  xStdArray(iV);
            xArrayNeg(iV) = xMeanArray(iV) -  xStdArray(iV);
            resistPos = flexure_total_ACI_MC(iDesignCase, xArrayPos(1), xArrayPos(2), xArrayPos(3), xArrayPos(4), xArrayPos(5), xArrayPos(6), xArrayPos(7));
            resistNeg = flexure_total_ACI_MC(iDesignCase, xArrayNeg(1), xArrayNeg(2), xArrayNeg(3), xArrayNeg(4), xArrayNeg(5), xArrayNeg(6), xArrayNeg(7));
            resistArray(iV) = 1/2 * (resistPos + resistNeg);
            covResistArray(iV) = (resistPos-resistNeg) / (resistPos+resistNeg);
        end
        resistMean = resist0 * prod(resistArray./resist0);
        resistCov = sqrt( prod(1+covResistArray.^2) - 1 );
        resistStd = resistMean * resistCov;
        resistSmp = 0;
    case {'hk' 'HK'}
        xMeanArray = [hBeamMean; fcMean; fctMean; EFrpMean; fFrpMean; fsMean; areaSteelMean];
        xStdArray = [hBeamStd; fcStd; fctStd; EFrpStd; fFrpStd; fsStd; areaSteelStd];
        resist0 = flexure_total_HK_MC(iDesignCase, hBeamMean, fcMean, fctMean, EFrpMean, fFrpMean, fsMean, areaSteelMean);
        nV = length( xMeanArray );
        resistArray = zeros( nV, 1);
        covResistArray = zeros( nV, 1 );
        for iV = 1:nV
            xArrayPos = xMeanArray;
            xArrayNeg = xMeanArray;
            xArrayPos(iV) = xMeanArray(iV) +  xStdArray(iV);
            xArrayNeg(iV) = xMeanArray(iV) -  xStdArray(iV);
            resistPos = flexure_total_HK_MC(iDesignCase, xArrayPos(1), xArrayPos(2), xArrayPos(3), xArrayPos(4), xArrayPos(5), xArrayPos(6), xArrayPos(7));
            resistNeg = flexure_total_HK_MC(iDesignCase, xArrayNeg(1), xArrayNeg(2), xArrayNeg(3), xArrayNeg(4), xArrayNeg(5), xArrayNeg(6), xArrayNeg(7));
            resistArray(iV) = 1/2 * (resistPos + resistNeg);
            covResistArray(iV) = (resistPos-resistNeg) / (resistPos+resistNeg);
        end
        resistMean = resist0 * prod(resistArray./resist0);
        resistCov = sqrt( prod(1+covResistArray.^2) - 1 );
        resistStd = resistMean * resistCov;
        resistSmp = 0;        
    case{'GB' 'gb'}
        xMeanArray = [hBeamMean; fcMean; fctMeanNoBrittleFactor; EFrpMean; fFrpMean; fsMean; areaSteelMean];
        xStdArray = [hBeamStd; fcStd; fctStdNoBrittleFactor; EFrpStd; fFrpStd; fsStd; areaSteelStd];
        resist0 = flexure_total_GB_MC(iDesignCase, hBeamMean, fcMean, fctMeanNoBrittleFactor, EFrpMean, fFrpMean, fsMean, areaSteelMean);
        nV = length( xMeanArray );
        resistArray = zeros( nV, 1);
        covResistArray = zeros( nV, 1 );
        for iV = 1:nV
            xArrayPos = xMeanArray;
            xArrayNeg = xMeanArray;
            xArrayPos(iV) = xMeanArray(iV) +  xStdArray(iV);
            xArrayNeg(iV) = xMeanArray(iV) -  xStdArray(iV);
            resistPos = flexure_total_GB_MC(iDesignCase, xArrayPos(1), xArrayPos(2), xArrayPos(3), xArrayPos(4), xArrayPos(5), xArrayPos(6), xArrayPos(7));
            resistNeg = flexure_total_GB_MC(iDesignCase, xArrayNeg(1), xArrayNeg(2), xArrayNeg(3), xArrayNeg(4), xArrayNeg(5), xArrayNeg(6), xArrayNeg(7));
            resistArray(iV) = 1/2 * (resistPos + resistNeg);
            covResistArray(iV) = (resistPos-resistNeg) / (resistPos+resistNeg);
        end
        resistMean = resist0 * prod(resistArray./resist0);
        resistCov = sqrt( prod(1+covResistArray.^2) - 1 );
        resistStd = resistMean * resistCov;
        resistSmp = 0;     
    otherwise
end

return
end