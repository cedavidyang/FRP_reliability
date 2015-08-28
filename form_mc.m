function beta = form_mc(modelErrorMean, modelErrorStd,...
                        resistMean, resistStd,...
                        liveMean, liveStd, ...
                        deadMean, deadStd)
% determine reliablity index beta with Monte Carlo-FORM method


modelErorCov = modelErrorStd / modelErrorMean;
modelErorLogStd = sqrt( log( modelErorCov^2 + 1 ) );
modelErorLogMean = log( modelErrorMean ) - .5*modelErorLogStd.^2;

%% D and R is correlated
rho_RD = 0;
% calculate Nataf coefficient R, since D and R are normal distribution, R=1
muX = [modelErrorMean; resistMean; deadMean; liveMean];
sigmaX = [modelErrorStd; resistStd; deadStd; liveStd];
aEv = sqrt(6)*sigmaX(4)/pi; uEv = -psi(1)*aEv-muX(4);
% muX1 = muX; sigmaX1 = sigmaX;
index = 10; deltaBeta = 10;
rhoX1 = [1 0 0 0; 0 1 rho_RD 0; 0 rho_RD 1 0; 0 0 0 1];
x = muX; normX = eps;
while abs(norm(x)-normX)/normX>1e-8 || deltaBeta > 2e-3
    normX = norm(x);
    g = x(1)*x(2)-(x(3)+x(4));
    gX = [x(2); x(1); -1; -1];
    cdfX = [logncdf(x(1), modelErorLogMean, modelErorLogStd);
            normcdf(x(2), resistMean, resistStd);
            normcdf(x(3), deadMean, deadStd);
            1-evcdf(-x(4), uEv, aEv)];
    pdfX = [lognpdf(x(1), modelErorLogMean, modelErorLogStd);
            normpdf(x(2), resistMean, resistStd);
            normpdf(x(3), deadMean, deadStd);        
            evpdf(-x(4), uEv, aEv)];  
    nc = norminv(cdfX);
    sigmaX1 = normpdf(nc)./pdfX;
    muX1 = x - nc .* sigmaX1;
    gs = gX.*sigmaX1;
    alphaX = -rhoX1*gs/sqrt(gs'*rhoX1*gs);
    bbeta = (g+gX'*(muX1-x)) / sqrt(gs'*rhoX1*gs);
    % deviation from previous bbeta
    deltaBeta = abs(bbeta - index);
    % update beta
    index = bbeta;
    % update x
    x = muX1 + bbeta*sigmaX1.*alphaX;
end

beta = index;

return
end