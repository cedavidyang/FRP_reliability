function [] = model_error_analysis(Mpre, Mexp)

ModelError = Mexp ./ Mpre;
figure;
plot(Mpre, Mexp, 'o', 'markerface', 'b')
xlim([0, 2000]); ylim([0, 2000]);
[cf, stats] = fit(Mpre, Mexp, 'poly1');
h2 = refline(cf.p1, cf.p2);
set(h2, 'color', 'g', 'linestyle','-.', 'linewidth', 2)
h1 = refline(1, 0);
set(h1, 'color', 'r', 'linewidth', 2)
xlim([0, 1400]); ylim([0, 1400]);
xlabel('M_{prediction} (kNm)')
ylabel('M_{experiment} (kNm)')
annotation(gcf, 'textbox', [.1 .1 .3 .3], 'String', sprintf('y = %.4f x + %.4f\n R^2 = %.4f', cf.p1, cf.p2, stats.rsquare)); 

figure; 
[param_wbl, param_log, param_gam] = dfit1(ModelError);
% Chi-square test
nbins = 10;
alpha = 0.05;
nparams = 2;
% [h_wbl p_wbl] = chi2test(ModelError, nbins, alpha, param_wbl, nparams, 'Weibull');
% [h_log p_log] = chi2test(ModelError, nbins, alpha, param_log, nparams, 'Lognormal');
% [h_gam p_gam] = chi2test(ModelError, nbins, alpha, param_gam, nparams, 'Gamma');
[h_wbl, p_wbl] = chi2gof(ModelError, 'cdf', fitdist(ModelError, 'Weibull'), 'alpha', alpha);
[h_log, p_log] = chi2gof(ModelError, 'cdf', fitdist(ModelError, 'Lognormal'), 'alpha', alpha);
[h_gam, p_gam] = chi2gof(ModelError, 'cdf', fitdist(ModelError, 'Gamma'), 'alpha', alpha);
% Kolmogorov-Smirnov One-Sample Test
% h_wbl_ks = Ukstest(ModelError, alpha, param_wbl, 'Weibull');
% h_log_ks = Ukstest(ModelError, alpha, param_log, 'Lognormal');
% h_gam_ks = Ukstest(ModelError, alpha, param_gam, 'Gamma');
h_wbl_ks = kstest(ModelError, fitdist(ModelError, 'wbl'), alpha, 'unequal');
h_log_ks = kstest(ModelError, fitdist(ModelError, 'logn'), alpha, 'unequal');
h_gam_ks = kstest(ModelError, fitdist(ModelError, 'gam'), alpha, 'unequal');

xlabel(sprintf('Model error (\\mu = %.4f, COV = %.2f%%)', mean(ModelError), std(ModelError, 1)./mean(ModelError)*100));

annotation(gcf,'textbox', [.1 .1 .3 .3], 'interpreter', 'tex', 'String',{['chi-square test (\alpha=5%)', sprintf('\n'), sprintf('h_{wbl}=%d, p_{wbl}=%6.4f\n', h_wbl, p_wbl), ...
                                               sprintf('h_{log}=%d, p_{log}=%6.4f\n', h_log, p_log), ...
                                               sprintf('h_{gam}=%d, p_{gam}=%6.4f\n', h_gam, p_gam), ...
                                               'Kolmogorov-Smirnov test(\alpha=5%)', sprintf('\n'), sprintf('h_{wbl}=%d, h_{log}=%d, h_{gam}=%d', h_wbl_ks, h_log_ks, h_gam_ks)]}, 'LineStyle','none');
return
end