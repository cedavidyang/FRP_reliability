clear;
%% user input and data processing
running_type = input(strcat('type of running: 1 for shear+side, 2 for shear+U,',...
'3 for shear+W, 4 for flexure+IC, 5 for flexure+rupture\n'));
fc_type = 1; fccov = 0.20;
get_bias = @(x) 1/(1-1.645*x);
DESIGN_CODE = lower(input('design guidelines/code\n', 's'));

switch running_type
    case 1
        run('../user_input.m')
        SUB_TEST_DATABASE_NAME ='shear+side';
    case 2
        run('../user_input.m')
        SUB_TEST_DATABASE_NAME ='shear+U';
    case 3
        run('../user_input.m')
        SUB_TEST_DATABASE_NAME ='shear+W';
    case 4
        run('../user_input.m')
        SUB_TEST_DATABASE_NAME ='flexure+IC';
    case 5
        run('../user_input.m')
        SUB_TEST_DATABASE_NAME ='flexure+rupture';
    otherwise
        fprintf('illegal running type');
        break;
end

switch lower(DESIGN_CODE)
    case {'aci', 'acinew'}
        psi_f = 1.00;
        FACTOR_FRP = (0.50:0.05:1.00)';
        TARGET_INDEX = 3.5;
    case {'hk', 'tr', 'fib'}
        psi_f = 0.5:0.05:1.00;
        FACTOR_FRP = (1.00:0.05:2.00)';
        TARGET_INDEX = 3.8;
    case {'gb', 'gbnew'}
        psi_f = 0.5:0.05:1.00;
        FACTOR_FRP = (1.00:0.05:2.00)';
        TARGET_INDEX = 3.2;
    case {'wu'}
        psi_f = 0.5:0.05:1.00;
        FACTOR_FRP = (1.00:0.05:2.00)';
        TARGET_INDEX = 3.2;        
end

FC_COV = fccov;
FC_BIAS = get_bias(FC_COV);
FCT_COV = fccov;
FCT_BIAS = get_bias(FCT_COV);
run('../preprocessing.m')

switch running_type
    case {1, 2, 3}
        load(strcat('me_', DESIGN_CODE, '_', 'shear_fccov0.1.mat'));
        indx = FRP_FORM_TEST_ARRAY==FRP_FORM_DESIGN;
    case {4,5}
        load(strcat('me_', DESIGN_CODE, '_', 'flexure_fccov0.1.mat'));
        indx = (FAIL_MODE_TEST_ARRAY == running_type-3);
%         if strcmpi(DESIGN_CODE, 'tr')
%             indx = (FAIL_MODE_TEST_ARRAY == running_type-3)&(~failModeFromPrediction);
%         else
%             indx = (FAIL_MODE_TEST_ARRAY == running_type-3);
%         end
    otherwise
        fprintf('illegal running type');
        break;
end

%% linear plot
figure1 = figure;
axes1 = axes('Parent',figure1,'FontSize',8,'FontName','Times New Roman');
box(axes1,'on');
hold(axes1,'all');

h1 = plot(axes1, resistanceFromPrediction(indx), resistanceFromTest(indx), 'b.');
xlim(axes1, [0, 1200]);
h = refline(1, 0); set(h, 'Color','r', 'LineStyle', '-');
set(h, 'Linewidth', 1.0);
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Lower = [-Inf 0];
opts.Upper = [Inf 0];
[cf, stats] = fit(resistanceFromPrediction(indx), resistanceFromTest(indx), 'poly1', opts);
h2 = refline(cf.p1, cf.p2); set(h2, 'Color', [0, 0.5, 0], 'LineStyle', '--');
switch running_type
    case {1, 2, 3}
        xtxt1 = xlabel('Predicted shear capacity, V_{pre} (kN)','FontSize',8,'FontName','Times New Roman', 'interpreter', 'tex');
        ytxt1 = ylabel('Test result, V_{exp} (kN)','FontSize',8,'FontName','Times New Roman', 'interpreter', 'tex');
    case {4,5}
        xtxt1 = xlabel('Predicted flexural capacity, M_{pre} (kNm)','FontSize',8,'FontName','Times New Roman', 'interpreter', 'tex');
        ytxt1 = ylabel('Test result, M_{exp} (kNm)','FontSize',8,'FontName','Times New Roman', 'interpreter', 'tex');
    otherwise
        fprintf('illegal running type');
        break;
end
ylim(axes1, [0, 1200]);
lgd = legend([h1, h, h2], {'Data', 'R_{exp} = R_{pre}', 'Linear regression'}, 'Location', 'NorthWest');

%% PDF plot and distribution fitting
figure2 = figure;
axes2 = axes('Parent',figure2,'FontSize',8,'FontName','Times New Roman');
box(axes2,'on');
hold(axes2,'all');

modelerror = resistanceFromTest(indx)./resistanceFromPrediction(indx);
fprintf(strcat('model error mean of ', SUB_TEST_DATABASE_NAME, ': %.5f\n'), mean(modelerror))
fprintf(strcat('model error std of ', SUB_TEST_DATABASE_NAME, ': %.5f\n'), std(modelerror))
fprintf(strcat('model error cov of ', SUB_TEST_DATABASE_NAME, ': %.5f\n'), std(modelerror)/mean(modelerror))


legh_ = []; legt_ = {};   % handles and text for legend
[F_,X_] = ecdf(modelerror,'Function','cdf');  % compute empirical cdf
Bin_.rule = 1;
[C_,E_] = dfswitchyard('dfhistbins',modelerror,[],[],Bin_,F_,X_);
[N_,C_] = ecdfhist(F_,X_,'edges',E_); % empirical pdf from cdf
h_ = bar(axes2, C_,N_,'hist');
set(h_,'FaceColor','none','EdgeColor',[0.75, 0, 0.75],'LineStyle','-', 'LineWidth',0.5);
xtxt2 = xlabel('Model error', 'FontSize',8,'FontName','Times New Roman', 'interpreter', 'tex');
ytxt2 = ylabel('Probablity density', 'FontSize',8,'FontName','Times New Roman', 'interpreter', 'tex');
legh_(end+1) = h_; legt_{end+1} = 'Model error';
% Nudge axis limits beyond data limits
xlim_ = get(axes2,'XLim');
if all(isfinite(xlim_))
   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
   set(axes2,'XLim',xlim_)
end
x_ = linspace(xlim_(1),xlim_(2),100);
t_ = ~isnan(modelerror);
Data_ = modelerror(t_);

% --- Create fit "Normal"
fprintf('------- Create fit "Normal" -------\n')
[param_norm(1), param_norm(2)] = normfit(Data_, 0.05);
y_ = normpdf(x_,param_norm(1), param_norm(2));
h_ = plot(x_,y_,'Color', [0, 0.5, 0], 'LineStyle','-.', 'LineWidth',1.5, 'Marker','none', 'MarkerSize',4);
legh_(end+1) = h_;
legt_{end+1} = 'Normal';

% --- Create fit "Lognormal"
fprintf('------- Create fit "Lognormal" -------\n')
param_log = lognfit(Data_, 0.05);
y_ = lognpdf(x_,param_log(1),param_log(2));
h_ = plot(x_,y_,'Color','b', 'LineStyle','-', 'LineWidth',1.5,'Marker','none', 'MarkerSize',4);
legh_(end+1) = h_;
legt_{end+1} = 'Lognormal';

% --- Create fit "Gumbel(max)"
fprintf('------- Create fit "Gumbel(max)" -------\n')
param_gbl = evfit(-Data_, 0.05); % since MATLAB extreme value distribution is for minimum, use negative
y_ = evpdf(-x_,param_gbl(1), param_gbl(2));
h_ = plot(x_,y_,'Color','r','LineStyle','--', 'LineWidth',1.5,'Marker','none', 'MarkerSize',4);
legh_(end+1) = h_;
legt_{end+1} = 'Gumbel';

leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'}; 
h_ = legend(axes2,legh_,legt_,leginfo_{:}, 'FontSize',8,'FontName','Times New Roman', 'Interpreter', 'tex');  % create legend

% output fitting results and hypothesis testing
test_cdf_norm = [Data_, normcdf(Data_, param_norm(1), param_norm(2))];
test_cdf_logn = [Data_, logncdf(Data_, param_log(1), param_log(2))];
test_cdf_gbl = [-Data_, evcdf(-Data_, param_gbl(1), param_gbl(2))];
[h_norm, p_norm] = kstest(Data_, 'CDF', test_cdf_norm);
[h_logn, p_logn] = kstest(Data_, 'CDF', test_cdf_logn);
[h_gbl, p_gbl] = kstest(-Data_, 'CDF', test_cdf_gbl);
like_norm = sum(log(normpdf(Data_, param_norm(1), param_norm(2))));
like_logn = sum(log(lognpdf(Data_, param_log(1), param_log(2))));
like_gbl = sum(log(evpdf(-Data_, param_gbl(1), param_gbl(2))));
disp('------- Hypothesis tests -------\n')
disp({'normal', 'lognormal', 'gumbel'; h_norm, h_logn, h_gbl})
disp('------- p-values -------\n')
disp({'normal', 'lognormal', 'gumbel'; p_norm, p_logn, p_gbl})
disp('------- loglikelihood -------\n')
disp({'normal', 'lognormal', 'gumbel'; like_norm, like_logn, like_gbl})

%% change the following properties from figure inspector
set(figure1, 'Units','inches', 'Position', [0, 0, 3.5, 2.625]);
set(axes1, 'Units', 'normalized', 'Position', [0.16, 0.15, 0.75, 0.79]);
% savefig(figure1, strcat('../figures/HKandACInew/me_', code, '_', frp, '.fig'));
% print(figure1, strcat('../figures/HKandACInew/me_', code, '_', frp, '.eps'), '-depsc');
% set(xtxt1, 'Units', 'normalized', 'Position', [-0.1, 0.5, 0]);
% set(ytxt1, 'Units', 'normalized', 'Position', [0.5, -0.1, 0]);

set(figure2, 'Units','inches', 'Position', [0, 0, 3.5, 2.625]);
set(axes2, 'Units', 'normalized', 'Position', [0.16, 0.15, 0.75, 0.79]);
% print(figure2, strcat('../figures/HKandACInew/mepdf_', code, '_', frp, '.eps'), '-depsc');
% savefig(figure2, strcat('../figures/HKandACInew/mepdf_', code, '_', frp, '.fig'))
% set(xtxt2, 'Units', 'normalized', 'Position', [-0.1, 0.5, 0]);
% set(ytxt2, 'Units', 'normalized', 'Position', [0.5, -0.1, 0]);