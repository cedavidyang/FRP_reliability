clear;
code = 'acinew';
% code = 'hk';

load(strcat('me_', code));
run('../user_input.m')
run('../preprocessing.m')
scheme=3;
switch scheme
    case 1
        frp = 'side';
    case 2
        frp = 'U';
    case 3
        frp = 'W';
    otherwise
end
indx = FRP_FORM_TEST_ARRAY==scheme;

% post processing of model error analysis
figure1 = figure;
axes1 = axes('Parent',figure1,'FontSize',8,'FontName','Times New Roman');
box(axes1,'on');
hold(axes1,'all');

h1 = plot(axes1, resistanceFromPrediction(indx), resistanceFromTest(indx), 'k.');
xlim(axes1, [0, 1200]);
h = refline(1, 0); set(h, 'Color','k', 'LineStyle', '-');
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Lower = [-Inf 0];
opts.Upper = [Inf 0];
[cf, stats] = fit(resistanceFromPrediction(indx), resistanceFromTest(indx), 'poly1', opts);
h = refline(cf.p1, cf.p2); set(h, 'Color','k', 'LineStyle', '--');
xtxt1 = xlabel('Predicted shear capacity, $V_{pre}$ (kN)','FontSize',8,'FontName','Times New Roman', 'interpreter', 'latex');
ytxt1 = ylabel('Test result $V_{exp}$ (kN)','FontSize',8,'FontName','Times New Roman', 'interpreter', 'latex');
ylim(axes1, [0, 1200]);

figure2 = figure;
axes2 = axes('Parent',figure2,'FontSize',8,'FontName','Times New Roman');
box(axes2,'on');
hold(axes2,'all');

modelerror = resistanceFromTest(indx)./resistanceFromPrediction(indx);

legh_ = []; legt_ = {};   % handles and text for legend
[F_,X_] = ecdf(modelerror,'Function','cdf');  % compute empirical cdf
Bin_.rule = 1;
[C_,E_] = dfswitchyard('dfhistbins',modelerror,[],[],Bin_,F_,X_);
[N_,C_] = ecdfhist(F_,X_,'edges',E_); % empirical pdf from cdf
h_ = bar(axes2, C_,N_,'hist');
set(h_,'FaceColor','none','EdgeColor','k','LineStyle','-', 'LineWidth',0.5);
xtxt2 = xlabel('Model error', 'FontSize',8,'FontName','Times New Roman', 'interpreter', 'latex');
ytxt2 = ylabel('Probablity density', 'FontSize',8,'FontName','Times New Roman', 'interpreter', 'latex');
legh_(end+1) = h_; legt_{end+1} = 'Model error';
% Nudge axis limits beyond data limits
xlim_ = get(axes2,'XLim');
if all(isfinite(xlim_))
   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
   set(axes2,'XLim',xlim_)
end
x_ = linspace(xlim_(1),xlim_(2),100);
% --- Create fit "Normal"
t_ = ~isnan(modelerror);
Data_ = modelerror(t_);
[p_(1), p_(2)] = normfit(Data_, 0.05);
y_ = normpdf(x_,p_(1), p_(2));
h_ = plot(x_,y_,'Color', 'k', 'LineStyle','-.', 'LineWidth',1, 'Marker','none', 'MarkerSize',4);
legh_(end+1) = h_;
legt_{end+1} = 'Normal';
param_norm = p_;
% --- Create fit "Lognormal"
t_ = ~isnan(modelerror);
Data_ = modelerror(t_);
p_ = lognfit(Data_, 0.05);
y_ = lognpdf(x_,p_(1), p_(2));
h_ = plot(x_,y_,'Color','k', 'LineStyle','-', 'LineWidth',1,'Marker','none', 'MarkerSize',4);
legh_(end+1) = h_;
legt_{end+1} = 'Lognormal';
param_log = p_;
% --- Create fit "Gamma"
t_ = ~isnan(modelerror);
Data_ = modelerror(t_);
p_ = gamfit(Data_, 0.05);
y_ = gampdf(x_,p_(1), p_(2));
h_ = plot(x_,y_,'Color','k','LineStyle','--', 'LineWidth',1,'Marker','none', 'MarkerSize',4);
legh_(end+1) = h_;
legt_{end+1} = 'Gamma';
param_gam = p_;

leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'}; 
h_ = legend(axes2,legh_,legt_,leginfo_{:}, 'FontSize',8,'FontName','Times New Roman', 'Interpreter', 'latex');  % create legend

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