clear;
FACTOR_FRP = (0.50:0.05:1.00)';
code = 'acinew';
% code = 'aci';
fc = 'fcdet'; fc='85'; %fc='427';
% fc = 'fc125';
% fc = 'cov10';
struct_s = load(strcat('data_', code, '_shear+side_', fc));
struct_u = load(strcat('data_', code, '_shear+U_', fc));
struct_w = load(strcat('data_', code, '_shear+W_', fc));
beta_T_50=3.5;
life_s = 50*log(normcdf(beta_T_50))./log(normcdf(struct_s.mean_RE)); life_s(life_s>75) = 75;
life_u = 50*log(normcdf(beta_T_50))./log(normcdf(struct_u.mean_RE)); life_u(life_u>75) = 75;
life_w = 50*log(normcdf(beta_T_50))./log(normcdf(struct_w.mean_RE)); life_w(life_w>75) = 75;

%% post processing of reliability analysis
figure1 = figure;
axes1 = axes('Parent',figure1,'FontSize',8,'FontName','Times New Roman');
box(axes1,'on');
hold(axes1,'all');

hs1 = plot(FACTOR_FRP, struct_s.mean_RE, 'ko-', 'LineWidth', 1, 'MarkerFace', 'k', 'MarkerSize', 4);
hu1 = plot(FACTOR_FRP, struct_u.mean_RE, 'k^-', 'LineWidth', 1, 'MarkerFace', 'k', 'MarkerSize', 4);
hw1 = plot(FACTOR_FRP, struct_w.mean_RE, 'kv-', 'LineWidth', 1, 'MarkerFace', 'k', 'MarkerSize', 4);

xtxt1 = xlabel('FRP strength reduction factor','FontSize',8,...
    'FontName','Times New Roman', 'Interpreter','latex');
ytxt1 = ylabel('Mean reliability index, $\bar{\beta}$','FontSize',8,...
    'FontName','Times New Roman', 'Interpreter','latex');

figure2 = figure;
axes2 = axes('Parent',figure2,'FontSize',8,'FontName','Times New Roman');
box(axes2,'on');
hold(axes2,'all');

hs2 = plot(FACTOR_FRP, life_s, 'ko-', 'LineWidth', 1, 'MarkerFace', 'k', 'MarkerSize', 4);
hu2 = plot(FACTOR_FRP, life_u, 'k^-', 'LineWidth', 1, 'MarkerFace', 'k', 'MarkerSize', 4);
hw2 = plot(FACTOR_FRP, life_w, 'kv-', 'LineWidth', 1, 'MarkerFace', 'k', 'MarkerSize', 4);

xtxt2 = xlabel('FRP strength reduction factor','FontSize',8,...
    'FontName','Times New Roman', 'Interpreter','latex');
ytxt2 = ylabel('Deisgn sevice life after strengthening (year)','FontSize',8,...
    'FontName','Times New Roman', 'Interpreter','latex');

%% change the following properties from figure inspector
set(figure1, 'Units','inches', 'Position', [0, 0, 3.5, 2.625]);
set(axes1, 'Units', 'normalized', 'Position', [0.16, 0.15, 0.75, 0.79]);
% savefig(figure1, strcat('../figures/HKandACInew/re_', code, '_', fc, '.fig'));
% set(xtxt1, 'Units', 'normalized', 'Position', [-0.1, 0.5, 0]);
% set(ytxt1, 'Units', 'normalized', 'Position', [0.5, -0.1, 0]);

set(figure2, 'Units','inches', 'Position', [0, 0, 3.5, 2.625]);
set(axes2, 'Units', 'normalized', 'Position', [0.16, 0.15, 0.75, 0.79]);
% savefig(figure2, strcat('../figures/HKandACInew/life_', code, '_', fc, '.fig'))
% set(xtxt2, 'Units', 'normalized', 'Position', [-0.1, 0.5, 0]);
% set(ytxt2, 'Units', 'normalized', 'Position', [0.5, -0.1, 0]);