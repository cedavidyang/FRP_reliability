clear;
MAX_YR = 50;
code = input('design code\n', 's');
switch lower(code)
    case {'hk','tr'}
        FACTOR_FRP = (1.00:0.05:2.00)';
        beta_T_50=3.8;
    case 'gb'
        FACTOR_FRP = (1.00:0.05:2.00)';
        beta_T_50=3.2;
    case{'aci', 'acinew'}
        FACTOR_FRP = (0.50:0.05:1.00)';
        beta_T_50=3.5;
    otherwise
end
fc_type = input('fc type: 1 for fcdet, 2 for 125 bias, 3 for cov10\n');
switch fc_type
    case 1
        fc = 'fcdet'; fccov=0.2;
    case 2
        fc = 'fc125'; fccov=0.2;
    case 3
        fc = 'cov10'; fccov=0.1;
    otherwise
end
struct_s = load(strcat('data_', code, '_shear+side_', fc, '_', num2str(fccov), '.mat'));
struct_u = load(strcat('data_', code, '_shear+U_', fc, '_', num2str(fccov), '.mat'));
struct_w = load(strcat('data_', code, '_shear+W_', fc, '_', num2str(fccov), '.mat'));
life_s = 50*log(normcdf(beta_T_50))./log(normcdf(struct_s.mean_RE)); life_s(life_s>MAX_YR) = MAX_YR;
life_u = 50*log(normcdf(beta_T_50))./log(normcdf(struct_u.mean_RE)); life_u(life_u>MAX_YR) = MAX_YR;
life_w = 50*log(normcdf(beta_T_50))./log(normcdf(struct_w.mean_RE)); life_w(life_w>MAX_YR) = MAX_YR;

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