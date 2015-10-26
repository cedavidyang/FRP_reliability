clear;
MAX_YR = 50;
fccov=0.10:0.01:0.20;
psif_s = 0.70; psif_u = 1.0; psif_w = 1.0;
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
% scheme_num = input('strengthening scheme: 1 shear+side, 2 shear+U, 3 shear+W\n', 's');
% switch scheme
%     case 1
%         scheme = 'shear+side';
%     case 2
%         scheme = 'shear+U';
%     case 3
%         scheme = 'shear+W';
%     otherwise
% end   
ncov = length(fccov);
re_s = zeros(ncov,1);
re_u = zeros(ncov,1);
re_w = zeros(ncov,1);
for icov=1:ncov 
    cov = fccov(icov);
    struct_s = load(strcat('data_', code, '_', 'shear+side', '_cov', num2str(100*cov), '_', num2str(cov), '.mat'));
    struct_u = load(strcat('data_', code, '_', 'shear+U', '_cov', num2str(100*cov), '_', num2str(cov), '.mat'));
    struct_w = load(strcat('data_', code, '_', 'shear+W', '_cov', num2str(100*cov), '_', num2str(cov), '.mat'));
    re_s(icov) = struct_s.mean_RE(FACTOR_FRP==psif_s);
    re_u(icov) = struct_u.mean_RE(FACTOR_FRP==psif_u);
    re_w(icov) = struct_w.mean_RE(FACTOR_FRP==psif_w);
end

%% post processing of reliability analysis
figure1 = figure;
axes1 = axes('Parent',figure1,'FontSize',8,'FontName','Times New Roman');
box(axes1,'on');
hold(axes1,'all');

hs1 = plot(fccov, re_s, 'ko-', 'LineWidth', 1, 'MarkerFace', 'k', 'MarkerSize', 4);
hu1 = plot(fccov, re_u, 'k^-', 'LineWidth', 1, 'MarkerFace', 'k', 'MarkerSize', 4);
hw1 = plot(fccov, re_w, 'kv-', 'LineWidth', 1, 'MarkerFace', 'k', 'MarkerSize', 4);

xtxt1 = xlabel('COV of concrete strength','FontSize',8,...
    'FontName','Times New Roman', 'Interpreter','latex');
ytxt1 = ylabel('Mean reliability index, $\bar{\beta}$','FontSize',8,...
    'FontName','Times New Roman', 'Interpreter','latex');

%% change the following properties from figure inspector
set(figure1, 'Units','inches', 'Position', [0, 0, 3.5, 2.625]);
set(axes1, 'Units', 'normalized', 'Position', [0.16, 0.15, 0.75, 0.79]);
% savefig(figure1, strcat('../figures/HKandACInew/re_', code, '_', fc, '.fig'));
% set(xtxt1, 'Units', 'normalized', 'Position', [-0.1, 0.5, 0]);
% set(ytxt1, 'Units', 'normalized', 'Position', [0.5, -0.1, 0]);