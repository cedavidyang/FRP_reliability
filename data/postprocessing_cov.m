clear;
MAX_YR = 50;
fccov=0.10:0.01:0.20;
psif_s = 0.70; psif_u = 1.0; psif_w = 1.0;
code = input('design code\n', 's');
switch lower(code)
    case {'hk','tr'}
        psi_f = (0.1:0.05:2.0)';
        FACTOR_FRP = (1.00:0.05:2.00)';
        beta_T_50=3.8;
    case 'gb'
        psi_f = (0.1:0.05:2.0)';
        FACTOR_FRP = 1.00;
        beta_T_50=3.2;
    case{'aci', 'acinew'}
        psi_f = 1.00;
        FACTOR_FRP = (0.1:0.05:2.0)';
        beta_T_50=3.5;
    otherwise
end
ncov = length(fccov);
re_s = zeros(ncov,1);
re_u = zeros(ncov,1);
re_w = zeros(ncov,1);
for icov=1:ncov 
    cov = fccov(icov);
    struct_s = load(strcat(code, '_', 'shear+side', '_fccov.mat'));
    struct_u = load(strcat(code, '_', 'shear+U', '_fccov.mat'));
    struct_w = load(strcat(code, '_', 'shear+W', '_fccov.mat'));
    re_s_smp = struct_s.re_data{icov}(:,FACTOR_FRP==psif_side,:); re_s_smp = re_s_smp(:);
    re_s(icov) = mean(re_s_smp);
    re_u_smp = struct_u.re_data{icov}(:,FACTOR_FRP==psif_u,:); re_u_smp = re_u_smp(:);
    re_u(icov) = mean(re_u_smp);
    re_w_smp = struct_w.re_data{icov}(:,FACTOR_FRP==psif_w,:); re_w_smp = re_w_smp(:);
    re_w(icov) = mean(re_w_smp);
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