clear;

code = input('design code\n', 's');
fail_mode= input('Failure mode: 1 for debonding, 2 for rupture\n');
switch lower(code)
    case {'aci', 'acinew'}
        psi_f = (0.10:0.05:2.00)';
        beta_T_50 = 3.5;
    case {'hk', 'tr'}
        psi_f = (0.10:0.05:2.00)';
        beta_T_50 = 3.8;
    case {'gb'}
        psi_f = (0.10:0.05:2.00)';
        beta_T_50 = 3.7;    
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

switch fail_mode
    case 1
        struct_s = load(strcat(code, '_', 'shear', '+side_', fc, num2str(fccov), '.mat'));
        struct_u = load(strcat(code, '_', 'shear', '+U_', fc, num2str(fccov), '.mat'));
        struct_ic = load(strcat(code, '_', 'flexure', '+IC_', fc, num2str(fccov), '.mat'));
        data_struct_array = {struct_s, struct_u, struct_ic}';
        lgd = {'Side bonding', 'U-jacketing', 'IC debonding'};
        xbound = [0.1, 1.0];
    case 2
        struct_w = load(strcat(code, '_', 'shear', '+W_', fc, num2str(fccov), '.mat'));
        struct_rup = load(strcat(code, '_', 'flexure', '+rupture_', fc, num2str(fccov), '.mat'));
        data_struct_array = {struct_w, struct_rup}';
        lgd = {'Complete wrapping', 'FRP rupture'};
        xbound = [0.5, 1.5];
end

%% debonding failure

figs = zeros(1,3);
axs = zeros(1,3);
for ifig=1:3
    figs(ifig) = figure;
    axs(ifig) = axes('Parent',figs(ifig),'FontSize',8,'FontName','Times New Roman');
    box(axs(ifig),'on');
    hold(axs(ifig),'all');
    xlim(axs(ifig),xbound);
end
colors = {'k', 'b', 'r'};
markers = {'o', '^', 'v'};

gamma = 1:0.05:2;

for idata = 1:length(data_struct_array)
    [~,~,~,~,~,norm_RE_yr30, norm_RE_yr40, norm_RE_yr50] = get_re_data( ...
        data_struct_array{idata}, psi_f, code, beta_T_50);
    [~, indx30] = min(norm_RE_yr30);
    [~, indx40] = min(norm_RE_yr40);
    [~, indx50] = min(norm_RE_yr50);
    frp_calibration_30 = psi_f(indx30);
    frp_calibration_40 = psi_f(indx40);
    frp_calibration_50 = psi_f(indx50);
    
    h= plot(axs(1), gamma*frp_calibration_30, gamma, 'Color', colors{idata}, 'marker', markers{idata}, ...
        'linewidth', 1, 'MarkerFace', colors{idata}, 'MarkerEdge', colors{idata}, 'MarkerSize', 4);
    l = legend(h, lgd(idata));
    set(l, 'visible', 'off')
    
    h= plot(axs(2), gamma*frp_calibration_40, gamma, 'Color', colors{idata}, 'marker', markers{idata}, ...
        'linewidth', 1, 'MarkerFace', colors{idata}, 'MarkerEdge', colors{idata}, 'MarkerSize', 4);
    l = legend(h, lgd(idata));
    set(l, 'visible', 'off')
    
    h= plot(axs(3), gamma*frp_calibration_50, gamma, 'Color', colors{idata}, 'marker', markers{idata}, ...
        'linewidth', 1, 'MarkerFace', colors{idata}, 'MarkerEdge', colors{idata}, 'MarkerSize', 4);
    l = legend(h, lgd(idata));
    set(l, 'visible', 'off')
    
    xtxt = xlabel('FRP reduction factor, \psi_f','FontSize',8,...
        'FontName','Times New Roman', 'Interpreter','tex');
    ytxt = ylabel('Partial safety factor, \gamma','FontSize',8,...
        'FontName','Times New Roman', 'Interpreter','tex');
end

%% change the following properties from figure inspector
nfig = length(figs);
for ifig = 1:nfig
    fig = figs(ifig);
    ax = axs(ifig);
    set(fig, 'Units','inches', 'Position', [0, 0, 3.5, 2.625]);
    set(ax, 'Units', 'normalized', 'Position', [0.16, 0.15, 0.75, 0.79]);
    % savefig(figure1, strcat('../figures/HKandACInew/re_', code, '_', fc, '.fig'));
    % set(xtxt1, 'Units', 'normalized', 'Position', [-0.1, 0.5, 0]);
    % set(ytxt1, 'Units', 'normalized', 'Position', [0.5, -0.1, 0]);
end