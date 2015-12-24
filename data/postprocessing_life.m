clear;

frp_form = input('frp for shear or flexure: 1 for shear, 2 for flexure\n');
code = input('design code\n', 's');
psi_f = (0.10:0.05:2.00)';
life = 10:5:50;
beta_T = [3.0, 3.2, 3.5, 3.7];

fc_type = input('fc type: 1 for fcdet, 2 for 125 bias, 3 for cov10\n');
switch fc_type
    case 1
        fc = 'fcdet'; fccov=0.2;
    case 2
        fc = 'fc125'; fccov=0.2;
    case 3
        fc = 'cov10'; fccov=0.1;
    case 5
        fc = 'cov15'; fccov=0.15;    
    otherwise
end

switch frp_form
    case 1
        struct_s = load(strcat(code, '_', 'shear', '+side_', fc, num2str(fccov), '.mat'));
        struct_u = load(strcat(code, '_', 'shear', '+U_', fc, num2str(fccov), '.mat'));
        struct_w = load(strcat(code, '_', 'shear', '+W_', fc, num2str(fccov), '.mat'));

        data_struct_array = {struct_s, struct_u, struct_w}';
        psi_array = zeros(4, 3, length(life));
        norm_RE = [];


        for ibeta = 1:4
            for idata = 1:3
                for iyear = 1:length(life)
                    ilife = life(iyear);
                    beta_T_life = life2re(ilife, beta_T(ibeta));
                    [~, ~, ~, ~, ~, ~, ~, norm_RE_array] = get_re_data( ...
                        data_struct_array{idata}, psi_f, code, beta_T_life);
                    cubespline = csapi(psi_f, norm_RE_array');
                    [norm_min, psi] = fnmin(cubespline, [0.1, 2.00]);
                    psi_array(ibeta, idata, iyear) = psi;
                end
            end
        end

    case 2
        struct_ic = load(strcat(code, '_', 'flexure', '+IC_', fc, num2str(fccov), '.mat'));
        struct_rup = load(strcat(code, '_', 'flexure', '+rupture_', fc, num2str(fccov), '.mat'));
        
        data_struct_array = {struct_ic, struct_rup}';
        psi_array = zeros(4, 3, length(life));
        norm_RE = [];
        
        for ibeta = 1:4
            for idata = 1:2
                for iyear = 1:length(life)
                    ilife = life(iyear);
                    beta_T_life = life2re(ilife, beta_T(ibeta));
                    [~, ~, ~, ~, ~, ~, ~, norm_RE_array] = get_re_data( ...
                        data_struct_array{idata}, psi_f, code, beta_T_life);
                    if idata==1
                        cubespline = csapi(psi_f, norm_RE_array');
                        [norm_min, psi] = fnmin(cubespline, [0.1, 2.00]);
                    else
                        cubespline = csapi(psi_f(psi_f<=1.00), norm_RE_array(psi_f<=1.00)');
                        [norm_min, psi] = fnmin(cubespline, [0.1, 1.00]);
                    end
                    psi_array(ibeta, idata, iyear) = psi;
                end
            end
        end
end     

%% post processing of shear
figs = {}; axs = {};
for ifig = 1:4
    figs{end+1} = figure;
    axs{end+1} = axes('Parent',figs{end},'FontSize',8,'FontName','Times New Roman');
    box(axs{end},'on');
    hold(axs{end},'all');
    
    switch frp_form
        case 1
            h1 = plot(axs{end}, life, squeeze(psi_array(ifig, 1, :)),...
                'o-', 'Color', 'b', 'LineWidth', 1, 'MarkerFace', 'b', 'MarkerEdge', 'b', 'MarkerSize', 4);
            h2 = plot(axs{end}, life, squeeze(psi_array(ifig, 2, :)),...
                '^-', 'Color', [0, 0.5, 0], 'LineWidth', 1, 'MarkerFace', [0, 0.5, 0], 'MarkerEdge', [0, 0.5, 0], 'MarkerSize', 4);
            h3 = plot(axs{end}, life, squeeze(psi_array(ifig, 3, :)),...
                'v-', 'Color', 'r', 'LineWidth', 1, 'MarkerFace', 'r', 'MarkerEdge', 'r', 'MarkerSize', 4);
            lgd = legend([h1, h2, h3], {'Side bonding', 'U-jacketing', 'Complete wrapping'});
        case 2
            h1 = plot(axs{end}, life, squeeze(psi_array(ifig, 1, :)),...
                'o-', 'Color', 'b', 'LineWidth', 1, 'MarkerFace', 'b', 'MarkerEdge', 'b', 'MarkerSize', 4);
            h2 = plot(axs{end}, life, squeeze(psi_array(ifig, 2, :)),...
                '^-', 'Color', [0, 0.5, 0], 'LineWidth', 1, 'MarkerFace', [0, 0.5, 0], 'MarkerEdge', [0, 0.5, 0], 'MarkerSize', 4);
            lgd = legend([h1, h2], {'Debonding', 'Rupture'});            
    end
    set(lgd, 'Color', 'None')
    legend boxoff
    xtxt = xlabel('Service life (year)','FontSize',8,...
        'FontName','Times New Roman', 'Interpreter','tex');
    switch lower(code)
        case {'aci', 'acinew'}
            ytxt = ylabel('FRP reduction factor, \psi_f','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
        case {'gb', 'fib', 'tr', 'gbnew'}
            ytxt = ylabel('FRP calibration factor, \phi_f','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
    end
end

%% change the following properties from figure inspector
nfig = length(figs);
for ifig = 1:nfig
    fig = figs{ifig};
    ax = axs{ifig};
    set(fig, 'Units','inches', 'Position', [0, 0, 3.5, 2.625]);
    set(ax, 'Units', 'normalized', 'Position', [0.16, 0.15, 0.75, 0.79]);
    % savefig(figure1, strcat('../figures/HKandACInew/re_', code, '_', fc, '.fig'));
    % set(xtxt1, 'Units', 'normalized', 'Position', [-0.1, 0.5, 0]);
    % set(ytxt1, 'Units', 'normalized', 'Position', [0.5, -0.1, 0]);
end