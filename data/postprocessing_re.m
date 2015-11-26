clear;

code = input('design code\n', 's');
frp_form = input('frp for shear or flexure: 1 for shear, 2 for flexure\n');
switch lower(code)
    case {'aci', 'acinew'}
        psi_f = (0.10:0.05:2.00)';
        beta_T_50 = 3.5;
    case {'hk', 'tr'}
        psi_f = (0.10:0.05:2.00)';
        beta_T_50 = 3.8;
    case {'gb'}
        psi_f = (0.10:0.05:2.00)';
        beta_T_50 = 3.2;    
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
        smp_RE = {};
        mean_RE = {};
        std_RE = {};
        upper_RE = {};
        lower_RE = {};
        norm_RE_yr30 = {};
        norm_RE_yr40 = {};
        norm_RE_yr50 = {};

        for idata = 1:3
            [re_smpi, mean_REi, std_REi, upper_REi, lower_REi, norm_RE_yr30i, ...
                norm_RE_yr40i, norm_RE_yr50i] = get_re_data( ...
                data_struct_array{idata}, psi_f, code, beta_T_50);
            smp_RE{end+1} = re_smpi;
            mean_RE{end+1} = mean_REi;
            std_RE{end+1} = std_REi;
            upper_RE{end+1} = upper_REi;
            lower_RE{end+1} = lower_REi;
            norm_RE_yr30{end+1} = norm_RE_yr30i;
            norm_RE_yr40{end+1} = norm_RE_yr40i;
            norm_RE_yr50{end+1} = norm_RE_yr50i;
        end

        %% post processing of shear
        figs = {}; axs = {};
        xlim_cell = {[0.1,2.0], [0.1,2.0], [0.1, 2.0]};
        for ifig = 1:3
            figs{end+1} = figure;
            axs{end+1} = axes('Parent',figs{end},'FontSize',8,'FontName','Times New Roman');
            box(axs{end},'on');
            hold(axs{end},'all');
            h1 = plot(psi_f, mean_RE{ifig}, 'k-', 'LineWidth', 1.5, 'MarkerFace', 'k', 'MarkerSize', 4);
            h2 = plot(psi_f, upper_RE{ifig}, 'b-', 'LineWidth', 1, 'MarkerSize', 4);
            h3 = plot(psi_f, lower_RE{ifig}, 'b-', 'LineWidth', 1, 'MarkerSize', 4);
            X=[psi_f', fliplr(psi_f')];                %#create continuous x value array for plotting
            Y=[upper_RE{ifig},fliplr(lower_RE{ifig})];              %#create continuous x value array for plotting
            h4 = fill(X, Y, [0.99, 0.92, 0.8]);
            uistack(h4,'bottom');
            set(axs{end}, 'Layer', 'top')
            legend([h1, h2, h3], {'mean', 'upper bound', 'load bound'})
            xtxt = xlabel('FRP calibration factor, \phi_f','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            ytxt = ylabel('Reliability index','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');

            figs{end+1} = figure;
            axs{end+1} = axes('Parent',figs{end},'FontSize',8,'FontName','Times New Roman');
            box(axs{end},'on');
            hold(axs{end},'all');
            h1 = plot(psi_f, norm_RE_yr30{ifig}, 'ko-', 'LineWidth', 1, 'MarkerFace', 'k', 'MarkerEdge', 'k', 'MarkerSize', 4);
            h2 = plot(psi_f, norm_RE_yr40{ifig}, 'b^-', 'LineWidth', 1, 'MarkerFace', 'b', 'MarkerEdge', 'b', 'MarkerSize', 4);
            h3 = plot(psi_f, norm_RE_yr50{ifig}, 'rv-', 'LineWidth', 1, 'MarkerFace', 'r', 'MarkerEdge', 'r', 'MarkerSize', 4);
            xtxt = xlabel('FRP calibration factor, \phi_f','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            ytxt = ylabel('Mean distance from target','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            legend([h1, h2, h3], {'30 years', '40 years', '50 years'})
            
            switch code
                case 'gb'
                    xlim(axs{end},xlim_cell{ifig});
                    xlim(axs{end-1},xlim_cell{ifig});
                otherwise
            end  
        end
    case 2
        struct_ic = load(strcat(code, '_', 'flexure', '+IC_', fc, num2str(fccov), '.mat'));
        struct_rup = load(strcat(code, '_', 'flexure', '+rupture_', fc, num2str(fccov), '.mat'));
        
        data_struct_array = {struct_ic, struct_rup}';
        smp_RE = {};
        mean_RE = {};
        std_RE = {};
        upper_RE = {};
        lower_RE = {};
        norm_RE_yr30 = {};
        norm_RE_yr40 = {};
        norm_RE_yr50 = {};
        
        for idata = 1:2
            [re_smpi, mean_REi, std_REi, upper_REi, lower_REi, norm_RE_yr30i, ...
                norm_RE_yr40i, norm_RE_yr50i] = get_re_data( ...
                data_struct_array{idata}, psi_f, code, beta_T_50);
            smp_RE{end+1} = re_smpi;
            mean_RE{end+1} = mean_REi;
            std_RE{end+1} = std_REi;
            upper_RE{end+1} = upper_REi;
            lower_RE{end+1} = lower_REi;
            norm_RE_yr30{end+1} = norm_RE_yr30i;
            norm_RE_yr40{end+1} = norm_RE_yr40i;
            norm_RE_yr50{end+1} = norm_RE_yr50i;
        end

        %% post processing of shear
        figs = {}; axs = {};
        xlim_cell = {[0.1,2.0], [0.1,2.0]};
        for ifig = 1:2
            figs{end+1} = figure;
            axs{end+1} = axes('Parent',figs{end},'FontSize',8,'FontName','Times New Roman');
            box(axs{end},'on');
            hold(axs{end},'all');
            h1 = plot(psi_f, mean_RE{ifig}, 'k-', 'LineWidth', 1.5, 'MarkerFace', 'k', 'MarkerSize', 4);
            h2 = plot(psi_f, upper_RE{ifig}, 'b-', 'LineWidth', 1, 'MarkerSize', 4);
            h3 = plot(psi_f, lower_RE{ifig}, 'b-', 'LineWidth', 1, 'MarkerSize', 4);
            X=[psi_f', fliplr(psi_f')];                %#create continuous x value array for plotting
            Y=[upper_RE{ifig},fliplr(lower_RE{ifig})];              %#create continuous x value array for plotting
            h4 = fill(X, Y, [0.99, 0.92, 0.8]);
            uistack(h4,'bottom');
            set(axs{end}, 'Layer', 'top')
            legend([h1, h2, h3], {'mean', 'upper bound', 'load bound'})
            xtxt = xlabel('FRP calibration factor, \phi_f','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            ytxt = ylabel('Reliability index','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');

            figs{end+1} = figure;
            axs{end+1} = axes('Parent',figs{end},'FontSize',8,'FontName','Times New Roman');
            box(axs{end},'on');
            hold(axs{end},'all');
            h1 = plot(psi_f, norm_RE_yr30{ifig}, 'ko-', 'LineWidth', 1, 'MarkerFace', 'k', 'MarkerEdge', 'k', 'MarkerSize', 4);
            h2 = plot(psi_f, norm_RE_yr40{ifig}, 'b^-', 'LineWidth', 1, 'MarkerFace', 'b', 'MarkerEdge', 'b', 'MarkerSize', 4);
            h3 = plot(psi_f, norm_RE_yr50{ifig}, 'rv-', 'LineWidth', 1, 'MarkerFace', 'r', 'MarkerEdge', 'r', 'MarkerSize', 4);
            xtxt = xlabel('FRP calibration factor, \phi_f','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            ytxt = ylabel('Mean distance from target','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            legend([h1, h2, h3], {'30 years', '40 years', '50 years'})
            
            switch code
                case 'gb'
                    xlim(axs{end},xlim_cell{ifig});
                    xlim(axs{end-1},xlim_cell{ifig});
                otherwise
            end  
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