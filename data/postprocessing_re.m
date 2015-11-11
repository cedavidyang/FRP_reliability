clear;
MAX_YR = 50;
code = input('design code\n', 's');
frp_form = input('frp for shear or flexure: 1 for shear, 2 for flexure\n');
switch lower(code)
    case {'aci', 'acinew'}
        psi_f = 1.00;
        FACTOR_FRP = (0.50:0.05:1.00)';
        beta_T_50 = 3.5;
    case {'hk', 'tr'}
        psi_f = (0.10:0.05:2.00)';
        FACTOR_FRP = 1.00;
        beta_T_50 = 3.8;
    case {'gb'}
        psi_f = (0.10:0.05:2.00)';
        FACTOR_FRP = 1.00;
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
    otherwise
end
beta_T_30 = life2re(30, beta_T_50);
beta_T_40 = life2re(40, beta_T_50);
switch frp_form
    case 1
        beam = 'shear';
        struct_s = load(strcat(code, '_', 'shear', '+side_', fc, num2str(fccov), '.mat'));
        struct_u = load(strcat(code, '_', 'shear', '+U_', fc, num2str(fccov), '.mat'));
        struct_w = load(strcat(code, '_', 'shear', '+W_', fc, num2str(fccov), '.mat'));
        
        npsi = length(psi_f);      
        for ipsi = 1:length(psi_f)            
            % side bonding
            s_RE_col = struct_s.re_data{ipsi}(:, 1, :);
            s_RE_col = s_RE_col(:);
            s_RE_col( isnan(s_RE_col) ) = [];
            norm_RE_s_yr30(ipsi) = mean((s_RE_col-beta_T_30).^2);
            norm_RE_s_yr40(ipsi) = mean((s_RE_col-beta_T_40).^2);
            norm_RE_s_yr50(ipsi) = mean((s_RE_col-beta_T_50).^2);
            mean_RE_s(ipsi) = mean(s_RE_col);
            std_RE_s(ipsi) = std(s_RE_col);
            upper_RE_s(ipsi) = max(s_RE_col);
            lower_RE_s(ipsi) = min(s_RE_col);
            
            % U jacketing
            u_RE_col = struct_u.re_data{ipsi}(:, 1, :);
            u_RE_col = u_RE_col(:);
            u_RE_col( isnan(u_RE_col) ) = [];
            norm_RE_u_yr30(ipsi) = mean((u_RE_col-beta_T_30).^2);
            norm_RE_u_yr40(ipsi) = mean((u_RE_col-beta_T_40).^2);
            norm_RE_u_yr50(ipsi) = mean((u_RE_col-beta_T_50).^2);
            mean_RE_u(ipsi) = mean(u_RE_col);
            std_RE_u(ipsi) = std(u_RE_col);
            upper_RE_u(ipsi) = max(u_RE_col);
            lower_RE_u(ipsi) = min(u_RE_col);
            
            % complete wrapping
            w_RE_col = struct_w.re_data{ipsi}(:, 1, :);
            w_RE_col = w_RE_col(:);
            w_RE_col( isnan(w_RE_col) ) = [];
            norm_RE_w_yr30(ipsi) = mean((w_RE_col-beta_T_30).^2);
            norm_RE_w_yr40(ipsi) = mean((w_RE_col-beta_T_40).^2);
            norm_RE_w_yr50(ipsi) = mean((w_RE_col-beta_T_50).^2);
            mean_RE_w(ipsi) = mean(w_RE_col);
            std_RE_w(ipsi) = std(w_RE_col);
            upper_RE_w(ipsi) = max(w_RE_col);
            lower_RE_w(ipsi) = min(w_RE_col);
        end

        %% post processing of shear
        figs = {}; axs = {};
        mean_RE = {mean_RE_s, mean_RE_u, mean_RE_w};
        upper_RE = {upper_RE_s, upper_RE_u, upper_RE_w};
        lower_RE = {lower_RE_s, lower_RE_u, lower_RE_w};
        norm_RE_yr30 = {norm_RE_s_yr30, norm_RE_u_yr30, norm_RE_w_yr30};
        norm_RE_yr40 = {norm_RE_s_yr40, norm_RE_u_yr40, norm_RE_w_yr40};
        norm_RE_yr50 = {norm_RE_s_yr50, norm_RE_u_yr50, norm_RE_w_yr50};
        xlim_cell = {[0.1,2.0], [0.1,2.0], [0.1, 2.0]};
        for ifig = 1:3
            figs{end+1} = figure;
            axs{end+1} = axes('Parent',figs{end},'FontSize',8,'FontName','Times New Roman');
            box(axs{end},'on');
            hold(axs{end},'all');
            plot(psi_f, mean_RE{ifig}, 'k-', 'LineWidth', 1.5, 'MarkerFace', 'k', 'MarkerSize', 4);
            plot(psi_f, upper_RE{ifig}, 'b--', 'LineWidth', 1, 'MarkerSize', 4);
            plot(psi_f, lower_RE{ifig}, 'b--', 'LineWidth', 1, 'MarkerSize', 4);
            xtxt = xlabel('FRP reduction factor','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            ytxt = ylabel('Reliability index','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');

            figs{end+1} = figure;
            axs{end+1} = axes('Parent',figs{end},'FontSize',8,'FontName','Times New Roman');
            box(axs{end},'on');
            hold(axs{end},'all');
            plot(psi_f, norm_RE_yr30{ifig}, 'ko-', 'LineWidth', 1.5, 'MarkerFace', 'k', 'MarkerSize', 4);
            plot(psi_f, norm_RE_yr40{ifig}, 'b^-', 'LineWidth', 1.5, 'MarkerSize', 4);
            plot(psi_f, norm_RE_yr50{ifig}, 'rv-', 'LineWidth', 1.5, 'MarkerSize', 4);
            xtxt = xlabel('FRP reduction factor','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            ytxt = ylabel('Mean distance from target','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            
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
        psi_f = struct_ic.psi_f;
        
        npsi = length(psi_f);
        
        for ipsi = 1:length(psi_f)            
         
            % IC debonding
            ic_RE_col = struct_ic.re_data{ipsi}(:, 1, :);
            ic_RE_col = ic_RE_col(:);
            ic_RE_col( isnan(ic_RE_col) ) = [];
            norm_RE_ic_yr30(ipsi) = mean((ic_RE_col-beta_T_30).^2);
            norm_RE_ic_yr40(ipsi) = mean((ic_RE_col-beta_T_40).^2);
            norm_RE_ic_yr50(ipsi) = mean((ic_RE_col-beta_T_50).^2);
            mean_RE_ic(ipsi) = mean(ic_RE_col);
            std_RE_ic(ipsi) = std(ic_RE_col);
            upper_RE_ic(ipsi) = max(ic_RE_col);
            lower_RE_ic(ipsi) = min(ic_RE_col);
            
            % FRP rupture
            rup_RE_col = struct_rup.re_data{ipsi}(:, 1, :);
            rup_RE_col = rup_RE_col(:);
            rup_RE_col( isnan(rup_RE_col) ) = [];
            norm_RE_rup_yr30(ipsi) = mean((rup_RE_col-beta_T_30).^2);
            norm_RE_rup_yr40(ipsi) = mean((rup_RE_col-beta_T_40).^2);
            norm_RE_rup_yr50(ipsi) = mean((rup_RE_col-beta_T_50).^2);
            mean_RE_rup(ipsi) = mean(rup_RE_col);
            std_RE_rup(ipsi) = std(rup_RE_col);
            upper_RE_rup(ipsi) = max(rup_RE_col);
            lower_RE_rup(ipsi) = min(rup_RE_col);
        end

        %% post processing of shear
        figs = {}; axs = {};
        mean_RE = {mean_RE_ic, mean_RE_rup};
        upper_RE = {upper_RE_ic, upper_RE_rup};
        lower_RE = {lower_RE_ic, lower_RE_rup};
        norm_RE_yr30 = {norm_RE_ic_yr30, norm_RE_rup_yr30};
        norm_RE_yr40 = {norm_RE_ic_yr40, norm_RE_rup_yr40};
        norm_RE_yr50 = {norm_RE_ic_yr50, norm_RE_rup_yr50};
        xlim_cell = {[0.1,2.0], [0.1,2.0]};
        for ifig = 1:2
            figs{end+1} = figure;
            axs{end+1} = axes('Parent',figs{end},'FontSize',8,'FontName','Times New Roman');
            box(axs{end},'on');
            hold(axs{end},'all');
            plot(psi_f, mean_RE{ifig}, 'k-', 'LineWidth', 1.5, 'MarkerFace', 'k', 'MarkerSize', 4);
            plot(psi_f, upper_RE{ifig}, 'b--', 'LineWidth', 1, 'MarkerSize', 4);
            plot(psi_f, lower_RE{ifig}, 'b--', 'LineWidth', 1, 'MarkerSize', 4);
            xtxt = xlabel('FRP reduction factor','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            ytxt = ylabel('Reliability index','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');

            figs{end+1} = figure;
            axs{end+1} = axes('Parent',figs{end},'FontSize',8,'FontName','Times New Roman');
            box(axs{end},'on');
            hold(axs{end},'all');
            plot(psi_f, norm_RE_yr30{ifig}, 'ko-', 'LineWidth', 1.5, 'MarkerFace', 'k', 'MarkerSize', 4);
            plot(psi_f, norm_RE_yr40{ifig}, 'b^-', 'LineWidth', 1.5, 'MarkerSize', 4);
            plot(psi_f, norm_RE_yr50{ifig}, 'rv-', 'LineWidth', 1.5, 'MarkerSize', 4);
            xtxt = xlabel('FRP reduction factor','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            ytxt = ylabel('Mean distance from target','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            
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