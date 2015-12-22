clear
rng(2)

code = input('design code\n', 's');
frp_form = input('frp for shear or flexure: 1 for shear, 2 for flexure\n');
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

switch frp_form
    case 1
        struct_s = load(strcat(code, '_', 'shear', '+side_', fc, num2str(fccov), '.mat'));
        struct_u = load(strcat(code, '_', 'shear', '+U_', fc, num2str(fccov), '.mat'));
        struct_w = load(strcat(code, '_', 'shear', '+W_', fc, num2str(fccov), '.mat'));
        
        rsmp_data = {struct_s.resistSmp, struct_u.resistSmp, struct_w.resistSmp};    
        [nsmp, ncase] = size(struct_s.resistSmp);
        %% post processing of shear
        figs = {}; axs = {};
        for ifig = 1:3
            figs{end+1} = figure;
            axs{end+1} = axes('Parent',figs{end},'FontSize',8,'FontName','Times New Roman');
            box(axs{end},'on');
            hold(axs{end},'all');
            
            indx = randi([1 ncase],1,1);
            rsmp = rsmp_data{1}(:, indx);
            cumMu = cumsum(rsmp)./((1:length(rsmp))');
            cumSigma = sqrt(cumsum(rsmp.^2)./((1:length(rsmp))') - (cumsum(rsmp)./((1:length(rsmp))')).^2);
            cumCov = cumSigma./cumMu;
            subplot(2, 1, 1);
            plot(cumMu, 'b-', 'LineWidth', 1, 'MarkerFace', 'k', 'MarkerSize', 4);
            ytxt = ylabel('Mean','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            subplot(2, 1, 2);
            plot(cumCov, 'b-', 'LineWidth', 1, 'MarkerFace', 'k', 'MarkerSize', 4);
            ytxt = ylabel('CoV','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            xtxt = xlabel('Sample number','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            
            if ifig == 1
                fig = figure;
                ax = axes('Parent',fig,'FontSize',8,'FontName','Times New Roman');
            end
            probplot(ax, 'normal',rsmp)
            ytxt = ylabel('Probablity','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            xtxt = xlabel('Data','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            set(get(ax,'Title'),'String','')
            box(ax,'on');
            hold(ax,'all');
        end
        figs{end+1} = fig;
        axs{end+1} = ax;
    case 2
        struct_ic = load(strcat(code, '_', 'flexure', '+IC_', fc, num2str(fccov), '.mat'));
        struct_rup = load(strcat(code, '_', 'flexure', '+rupture_', fc, num2str(fccov), '.mat'));
        
        rsmp_data = {struct_ic.resistSmp, struct_rup.resistSmp};    
        [nsmp, ncase] = size(struct_ic.resistSmp);
        %% post processing of shear
        figs = {}; axs = {};
        for ifig = 1:2
            figs{end+1} = figure;
            axs{end+1} = axes('Parent',figs{end},'FontSize',8,'FontName','Times New Roman');
            box(axs{end},'on');
            hold(axs{end},'all');
            
            indx = randi([1 ncase],1,1);
            rsmp = rsmp_data{1}(:, indx);
            cumMu = cumsum(rsmp)./((1:length(rsmp))');
            cumSigma = sqrt(cumsum(rsmp.^2)./((1:length(rsmp))') - (cumsum(rsmp)./((1:length(rsmp))')).^2);
            cumCov = cumSigma./cumMu;
            subplot(2, 1, 1);
            plot(cumMu, 'b-', 'LineWidth', 1, 'MarkerFace', 'k', 'MarkerSize', 4);
            ytxt = ylabel('Mean','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            subplot(2, 1, 2);
            plot(cumCov, 'b-', 'LineWidth', 1, 'MarkerFace', 'k', 'MarkerSize', 4);
            ytxt = ylabel('CoV','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            xtxt = xlabel('Sample number','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            
            if ifig == 1
                fig = figure;
                ax = axes('Parent',fig,'FontSize',8,'FontName','Times New Roman');
            end
            probplot(ax, 'normal',rsmp)
            ytxt = ylabel('Probablity','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            xtxt = xlabel('Data','FontSize',8,...
                'FontName','Times New Roman', 'Interpreter','tex');
            set(get(ax,'Title'),'String','')
            box(ax,'on');
            hold(ax,'all');
        end
        figs{end+1} = fig;
        axs{end+1} = ax;
    otherwise
end

%% change the following properties from figure inspector
nfig = length(figs);
for ifig = 1:nfig
    fig = figs{ifig};
    ax = axs{ifig};
    set(fig, 'Units','inches', 'Position', [0, 0, 3.5, 2.625]);
%     set(ax, 'Units', 'normalized', 'Position', [0.16, 0.15, 0.75, 0.79]);
    % savefig(figure1, strcat('../figures/HKandACInew/re_', code, '_', fc, '.fig'));
    % set(xtxt1, 'Units', 'normalized', 'Position', [-0.1, 0.5, 0]);
    % set(ytxt1, 'Units', 'normalized', 'Position', [0.5, -0.1, 0]);
end