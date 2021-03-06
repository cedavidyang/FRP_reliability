% post processing of reliability analysis
IndexHat = TARGET_INDEX;
nCase = N_DESIGN_CASE*nLoadRatio;
norm_RE = zeros(nFactorFrp,1);
mean_RE = zeros(nFactorFrp,1);
std_RE = zeros(nFactorFrp,1);
upper_RE = zeros(nFactorFrp,1);
lower_RE = zeros(nFactorFrp,1);
norm_RE_ductile = zeros(nFactorFrp,1);
norm_RE_brittle = zeros(nFactorFrp, 1);
mean_RE_ductile = zeros(nFactorFrp,1);
std_RE_ductile = zeros(nFactorFrp,1);
upper_RE_ductile = zeros(nFactorFrp,1);
lower_RE_ductile = zeros(nFactorFrp,1);
mean_RE_brittle = zeros(nFactorFrp,1);
std_RE_brittle = zeros(nFactorFrp,1);
upper_RE_brittle = zeros(nFactorFrp,1);
lower_RE_brittle = zeros(nFactorFrp,1);

for i_factor = 1:nFactorFrp
    RE_col = reliabilityResults(:, i_factor, :);
    RE_col = RE_col(:);
    RE_col( isnan(RE_col) ) = [];
    
    norm_RE(i_factor) = mean((RE_col-IndexHat).^2);
    mean_RE(i_factor) = mean(RE_col);
    std_RE(i_factor) = std(RE_col);
%     upper_RE(i_factor) = prctile(RE_col, 95);
%     lower_RE(i_factor) = prctile(RE_col, 5);
    upper_RE(i_factor) = max(RE_col);
    lower_RE(i_factor) = min(RE_col);  
    switch SUB_TEST_DATABASE_NAME
        case {'flexure+IC', 'flexure+ic', 'flexure+rup', 'flexure+rupture'}
            tmpIsSteelYield = [isSteelYielding(:, i_factor);isSteelYielding(:, i_factor)];
            isDuctile = tmpIsSteelYield==1;
            if sum(isDuctile) ~=0
                mean_RE_ductile(i_factor) = mean(RE_col(isDuctile));
                std_RE_ductile(i_factor) = std(RE_col(isDuctile));
                upper_RE_ductile(i_factor) = max(RE_col(isDuctile));
                lower_RE_ductile(i_factor) = min(RE_col(isDuctile));
                norm_RE_ductile(i_factor) = mean((RE_col(isDuctile)-IndexHat).^2);
            end
            if sum(~isDuctile) ~= 0
                mean_RE_brittle(i_factor) = mean(RE_col(~isDuctile));
                std_RE_brittle(i_factor) = std(RE_col(~isDuctile));
                upper_RE_brittle(i_factor) = max(RE_col(~isDuctile));
                lower_RE_brittle(i_factor) = min(RE_col(~isDuctile));
                norm_RE_brittle(i_factor) = mean((RE_col(isDuctile)-(IndexHat+0.5)).^2);
            end
        otherwise
    end
    nBeta = length(RE_col);
    if i_factor == 1
        betaData = RE_col;
        iFactor = FACTOR_FRP(i_factor)*ones(nBeta, 1);
    else
        betaData = [betaData; RE_col];
        iFactor = [iFactor; FACTOR_FRP(i_factor)*ones(nBeta, 1)];
    end
end

figure1 = figure;
axes1 = axes('Parent',figure1,'FontSize',12,'FontName','Times New Roman');
box(axes1,'on');
hold(axes1,'all');

h1 = plot(FACTOR_FRP, mean_RE, 'MarkerFaceColor','b','Marker','o',...
          'LineWidth', 2, 'Color', 'b', 'LineStyle', '-');
h3 = plot(FACTOR_FRP, upper_RE, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
plot(FACTOR_FRP, lower_RE, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
h2 = plot(FACTOR_FRP, mean_RE+std_RE, 'LineWidth', 1, 'Color', [0 0.498039215803146 0], 'LineStyle', '-.');
plot(FACTOR_FRP, mean_RE-std_RE, 'LineWidth', 1, 'Color', [0 0.498039215803146 0], 'LineStyle', '-.');
% plot(iFactor, betaData, '.')
if strcmpi(DESIGN_CODE, 'aci')
    xlabel('Strength reduction factor for FRP contribution','FontSize',12,...
        'FontName','Times New Roman');
else
    xlabel('Partial safety factor for FRP contribution','FontSize',12,...
        'FontName','Times New Roman');
end
ylabel('Reliability index,  \beta','FontSize',12,...
    'FontName','Times New Roman');
legend([h1, h2, h3], 'Average reliability', '\beta \pm std', '\beta min & max');
% Create textarrow
annotation(figure1,'textarrow',[0.730357142857143 0.821428571428571],...
    [0.601380952380953 0.323809523809524],'TextEdgeColor','none','FontSize',12,...
    'FontName','Times New Roman',...
    'String',{'\psi_f = 0.95','\beta = 3.613'});
savefig('./figures/fig1.fig');
print('./figures/fig1', '-depsc');

figure;
plot(FACTOR_FRP, norm_RE, 'bo-', 'linewidth', 2, 'markerfacecolor', 'b')
if strcmpi(DESIGN_CODE, 'aci')
    xlabel('Strength reduction factor for FRP contribution','FontSize',12,...
        'FontName','Times New Roman');
else
    xlabel('Partial safety factor for FRP contribution','FontSize',12,...
        'FontName','Times New Roman');
end
ylabel('Norm of Reliability Index','FontSize',12,...
    'FontName','Times New Roman');
savefig('./figures/fig2.fig');

switch SUB_TEST_DATABASE_NAME
    case {'flexure+IC', 'flexure+ic', 'flexure+rup', 'flexure+rupture'}
        figure1 = figure;
        axes1 = axes('Parent',figure1,'FontSize',12,'FontName','Times New Roman');
        box(axes1,'on');
        hold(axes1,'all');
        h1 = plot(FACTOR_FRP, mean_RE_ductile, 'MarkerFaceColor','b','Marker','o',...
            'LineWidth', 2, 'Color', 'b', 'LineStyle', '-');
        h3 = plot(FACTOR_FRP, upper_RE_ductile, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
        plot(FACTOR_FRP, lower_RE_ductile, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
        h2 = plot(FACTOR_FRP, mean_RE_ductile+std_RE_ductile, 'LineWidth', 1, 'Color',[0 0.498039215803146 0], 'LineStyle', '-.');
        plot(FACTOR_FRP, mean_RE_ductile-std_RE_ductile, 'LineWidth', 1, 'Color', [0 0.498039215803146 0], 'LineStyle', '-.');
        xlabel('Partial Safety Factor for FRP contribution');
        if strcmpi(DESIGN_CODE, 'aci')
            xlabel('Strength reduction factor for FRP contribution','FontSize',12,...
                'FontName','Times New Roman');
        else
            xlabel('Partial safety factor for FRP contribution','FontSize',12,...
                'FontName','Times New Roman');
        end
        ylabel('Reliability index,  \beta','FontSize',12,...
            'FontName','Times New Roman');
        legend([h1, h2, h3], 'Average reliability', '\beta \pm std', '\beta min & max');
        % Create textarrow
        annotation(figure1,'textarrow',[0.730357142857143 0.821428571428571],...
            [0.601380952380953 0.323809523809524],'TextEdgeColor','none','FontSize',12,...
            'FontName','Times New Roman',...
            'String',{'\psi_f = 0.95','\beta = 3.613'});
        savefig('./figures/fig11.fig');
        print('./figures/fig11', '-depsc');
        
        figure;
        plot(FACTOR_FRP, norm_RE_ductile, 'bo-', 'linewidth', 2, 'markerfacecolor', 'b')
        if strcmpi(DESIGN_CODE, 'aci')
            xlabel('Strength reduction factor for FRP contribution','FontSize',12,...
                'FontName','Times New Roman');
        else
            xlabel('Partial safety factor for FRP contribution','FontSize',12,...
                'FontName','Times New Roman');
        end
        ylabel('Norm of Reliability Index','FontSize',12,...
            'FontName','Times New Roman');
        hold on
        plot(FACTOR_FRP, norm_RE_brittle, 'go-', 'linewidth', 2, 'markerfacecolor', 'g')
        savefig('./figures/fig21.fig');

        figure1 = figure;
        axes1 = axes('Parent',figure1,'FontSize',12,'FontName','Times New Roman');
        box(axes1,'on');
        hold(axes1,'all');
        h1 = plot(FACTOR_FRP, mean_RE_brittle, 'MarkerFaceColor','b','Marker','o',...
            'LineWidth', 2, 'Color', 'b', 'LineStyle', '-');
        h3 = plot(FACTOR_FRP, upper_RE_brittle, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
        plot(FACTOR_FRP, lower_RE_brittle, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
        h2 = plot(FACTOR_FRP, mean_RE_brittle+std_RE_brittle, 'LineWidth', 1, 'Color', [0 0.498039215803146 0], 'LineStyle', '-.');
        plot(FACTOR_FRP, mean_RE_brittle-std_RE_brittle, 'LineWidth', 1, 'Color', [0 0.498039215803146 0], 'LineStyle', '-.');
        if strcmpi(DESIGN_CODE, 'aci')
            xlabel('Strength reduction factor for FRP contribution','FontSize',12,...
                'FontName','Times New Roman');
        else
            xlabel('Partial safety factor for FRP contribution','FontSize',12,...
                'FontName','Times New Roman');
        end
        ylabel('Reliability index,  \beta','FontSize',12,...
            'FontName','Times New Roman');
        legend([h1, h2, h3], 'Average reliability', '\beta \pm std', '\beta min & max');        
        figure;
        plot(FACTOR_FRP, norm_RE_brittle, 'bo-', 'linewidth', 2, 'markerfacecolor', 'b')
        if strcmpi(DESIGN_CODE, 'aci')
            xlabel('Strength reduction factor for FRP contribution','FontSize',12,...
                'FontName','Times New Roman');
        else
            xlabel('Partial safety factor for FRP contribution','FontSize',12,...
                'FontName','Times New Roman');
        end
        ylabel('Norm of Reliability Index','FontSize',12,...
            'FontName','Times New Roman');
        hold on
        plot(FACTOR_FRP, norm_RE_brittle, 'go-', 'linewidth', 2, 'markerfacecolor', 'g')
    otherwise
end

% noFactorFrp = find(FACTOR_FRP==1);
% betaExp = reliabilityResults(:, noFactorFrp, 1);
% figure; boxplot(betaExp, H_DESIGN_ARRAY_MM);

% switch SUB_TEST_DATABASE_NAME
%     case {'shear+side', 'shear+U', 'shear+W'}
%         attention = find( isOverReinforce(:,noFactorFrp)==0 & betaExp<2.8 );
%         attention = attention(:);
%     case {'flexure+IC', 'flexure+ic', 'flexure+rup', 'flexure+rupture'}
%     otherwise
% end

%% nine typical PDF of resistance
% figure;
% nRepresentCase = 9;
% representCase = randi(N_DESIGN_CASE, nRepresentCase, 1);
% resistRepresentSmp = resistSmp(:, representCase);
% for iRepresentCase = 1:nRepresentCase
%     subplot( ceil(sqrt(nRepresentCase)), ceil(sqrt(nRepresentCase)), iRepresentCase);
%     tmpSmp = resistRepresentSmp(:,iRepresentCase);
%     probplot(tmpSmp); box on;
%     xlabel(['iCase = ', num2str( representCase(iRepresentCase) )]);
% end

%% cov of resistance
figure;
plot(resistStd./resistMean*100, 'o');
savefig('./figures/fig3.fig');

% convergMean = zeros(N_MC, 1);
% convergStd = zeros(N_MC, 1);
% for iMC = 1:N_MC
%     convergMean(iMC) = mean(tmpSmp(1:iMC));
%     convergStd(iMC) = std(tmpSmp(1:iMC));
% end
% figure;
% plot(convergMean);
% figure;
% plot(convergStd);