% post processing of Monte Carlo simulations
resistSmp = tmpResistSmp;
nMcSim = length(tmpResistSmp);
% figure setting
f_ = clf;
figure(f_);
set(f_,'Units','Pixels','Position',[102 175 688 475.975]);

%% Convergence of Monte Carlo;
resistMeanConv = zeros(nMcSim, 1);
resistStdConv = zeros(nMcSim, 1);
for iMc = 1:nMcSim
    resistMeanConv(iMc) = mean(resistSmp(1:iMc));
    resistStdConv(iMc) = std(resistSmp(1:iMc));
end
subplot(2,2,1)
plot(resistMeanConv, 'linewidth', 2);
xlabel('Number of simulations')
ylabel('Mean resistance (kN)')
subplot(2,2,2)
plot(resistStdConv, 'linewidth', 2);
xlabel('Number of simulations')
ylabel('Std of resistance (kN)')

%% Distributio fit
subplot(2,2,3);
legh_ = []; legt_ = {};   % handles and text for legend
ax_ = gca;
set(ax_,'Box','on');
hold on;

% --- Plot data originally in dataset "tmpResistSmp data"
Data_ = resistSmp;
[F_,X_] = ecdf(Data_,'Function','cdf'); % compute empirical cdf
Bin_.rule = 1;
[C_,E_] = dfswitchyard('dfhistbins',Data_,[],[],Bin_,F_,X_);
[N_,C_] = ecdfhist(F_,X_,'edges',E_); % empirical pdf from cdf
h_ = bar(C_,N_,'hist');
set(h_,'FaceColor','none','EdgeColor',[0.333333 0 0.666667],...
       'LineStyle','-', 'LineWidth',1);
xlabel('Shear resistance (kN)');
ylabel('Density')
legh_(end+1) = h_;
legt_{end+1} = 'MC simulation';   

% --- Create fit "Normal"
% Nudge axis limits beyond data limits
xlim_ = get(ax_,'XLim');
if all(isfinite(xlim_))
   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
   set(ax_,'XLim',xlim_)
end

x_ = linspace(xlim_(1),xlim_(2),100);
pargs_ = cell(1,2);
[pargs_{:}] = normfit(Data_, 0.05);
p_ = [pargs_{:}];
y_ = normpdf(x_,p_(1), p_(2));

h_ = plot(x_,y_,'Color',[1 0 0],...
          'LineStyle','-', 'LineWidth',2,...
          'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_;
legt_{end+1} = 'Normal';

hold off;
leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'}; 
h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
set(h_,'Interpreter','none');

subplot(2,2,4)
probplot('normal');
legh_ = []; legt_ = {};   % handles and text for legend
ax_ = gca;
set(ax_,'Box','on');
hold on;

h_ = probplot(ax_,Data_,[],[],'noref'); % add to probability plot
set(h_,'Color',[0.333333 0 0.666667],'Marker','o', 'MarkerSize',6);
xlabel('Shear resistance (kN)');
ylabel('Density')
legh_(end+1) = h_;
legt_{end+1} = 'MC simulation';   

h_ = probplot(ax_,@normcdf,p_);
set(h_,'Color',[1 0 0],'LineStyle','-', 'LineWidth',2);
legh_(end+1) = h_;
legt_{end+1} = 'Normal';

hold off;
leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'}; 
h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
set(h_,'Interpreter','none');