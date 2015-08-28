function [meanConv, stdConv] = monte_carlo_conv(mcData)
% postprocessing monte carlo data
nMc = length(mcData);
meanConv = zeros(nMc, 1);
stdConv = zeros(nMc, 1);
for iMc = 1:nMc
    meanConv(iMc) = mean(mcData(1:iMc));
    stdConv(iMc) = std(mcData(1:iMc));
end

return
end