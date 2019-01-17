function zScorePeriHist(data, coveredBinIdx, binEdges, figHandle)
% function zScorePeriHist(data, coveredBinIdx, binEdges, figHandle)
% Displays a peristimulus z Score of spike counts in each bin on a given figure.
% INPUT
%   data = a struct with each element containing z score in each bin
%   coveredBinIdx = a struct with each element marking the bins covered by
%    the traces
%   binEdges = timestampes of bin edges in seconds
%   figHandle = figure handle showing where to plot the PSTH
% OUTPUT
%   A z Score figure showing z value in each bin
zscorseInBin = zeros(1,length(binEdges(1:end-1)));
binsCoveredByData = zeros(1,length(binEdges(1:end-1)));
for i = 1:length(data)
    zscorseInBin = zscorseInBin + data(i).zScore;
    binsCoveredByData = binsCoveredByData + (coveredBinIdx(i).coveredIdx > 0);
end

coveredIdx = binsCoveredByData > 0;
binSize = mean(diff(binEdges));
b = bar(figHandle,binEdges([coveredIdx false])+binSize/2 , zscorseInBin(:,coveredIdx)./binsCoveredByData(coveredIdx), 1);
text(0.9,0.9,['n = ' num2str(length(data))],'Units','normalized','HorizontalAlignment','right')
set(gca,'FontName','Arial','fontsize',10)