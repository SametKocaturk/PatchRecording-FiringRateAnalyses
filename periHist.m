function periHist(data, coveredBinIdx, binEdges, figHandle)
% function periHist(data, coveredBinIdx, binEdges, figHandle)
% Displays a peristimulus time histogram on a given figure.
% INPUT
%   data = a struct with each element containing spike counts in each bin
%   coveredBinIdx = a struct with each element marking the bins covered by
%    the traces
%   binEdges = timestampes of bin edges in seconds
%   figHandle = figure handle showing where to plot the PSTH
% OUTPUT
%   A PSTH figure showing firing rate (Hz) in each bin
spikeCountsInBin = zeros(1,length(binEdges(1:end-1)));
binsCoveredByData = zeros(1,length(binEdges(1:end-1)));
for i = 1:length(data)
    spikeCountsInBin = spikeCountsInBin + data(i).spikeCountsInBin;
    binsCoveredByData = binsCoveredByData + coveredBinIdx(i).coveredIdx;
end

coveredIdx = binsCoveredByData ~= 0;
binSize = mean(diff(binEdges));
b = bar(figHandle,binEdges([coveredIdx false])+binSize/2 , spikeCountsInBin(:,coveredIdx)./binsCoveredByData(:,coveredIdx)*1/binSize, 1);
text(0.9,0.9,['n = ' num2str(length(data))],'Units','normalized','HorizontalAlignment','right')
set(gca,'FontName','Arial','fontsize',10)
end