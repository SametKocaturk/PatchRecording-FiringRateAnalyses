close all;clear all;clc

load rasterData
rasterData = rasterData(2:end);

binSize = 200e-3;
negBinEdges = flip(-binSize:-binSize:-4);
positiveBinEdges = 0:binSize:13;
binEdges = [negBinEdges positiveBinEdges];

for dataNo = 1:length(rasterData)
    spikeCountsInBin = zeros(length(rasterData(dataNo).EpRasterData), length(binEdges)-1);
    binsCoveredByTrace = nan(length(rasterData(dataNo).EpRasterData), length(binEdges));
    traceHistNo = 0;
    for traceNo = 1: length(rasterData(dataNo).EpRasterData)
        if ~isempty(rasterData(dataNo).EpRasterData(traceNo).spikeTime)
            traceHistNo = traceHistNo +1;
            trSpikeTimes = rasterData(dataNo).EpRasterData(traceNo).spikeTime;
            spikeCountsInBin(traceHistNo,:) = histcounts(trSpikeTimes - rasterData(dataNo).stimTime,binEdges);
            binsCoveredByTrace(traceHistNo,:) = ( (binEdges >= -rasterData(dataNo).stimTime) ) & ( (binEdges <= rasterData(dataNo).traceDur - rasterData(dataNo).stimTime) );
        end
    end
    rasterData(dataNo).spikeCountsInBin = sum(spikeCountsInBin);
    rasterData(dataNo).binsCoveredByData = nansum(binsCoveredByTrace(:,1:end-1));
    coveredIdx = nansum(binsCoveredByTrace) ~= 0;
    baseLineIdx = (binEdges < 0) & coveredIdx;%bins up to the stimulation time
    stimIdx = (binEdges >= 0) & (binEdges <  rasterData(dataNo).stimDur/1000);
    avgBaseline = mean( rasterData(dataNo).spikeCountsInBin(baseLineIdx(1:end-1)) );
    stdBaseline = std( rasterData(dataNo).spikeCountsInBin(baseLineIdx(1:end-1)) );
    stimActivity = nanmean( rasterData(dataNo).spikeCountsInBin(stimIdx(1:end-1)) );
    rasterData(dataNo).zScore = (rasterData(dataNo).spikeCountsInBin - avgBaseline) ./ stdBaseline;
    
    if stimActivity < avgBaseline - stdBaseline
        rasterData(dataNo).responseType = 'Decreasing';
    elseif stimActivity > avgBaseline + stdBaseline
        rasterData(dataNo).responseType = 'Increasing';
    else
        rasterData(dataNo).responseType = 'No Change';
    end
    
    
end


%%

attachIdx = arrayfun(@(x) strcmp(x.recordingType,'Attach Mode'),rasterData);attachIdx = find(attachIdx == 1);
wholeIdx = arrayfun(@(x) strcmp(x.recordingType,'Whole Cell'),rasterData);wholeIdx = find(wholeIdx == 1);
increasedResponseIdx = arrayfun(@(x) strcmp(x.responseType,'Increasing'),rasterData);increasedResponseIdx = find(increasedResponseIdx == 1);
decreasedResponseIdx = arrayfun(@(x) strcmp(x.responseType,'Decreasing'),rasterData);decreasedResponseIdx = find(decreasedResponseIdx == 1);
nochangeResponseIdx = arrayfun(@(x) strcmp(x.responseType,'No Change'),rasterData);nochangeResponseIdx = find(nochangeResponseIdx == 1);

%Spike counts histogram
subplotH = subplot(2,3,1);
periHistByRecordingTypeAndReponseType(rasterData, attachIdx, decreasedResponseIdx, binEdges, binSize, subplotH);
title(['Attach Mode+' 'Decreasing']);
subplotH = subplot(2,3,2);
periHistByRecordingTypeAndReponseType(rasterData, attachIdx, increasedResponseIdx, binEdges, binSize, subplotH);
title(['Attach Mode+' 'Increasing'])
subplotH = subplot(2,3,3);
periHistByRecordingTypeAndReponseType(rasterData, attachIdx, nochangeResponseIdx, binEdges, binSize, subplotH);
title(['Attach Mode+' 'No Change'])

subplotH = subplot(2,3,4);
periHistByRecordingTypeAndReponseType(rasterData, wholeIdx, decreasedResponseIdx, binEdges, binSize, subplotH);
title(['Whole Cell+' 'Decreasing'])
subplotH = subplot(2,3,5);
periHistByRecordingTypeAndReponseType(rasterData, wholeIdx, increasedResponseIdx, binEdges, binSize, subplotH);
title(['Whole Cell+' 'Increasing'])
subplotH = subplot(2,3,6);
periHistByRecordingTypeAndReponseType(rasterData, wholeIdx, nochangeResponseIdx, binEdges, binSize, subplotH);
title (['Whole Cell+' 'No Change'])

%z score histogram

figure
subplotH = subplot(2,3,1);
periHistZscoreByRecordingTypeAndReponseType(rasterData, attachIdx, decreasedResponseIdx, binEdges, binSize, subplotH);
title(['Attach Mode+' 'Decreasing']);
subplotH = subplot(2,3,2);
periHistZscoreByRecordingTypeAndReponseType(rasterData, attachIdx, increasedResponseIdx, binEdges, binSize, subplotH);
title(['Attach Mode+' 'Increasing'])
subplotH = subplot(2,3,3);
periHistZscoreByRecordingTypeAndReponseType(rasterData, attachIdx, nochangeResponseIdx, binEdges, binSize, subplotH);
title(['Attach Mode+' 'No Change'])

subplotH = subplot(2,3,4);
periHistZscoreByRecordingTypeAndReponseType(rasterData, wholeIdx, decreasedResponseIdx, binEdges, binSize, subplotH);
title(['Whole Cell+' 'Decreasing'])
subplotH = subplot(2,3,5);
periHistZscoreByRecordingTypeAndReponseType(rasterData, wholeIdx, increasedResponseIdx, binEdges, binSize, subplotH);
title(['Whole Cell+' 'Increasing'])
subplotH = subplot(2,3,6);
periHistZscoreByRecordingTypeAndReponseType(rasterData, wholeIdx, nochangeResponseIdx, binEdges, binSize, subplotH);
title (['Whole Cell+' 'No Change'])

%% 
function periHistByRecordingTypeAndReponseType(rasterData, recordingTypeIdx, responseTypeId, binEdges, binSize, subplotHandle)

spikeCountsInBin = zeros(1,length(binEdges(1:end-1)));
binsCoveredByData = zeros(1,length(binEdges(1:end-1)));
for i = intersect(recordingTypeIdx, responseTypeId)
    spikeCountsInBin = spikeCountsInBin + rasterData(i).spikeCountsInBin;
    binsCoveredByData = binsCoveredByData + rasterData(i).binsCoveredByData;
end

coveredIdx = binsCoveredByData ~= 0;
b = bar(subplotHandle,binEdges([coveredIdx false])+binSize/2 , spikeCountsInBin(:,coveredIdx)./binsCoveredByData(:,coveredIdx)*1/binSize, 1);
text(0.9,0.9,['n = ' num2str(length(intersect(recordingTypeIdx, responseTypeId)))],'Units','normalized','HorizontalAlignment','right')
ylabel('firing rate (spikes/sec)')
end

function periHistZscoreByRecordingTypeAndReponseType(rasterData, recordingTypeIdx, responseTypeId, binEdges, binSize, subplotHandle)

zscorseInBin = zeros(1,length(binEdges(1:end-1)));
binsCoveredByData = zeros(1,length(binEdges(1:end-1)));
for i = intersect(recordingTypeIdx, responseTypeId)
    zscorseInBin = zscorseInBin + rasterData(i).zScore;
    binsCoveredByData = binsCoveredByData + rasterData(i).binsCoveredByData;
end

coveredIdx = binsCoveredByData ~= 0;
b = bar(subplotHandle,binEdges([coveredIdx false])+binSize/2 , zscorseInBin(:,coveredIdx)./length(intersect(recordingTypeIdx, responseTypeId)), 1);
text(0.9,0.9,['n = ' num2str(length(intersect(recordingTypeIdx, responseTypeId)))],'Units','normalized','HorizontalAlignment','right')
ylabel('z score')
end











