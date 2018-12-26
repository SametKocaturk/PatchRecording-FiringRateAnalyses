% This script uses rasterData saved by "main" script. Calculates firing
% rate and z score peristimulus time histogram(PSTH) with given binSize.
% Sorts them based on stimulation response during the stimulation duration.
% Runs kruska wallis analysis on the z score PSTH and marks significant 
% bins. 
% As outputs, it plots firing rate PSTH and z score PSTH.
%
close all;clear all;clc

load rasterData
rasterData = rasterData(2:end);

binSize = 200e-3;
negBinEdges = flip(-binSize:-binSize:-10);
positiveBinEdges = 0:binSize:13;
binEdges = [negBinEdges positiveBinEdges];% stimulation time is 0. Creates longer than needed in order to cover different length of trials. It will later trancated down to a coverage by trials

for dataNo = 1:length(rasterData)
    spikeCountsInBin = nan(length(rasterData(dataNo).EpRasterData), length(binEdges)-1);
    binsCoveredByTrace = zeros(length(rasterData(dataNo).EpRasterData), length(binEdges));
    traceHistNo = 0;
    for traceNo = 1: length(rasterData(dataNo).EpRasterData)
        if ~isempty(rasterData(dataNo).EpRasterData(traceNo).spikeTime)% skips discarded trials
            traceHistNo = traceHistNo +1; % because of the discarded trials, the number of histogram is different than the number of trials
            trSpikeTimes = rasterData(dataNo).EpRasterData(traceNo).spikeTime;
            spikeCountsInBin(traceHistNo,:) = histcounts(trSpikeTimes - rasterData(dataNo).stimTime,binEdges);
            binsCoveredByTrace(traceHistNo,:) = ( (binEdges > (-rasterData(dataNo).stimTime-binSize)) ) ...
                & ( (binEdges < (rasterData(dataNo).traceDur - rasterData(dataNo).stimTime)) );
        end
    end
    rasterData(dataNo).spikeCountsInBin = nansum(spikeCountsInBin);
    rasterData(dataNo).binsCoveredByData = sum(binsCoveredByTrace(:,1:end-1));
    coveredIdx = sum(binsCoveredByTrace) ~= 0;
    baseLineIdx = (binEdges < 0) & coveredIdx;%bins up to the stimulation time
    stimIdx = (binEdges >= 0) & (binEdges <  rasterData(dataNo).stimDur/1000);%bins during stimulation
    avgBaseline = mean( rasterData(dataNo).spikeCountsInBin(baseLineIdx(1:end-1)) );
    stdBaseline = std( rasterData(dataNo).spikeCountsInBin(baseLineIdx(1:end-1)) );
    if stdBaseline == 0; stdBaseline =1;end
    stimActivity = nanmean( rasterData(dataNo).spikeCountsInBin(stimIdx(1:end-1)) );
    rasterData(dataNo).zScore = (rasterData(dataNo).spikeCountsInBin - avgBaseline) ./ stdBaseline;
    rasterData(dataNo).zScore(~coveredIdx(1:end-1)) = 0;
    
    if stimActivity < avgBaseline - stdBaseline% compares the firing rate between before the stimulation and during the stimulation
        rasterData(dataNo).responseType = 'Decreasing';
    elseif stimActivity > avgBaseline + stdBaseline
        rasterData(dataNo).responseType = 'Increasing';
    else
        rasterData(dataNo).responseType = 'No Change';
    end
    
    
end


%% Plots

% Sort the episodes in _rasterData_ by recording type and response type
attachIdx = arrayfun(@(x) strcmp(x.recordingType,'Attach Mode'),rasterData);attachIdx = find(attachIdx == 1);
wholeIdx = arrayfun(@(x) strcmp(x.recordingType,'Whole Cell'),rasterData);wholeIdx = find(wholeIdx == 1);
increasedResponseIdx = arrayfun(@(x) strcmp(x.responseType,'Increasing'),rasterData);increasedResponseIdx = find(increasedResponseIdx == 1);
decreasedResponseIdx = arrayfun(@(x) strcmp(x.responseType,'Decreasing'),rasterData);decreasedResponseIdx = find(decreasedResponseIdx == 1);
nochangeResponseIdx = arrayfun(@(x) strcmp(x.responseType,'No Change'),rasterData);nochangeResponseIdx = find(nochangeResponseIdx == 1);

%Figure settings
nrFigures = 2;
for i = 1:nrFigures
    fH(i) = figure;
    width = 17.6;%cm
    height = width*0.4260;
    fH(i).PaperUnits = 'centimeters';
    fH(i).PaperOrientation ='portrait';
    fH(i).PaperPosition = [0 0 width height];
    fH(i).PaperPositionMode = 'auto';
end

fontName = 'Arial';
fontSize = 10;
%----------------Spike counts histogram-----------------
figure(fH(1))
subplotH = subplot(2,3,1);
periHistByRecordingTypeAndReponseType(rasterData, attachIdx, decreasedResponseIdx, binEdges, binSize, subplotH);
title(['Attach Mode+' 'Decreasing'],'FontName',fontName,'FontWeight','Bold','FontSize',fontSize);
subplotH = subplot(2,3,2);
periHistByRecordingTypeAndReponseType(rasterData, attachIdx, increasedResponseIdx, binEdges, binSize, subplotH);
title(['Attach Mode+' 'Increasing'],'FontName',fontName,'FontWeight','Bold','FontSize',fontSize);
subplotH = subplot(2,3,3);
periHistByRecordingTypeAndReponseType(rasterData, attachIdx, nochangeResponseIdx, binEdges, binSize, subplotH);
title(['Attach Mode+' 'No Change'],'FontName',fontName,'FontWeight','Bold','FontSize',fontSize);
subplotH = subplot(2,3,4);
periHistByRecordingTypeAndReponseType(rasterData, wholeIdx, decreasedResponseIdx, binEdges, binSize, subplotH);
title(['Whole Cell+' 'Decreasing'],'FontName',fontName,'FontWeight','Bold','FontSize',fontSize);
subplotH = subplot(2,3,5);
periHistByRecordingTypeAndReponseType(rasterData, wholeIdx, increasedResponseIdx, binEdges, binSize, subplotH);
title(['Whole Cell+' 'Increasing'],'FontName',fontName,'FontWeight','Bold','FontSize',fontSize);
subplotH = subplot(2,3,6);
periHistByRecordingTypeAndReponseType(rasterData, wholeIdx, nochangeResponseIdx, binEdges, binSize, subplotH);
title (['Whole Cell+' 'No Change'],'FontName',fontName,'FontWeight','Bold','FontSize',fontSize);
axes('Parent', fH(1), ... % common x and y label for all subplots
    'Units', 'normalized', ...
    'Position', [0, 0, 1, 1], ...
    'Visible', 'off', ...
    'XLim', [0, 1], ...
    'YLim', [0, 1], ...
    'NextPlot', 'add');
text(0.05,0.5, 'Firing Rate (spikes/sec)','FontName',fontName,'FontSize',fontSize, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top',...
    'Rotation',90);
text(0.5,0, 'Time (s)','FontName',fontName,'FontSize',fontSize, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom');




%----------------z score histogram------------------

figure(fH(2))
subplotH = subplot(2,3,1);
periHistZscoreByRecordingTypeAndReponseType(rasterData, attachIdx, decreasedResponseIdx, binEdges, binSize, subplotH);
title(['Attach Mode+' 'Decreasing'],'FontName',fontName,'FontWeight','Bold','FontSize',fontSize);
subplotH = subplot(2,3,2);
periHistZscoreByRecordingTypeAndReponseType(rasterData, attachIdx, increasedResponseIdx, binEdges, binSize, subplotH);
title(['Attach Mode+' 'Increasing'],'FontName',fontName,'FontWeight','Bold','FontSize',fontSize);
subplotH = subplot(2,3,3);
periHistZscoreByRecordingTypeAndReponseType(rasterData, attachIdx, nochangeResponseIdx, binEdges, binSize, subplotH);
title(['Attach Mode+' 'No Change'],'FontName',fontName,'FontWeight','Bold','FontSize',fontSize);
subplotH = subplot(2,3,4);
periHistZscoreByRecordingTypeAndReponseType(rasterData, wholeIdx, decreasedResponseIdx, binEdges, binSize, subplotH);
title(['Whole Cell+' 'Decreasing'],'FontName',fontName,'FontWeight','Bold','FontSize',fontSize);
subplotH = subplot(2,3,5);
periHistZscoreByRecordingTypeAndReponseType(rasterData, wholeIdx, increasedResponseIdx, binEdges, binSize, subplotH);
title(['Whole Cell+' 'Increasing'],'FontName',fontName,'FontWeight','Bold','FontSize',fontSize);
subplotH = subplot(2,3,6);
periHistZscoreByRecordingTypeAndReponseType(rasterData, wholeIdx, nochangeResponseIdx, binEdges, binSize, subplotH);
title (['Whole Cell+' 'No Change'],'FontName',fontName,'FontWeight','Bold','FontSize',fontSize);

axes('Parent', fH(2), ... % common x and y label for all subplots
    'Units', 'normalized', 'Position', [0, 0, 1, 1], ...
    'Visible', 'off', 'NextPlot', 'add');
text(0.05,0.5, 'z Score','FontName',fontName,'FontSize',fontSize, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top',...
    'Rotation',90);
text(0.5,0, 'Time (s)','FontName',fontName,'FontSize',fontSize, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom');

%% Firing rate histogram function
function periHistByRecordingTypeAndReponseType(rasterData, recordingTypeIdx, responseTypeId, binEdges, binSize, subplotHandle)

spikeCountsInBin = zeros(1,length(binEdges(1:end-1)));
binsCoveredByData = zeros(1,length(binEdges(1:end-1)));
for dataNo = intersect(recordingTypeIdx, responseTypeId)
    spikeCountsInBin = spikeCountsInBin + rasterData(dataNo).spikeCountsInBin;
    binsCoveredByData = binsCoveredByData + rasterData(dataNo).binsCoveredByData;
end

coveredIdx = binsCoveredByData ~= 0;
b = bar(subplotHandle,binEdges([coveredIdx false])+binSize/2 , spikeCountsInBin(:,coveredIdx)./binsCoveredByData(:,coveredIdx)*1/binSize, 1);
text(0.9,0.9,['n = ' num2str(length(intersect(recordingTypeIdx, responseTypeId)))],'Units','normalized','HorizontalAlignment','right')
set(gca,'FontName','Arial','fontsize',10)
end
%% z socer function and kruskawallis test for detecting significant bins
function periHistZscoreByRecordingTypeAndReponseType(rasterData, recordingTypeIdx, responseTypeId, binEdges, binSize, subplotHandle)

zscorseInBin = zeros(1,length(binEdges(1:end-1)));
binsCoveredByData = zeros(1,length(binEdges(1:end-1)));
for dataNo = intersect(recordingTypeIdx, responseTypeId)
    zscorseInBin = zscorseInBin + rasterData(dataNo).zScore;
    binsCoveredByData = binsCoveredByData + (rasterData(dataNo).binsCoveredByData > 0);
end

coveredIdx = binsCoveredByData > 0;
b = bar(subplotHandle,binEdges([coveredIdx false])+binSize/2 , zscorseInBin(:,coveredIdx)./binsCoveredByData(coveredIdx), 1);
text(0.9,0.9,['n = ' num2str(length(intersect(recordingTypeIdx, responseTypeId)))],'Units','normalized','HorizontalAlignment','right')
set(gca,'FontName','Arial','fontsize',10)
hold on

%kruskawallis statistics and marks significantly different bins
for dataNo = intersect(recordingTypeIdx, responseTypeId)
    baseLineIdx = (binEdges(1:end-1) < 0) & rasterData(dataNo).binsCoveredByData;
    baseLineIdx = find(baseLineIdx == true);
    baseLineZ = rasterData(dataNo).zScore(:,baseLineIdx);
    baseLineZ = baseLineZ(:);
    afterStimIdx = (binEdges(1:end-1) >= 0) & rasterData(dataNo).binsCoveredByData;
    afterStimIdx = find(afterStimIdx == true);
    kwSubTable = nan (length(baseLineZ),size(afterStimIdx,2)+1);%prestimulation in a column + a column for each bin after stimulation
    kwSubTable(:,1) = baseLineZ;
    kwSubTable(1,2:end) = rasterData(dataNo).zScore(:,afterStimIdx);
    if ~exist('kwTableAll','var')% initiates z score matrix for kruskal wallis analysis. Once its initiated this part isn't executed any more.
        kwTableAll = kwSubTable;
    else
        if size(kwTableAll,2) < size(kwSubTable,2)% adds columns to total table if new neuron has longer trial
            nrRequiredCol = size(kwSubTable,2) - size(kwTableAll,2);
            kwTableAll = [kwTableAll nan(size(kwTableAll,1),nrRequiredCol)];
        elseif size(kwTableAll,2) > size(kwSubTable,2)
            nrRequiredCol = size(kwTableAll,2) - size(kwSubTable,2);
            kwSubTable = [kwSubTable nan(size(kwSubTable,1),nrRequiredCol)];
        end
        kwTableAll = [kwTableAll; kwSubTable];
    end
    
    
    
    
end

if exist('kwTableAll','var')% when no data for given condition, data for kruskawallis doesn't exist
    kwTableAll = [kwTableAll; kwTableAll];
    [p,tbl,stats] = kruskalwallis(kwTableAll,[],'off');
    c = multcompare(stats,'Ctype', 'dunn-sidak','Display','off');
    comp1vsOthersIdx = find (c(:,1) ==1);
    signifIdx = find(c(comp1vsOthersIdx,end) <= 0.05);% finds significantly different bins in comparison to the first bin which consists of baseline
    signfBinsXval = binSize * signifIdx - binSize/2;%significance mark x axis values
    yLims = get(gca,'ylim');
    signfMarkYval = yLims(2) - (yLims(2)-yLims(1)) * 0.05;
    if ~isempty(signfBinsXval)% if there is a significant bin
        plot(signfBinsXval,signfMarkYval, '*r','MarkerSize',5)
    end
end
end% periHistZscoreByRecordingTypeAndReponseType








