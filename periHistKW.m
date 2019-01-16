% This script uses rasterData saved by "main.m" script. Calculates firing
% rate and z score peristimulus time histogram(PSTH) with given binSize.
% Sorts them based on stimulation response during the stimulation duration.
% Runs kruska wallis analysis on the z score PSTH and marks significant
% bins. The baseline activity is compared with the rest of the bins.
% As outputs, it plots firing rate PSTH and z score PSTH in a formated
% figure.
%
close all;clear all;clc

load rasterData
rasterData = rasterData(2:end);

binSize = 200e-3;

for dataNo = 1:length(rasterData)
    [spikeCountsInBin,zScore,avgBaseline,stdBaseline,binsCoveredByTrace,binEdges] = binningZscoring...
        (binSize,rasterData(dataNo).EpRasterData,rasterData(dataNo).stimTime,rasterData(dataNo).traceDur);
    rasterData(dataNo).spikeCountsInBin = spikeCountsInBin;
    rasterData(dataNo).binsCoveredByData = sum(binsCoveredByTrace(:,1:end-1));
    rasterData(dataNo).zScore = zScore;
    stimIdx = (binEdges >= 0) & (binEdges <  rasterData(dataNo).stimDur/1000);%bins during stimulation
    stimActivity = nanmean( rasterData(dataNo).spikeCountsInBin(stimIdx(1:end-1)) );
    
    
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

%Figure properties
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

conditionIdx = { attachIdx, decreasedResponseIdx;...
    attachIdx, increasedResponseIdx;...
    attachIdx, nochangeResponseIdx;...
    wholeIdx, decreasedResponseIdx;...
    wholeIdx, increasedResponseIdx;...
    wholeIdx, nochangeResponseIdx };
titlesSubFig = {['Attach Mode+' 'Decreasing'];...
    ['Attach Mode+' 'Increasing'];...
    ['Attach Mode+' 'No Change'];...
    ['Whole Cell+' 'Decreasing'];...
    ['Whole Cell+' 'Increasing'];...
    ['Whole Cell+' 'No Change'] };
for i = 1:6 % six subplots in PSTH and z score
    %---------PSTH-----
    figure(fH(1))
    subplotH = subplot(2,3,i);
    dataNoTemp = intersect(conditionIdx{i,1}, conditionIdx{i,2});
    ii = 0;
    for dataNo = dataNoTemp
        ii = ii +1;
        dataForPSTH(ii).spikeCountsInBin = rasterData(dataNo).spikeCountsInBin;
        coveredBinsForPSTH(ii).coveredIdx = rasterData(dataNo).binsCoveredByData;
    end
    if exist('dataForPSTH','var')% skipped if there is no recording for given response type and recording mode
        periHist(dataForPSTH, coveredBinsForPSTH, binEdges, subplotH);
    end
    title(titlesSubFig{i},'FontName',fontName,'FontWeight','Bold','FontSize',fontSize);
    %-------Z Score-----
    figure(fH(2))
    subplotH = subplot(2,3,i);
    ii = 0;
    for dataNo = dataNoTemp
        ii = ii +1;
        dataForZScore(ii).zScore = rasterData(dataNo).zScore;
    end
    if exist('dataForPSTH','var')% skipped if there is no recording for given response type and recording mode
        zScorePeriHist(dataForZScore, coveredBinsForPSTH, binEdges, subplotH);
    end
    title(titlesSubFig{i},'FontName',fontName,'FontWeight','Bold','FontSize',fontSize);
    clear dataForPSTH coveredBinsForPSTH dataForZScore
    
    hold on
    
    %-----Creates a matrix for Kruskawallis analysis.
    for dataNo = dataNoTemp% Episodes of given condition
        baseLineIdx = (binEdges(1:end-1) < 0) & rasterData(dataNo).binsCoveredByData;
        baseLineIdx = find(baseLineIdx == true);
        baseLineZ = rasterData(dataNo).zScore(:,baseLineIdx);
        baseLineZ = baseLineZ(:);
        afterStimIdx = (binEdges(1:end-1) >= 0) & rasterData(dataNo).binsCoveredByData;
        afterStimIdx = find(afterStimIdx == true);
        kwSubTable = nan (length(baseLineZ),size(afterStimIdx,2)+1);%prestimulation in a column + a column for each bin after stimulation
        kwSubTable(:,1) = baseLineZ;%baseline at first coloumn
        kwSubTable(1,2:end) = rasterData(dataNo).zScore(:,afterStimIdx);
        if ~exist('kwTableAll','var')% initiates a z score matrix for kruskal wallis analysis. Once its initiated this part isn't executed any more.
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
    %-----kruskawallis statistics and marks significantly different bins
    if exist('kwTableAll','var')% when no data for given condition, data for kruskawallis doesn't exist
        kwTableAll = [kwTableAll; kwTableAll];% this part will be removed when the number of samples is increased
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
    clear kwTableAll kwSubTable
end% for i = 1:6

% ---------common x and y label for all subplots in PSTH figure
figure(fH(1))
axes('Parent', fH(1), ... 
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
% -------------common x and y label for all subplots in z Score figure
figure(fH(2))
axes('Parent', fH(2), ... 
    'Units', 'normalized', 'Position', [0, 0, 1, 1], ...
    'Visible', 'off', 'NextPlot', 'add');
text(0.05,0.5, 'z Score','FontName',fontName,'FontSize',fontSize, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top',...
    'Rotation',90);
text(0.5,0, 'Time (s)','FontName',fontName,'FontSize',fontSize, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom');
