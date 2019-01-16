function [spikeCountsInBin,...
    zScore,...
    avgBaseline,...
    stdBaseline,...
    binsCoveredByTrace,...
    binEdges] =...
    binningZscoring...
    (binSize,EpRasterData,stimTime,traceDur)
% function [spikeCountsInBin, zScore,avgBaseline,stdBaseline,
% binsCoveredByTrace,binEdges] = binningZscoring (binSize,EpRasterData,
% stimTime,traceDur)
% This function calculates the PSTH and z score from spike times.
% INPUT
%   binSize = Binning size in seconds
%   EpRasterData = A struct with 'spikeTime' field and each element is an
%    trace of the episode
%   stimTime = stimulation time in seconds
%   traceDur = trace duration in seconds
% OUTPUT
%   spikeCounts = a vector containing the number of spikes in each bin. It
%    is oversized [-10 13 seconds] in case a long trace is given
%   zScore = z score calculated from spikeCounts. Baseline is up to the
%    stimulation.
%   avgBaseline = mean spike counts of baseline bins(up to stimulation)
%   stdBaseline = standard deviation of baseline bins. If there is only one
%    bin in baseline, the standard deviation is one.
%   binsCoveredByTrace = each row is a trace excluding empty elements of
%    EpRasterData. Each column corresponds a bin where the value is 1 if it
%    is spanned by the trace.
%   binEdges = starting and end time point of each bin in seconds
%
negBinEdges = flip(-binSize:-binSize:-10);
positiveBinEdges = 0:binSize:13;
binEdges = [negBinEdges positiveBinEdges];% stimulation time is 0. Creates longer than needed in order to cover different length of trials. It will later trancated down to a coverage by trials

spikeCountsInBin = nan(length(EpRasterData), length(binEdges)-1);
binsCoveredByTrace = zeros(length(EpRasterData), length(binEdges));
traceHistNo = 0;
for traceNo = 1: length(EpRasterData)
    if ~isempty(EpRasterData(traceNo).spikeTime)% skips discarded trials
        traceHistNo = traceHistNo +1; % because of the discarded trials, the number of histogram is different than the number of trials
        trSpikeTimes = EpRasterData(traceNo).spikeTime;
        spikeCountsInBin(traceHistNo,:) = histcounts(trSpikeTimes - stimTime,binEdges);
        binsCoveredByTrace(traceHistNo,:) = ( (binEdges > (-stimTime-binSize)) ) ...
            & ( (binEdges < (traceDur - stimTime)) );
    end
end

coveredIdx = sum(binsCoveredByTrace) ~= 0;
baseLineIdx = (binEdges < 0) & coveredIdx;%bins up to the stimulation time

spikeCountsInBin = nansum(spikeCountsInBin);
avgBaseline = mean( spikeCountsInBin(baseLineIdx(1:end-1)) );
stdBaseline = std(spikeCountsInBin(baseLineIdx(1:end-1)) );
if stdBaseline == 0; stdBaseline =1;end %if the baseline consits of one bin sd = 0 so sd converted to 1
zScore = (spikeCountsInBin - avgBaseline) ./ stdBaseline;
zScore(~coveredIdx(1:end-1)) = 0;