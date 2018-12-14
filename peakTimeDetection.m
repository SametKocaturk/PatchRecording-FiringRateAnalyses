function [spikeTime,...
    spikeIndicesAll,...
    recordingType,...
    spikeCountMismatchFlag,...
    amplitudeMisMatchFlag] = peakTimeDetection(dataValue, dataTime, newThreshold, recordingType)
% This function detects the peaks in the traces. When it's called for whole
% episode visualization, it is called with 2 inputs. Then it calls
% recordingTypeDetection and peakDetection functions. When it is called for
% trace by trace analysis, it is called with 4 inputs. Then it calls
% peakDetection function to detect spikes.
% It compares the spike counts in raw data and bandpass filtered data and
% raises countMisMatchFlag if they are not equal.
% If a spike amplitudes is too different (> 1 standard deviation),
% amplitudeMisMatchFlag is raised.
% INPUT
%  dataValue = rows containing each trace when it is called for whole
%    episode or only the trace data when it is called for trace analysis
%  dataTime = rows containing each trace when it is called for whole
%    episode or only the trace data when it is called for trace analysis
%  newThreshold = the threshold manuelly inputted in single trace analysis
%  recordingType = Sting of recording type fed when it is called for
%    single trace analysis
% OUTPUT
%  spikeTime = cell array with each element containing spike times in each
%    trace
%  spikeIndicesAll = cell array with each element containing spike indices
%    in each trace
%  recordingType = string of recording type ('Whole Cell' or 'Attach Mode')
%    detected by recordingType function
%  spikeCountMismatchFlag = logical 1 when the number of spikes in the
%    filterred and raw data do not match
%  amplitudeMisMatchFlag = logical 1 when a spike amplitude is more than 1
%    standard deviation

if nargin == 2 %when it is called for whole episode but not trace by trace visual
    disp('WAIT!!!')
    disp('Loading Filter')
    Fs = 1 / mean(mean(diff(dataTime)));
    cutOffFreq = [10 3000];
    if Fs/2 < cutOffFreq(2)
        filtDataValue = highpass(dataValue, cutOffFreq(1), Fs);%filters baseline drifts and noise. The # of spikes will be compared.
    else
        filtDataValue = bandpass(dataValue, cutOffFreq, Fs);%filters baseline drifts and noise. The # of spikes will be compared.
    end
    disp('Filter Loaded')
    recordingType = recordingTypeDetection(dataValue);
    if strcmp(recordingType, 'Attach Mode')%when recording type is in attach mode, the signal is inverted
        dataValue = -dataValue;
        filtDataValue = -filtDataValue;
    end
    
    avgDataValue = mean(dataValue);
    stdDataValue = std(dataValue);
    avgFiltDataValue = mean(filtDataValue);
    stdFiltDataValue = std(filtDataValue);
    
    threshold = avgDataValue + 10*stdDataValue;% 10 standard deviation for threshold
    filtThreshold = avgFiltDataValue + 10*stdFiltDataValue;
    
    [filtSpikeTime, ~] = peakDetection(filtDataValue, dataTime, filtThreshold);
    [spikeTime, spikeIndicesAll] = peakDetection(dataValue, dataTime, threshold);
    
    if nargout == 5 %when it is called for whole episode but not trace by trace analysis
        amplitudeMisMatchFlag = false(size(spikeIndicesAll));
        for tr = 1: length(spikeIndicesAll)
            spikeAmplitude = dataValue(spikeIndicesAll{tr},tr);
            meanSpikeAmplitude = mean(spikeAmplitude);
            stdSpikeAmplitude = std(spikeAmplitude);
            if any(     abs(spikeAmplitude) > ( abs(meanSpikeAmplitude) + 2*stdSpikeAmplitude )  |  abs(spikeAmplitude) < ( abs(meanSpikeAmplitude) - 2*stdSpikeAmplitude )      )
                amplitudeMisMatchFlag(tr) = true;
            end
        end
        
        spikeCountMismatchFlag = cellfun('length',filtSpikeTime) ~= cellfun('length',spikeTime);%if the number of spikes doesn't match, a warning flag is set for that trace
    end
elseif nargin ==4 && strcmp(recordingType, 'Attach Mode')%when this function called for trace by trace analysis and recordingType = Attach Mode
    dataValue = -dataValue;% In attach mode, the traces are flipped up side down
    threshold = -newThreshold;% In attach mode, the threshold is flipped up side down because peakDetection looks for max points
    [spikeTime, spikeIndicesAll] = peakDetection(dataValue, dataTime, threshold);
else %when this function is called for trace by trace analysis and recordingType = Whole Cell
    threshold = newThreshold;
    [spikeTime, spikeIndicesAll] = peakDetection(dataValue, dataTime, threshold); 
end
end % function peakTimeDetection

%%
function [spikeTime, spikeIndicesAll] = peakDetection(dataValue, dataTime, threshold)
% Detects spike peaks
spikeTime = cell(size(dataValue,2),1);
spikeIndicesAll = cell(size(dataValue,2),1);
for tr = 1 : size(dataValue,2)
    traceValues = dataValue(:,tr);
    peakIndices = find( traceValues >= threshold(tr) );
    if ~isempty(peakIndices)
        tmpPeakBeginEndIndices = find(diff(peakIndices) > 10); %find the beginning and end of peaks above the threshold
        peakBeginEndIndices = [   1; [tmpPeakBeginEndIndices;length(peakIndices)]+1   ];
        spikeIndices = [];
        for ii = 1 : length(peakBeginEndIndices)-1
            [~,spikeIdx] = max(     traceValues( peakIndices(peakBeginEndIndices(ii)) : peakIndices(peakBeginEndIndices(ii+1)-1) )    );
            spikeIdx = spikeIdx-1 + peakIndices(peakBeginEndIndices(ii));
            spikeIndices = [spikeIndices spikeIdx];
        end
        spikeTime{tr} = dataTime(spikeIndices,tr);
        spikeIndicesAll{tr} =  spikeIndices;
    end
end

end % function peakDetection
%%
function recordingType = recordingTypeDetection(dataValue)

maxValue = max(dataValue);
minValue = min(dataValue);
avgValue = mean(dataValue);

if abs(maxValue - avgValue) < abs(minValue - avgValue)
    recordingType = 'Attach Mode';
else
    recordingType = 'Whole Cell';
end

end%function recordingTypeDetection