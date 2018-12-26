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
    stdFiltDataValue = stdfunction [EPRaster,recordingType] = spikeVisualize(data)
% function [EPRaster,recordingType] = spikeVisualize(data)
% This function creates a plot of whole episode then an interactive trace by
% trace figure.  It calls the peakTimeDetection function to find spike 
% peaks and then marks the spikes on the raw data.In trace by trace mode, 
% you can toggle to discard mode by keyboard input. You can manually peak a
% threshold by clicking a trace.
% INPUT
%  data = struct containing exported data
% OUTPUT
%  EPRaster = struct with field spikeTime. Each element stands for a trace.
%   spikeTime containing spike times as a coloumn vector
%  recordingType = a string returned by peakTimeDetection function('Whole
%   Cell' or 'Attach Mode')
global selectedPlot  selectionFlag doneFlag nextFlag
opts.Interpreter = 'tex';

%% -----------------------------WHOLE EPISODE------------------------------
% For whole Episode. General threshold. Spikes. Raster Table
if ~exist('fh1','var'),hF1 = figure;end %If user skips an episode which was analysed previously
hF1.Units = 'normalized';
hF1.OuterPosition = [0.1 0.1 0.9 0.9];
hold off
plot(data.time,data.values);% plot raw data
axh = gca;
[spikeTime, spikeIdx, recordingType, spikeCountMismatchFlag, amplitudeMisMatchFlag] = peakTimeDetection(data.values,data.time);
hold on
for tr = 1:size(data.values,2)% the number of trials
    plot(data.time(spikeIdx{tr},tr),data.values(spikeIdx{tr},tr),'*r' )
    EPRaster(tr).spikeTime = spikeTime{tr};
end

if any(spikeCountMismatchFlag)
    opts.Default = 'Gotcha';
    spikeCountWarning = questdlg('\fontsize{15} Spike count seems problematic in at least one trace?',...
        'Warning?',...
        'Gotcha',opts);
    if ~ strcmp(spikeCountWarning,'Gotcha')
        error('Wrong selection')
    end
end

if any(amplitudeMisMatchFlag)
    opts.Default = 'Gotcha';
    amplitudeMisMatchWarn = questdlg('\fontsize{15} Spike amplitudes seem problematic in at least one trace?',...
        'Warning?',...
        'Gotcha',opts);% if 'Yes', it can escape from the while loop.
    if ~ strcmp(amplitudeMisMatchWarn,'Gotcha')
        error('Wrong selection')
    end
end
%% -----------------------------TRACE BY TRACE-----------------------------
opts.Default = 'No';
traceBytrace = questdlg('\fontsize{15} Do you want to analyze \bftrace by trace?',...
    'Redo?',...
    'Yes','No',opts);
if ~(     strcmp(traceBytrace,'Yes') || strcmp(traceBytrace,'No')     )
    error('Wrong selection')
end

nrTraces = size(data.values,2);
if strcmp(traceBytrace,'Yes')
    xlimits = axh.XLim;
    ylimits = axh.YLim;
    nrSubFigPages = fix(nrTraces / 8) + (nrTraces / 8 ~=0) - 1;%the number of pages required for all sub plots(starting from 0)
    subFigPage = 0;
    subFigSize = [2,4];% 8 subplots per page
    doneFlag = false;
    nextFlag = false;
    discardList = [];
    while doneFlag == false
        if nextFlag == true && (subFigPage+1 <= nrSubFigPages)
            subFigPage = subFigPage +1;
            discardList = [];
            nextFlag = false;
        end
        if nrTraces >= (subFigPage+1) *subFigSize(1)* subFigSize(2)%if the number of left over traces bigger than the number of subplots fitting into a page
            nrSubPlots = subFigSize(1)* subFigSize(2);
        else
            nrSubPlots = mod(nrTraces,subFigSize(1)* subFigSize(2));
        end
        
        subH = zeros(nrSubPlots,1);
        clf(hF1)
        for  i = 1:nrSubPlots
            subH(i) = subplot(subFigSize(1),subFigSize(2),i,'Tag',num2str(i),'ButtonDownFcn', @(src,evnt) subFigSelectionCallback(src,evnt));%clickable subplots
            traceNo = i + (subFigPage *subFigSize(1)* subFigSize(2));
            hold all
            plot(data.time(:,traceNo),data.values(:,traceNo),'HitTest','off')
            plot(EPRaster(traceNo).spikeTime,...
                data.values(spikeIdx{traceNo},traceNo), '*r','HitTest','off')%spike peaks
            set(gca,'xlim',xlimits)
            set(gca,'ylim',ylimits)
            if any(ismember(i,discardList))
                plot(xlimits,ylimits,'r')%draw a cross on discarded traces
            end
            if spikeCountMismatchFlag(traceNo)
                title('Spike Count Mismatch')
            end
            if amplitudeMisMatchFlag(traceNo)
                tH = get(gca,'title');
                tH.String = {tH.String, 'Spike Amplitude Problem'};% adds the warnin into the title
            end
            hold off
        end
        uicontrol('Parent',hF1,'Style','pushbutton','String','DONE','Units','normalized','Position',[0.5 0.02 0.1 0.05],'Visible','on','Enable','inactive','ButtonDownFcn',@(src,evnt) doneButtonCallback(src,evnt));
        uicontrol('Parent',hF1,'Style','pushbutton','String','NEXT','Units','normalized','Position',[0.3 0.02 0.1 0.05],'Visible','on','Enable','inactive','ButtonDownFcn',@(src,evnt) nextButtonCallback(src,evnt));
        discardModeFlag = false;
        selectionFlag = false;
        while selectionFlag == false%a subplot must be selected
            button = waitforbuttonpress; %mouse click 0, keyboard 1
            if button == 1 % toggles the discard mode when a keyboard key is pressed
                discardModeFlag = ~discardModeFlag;
            end
            for bx = 1:nrSubPlots % the boxes turns red in discard mode
                if discardModeFlag == true
                    set(subH(bx),'XColor','red')
                    set(subH(bx),'YColor','red')
                else
                    set(subH(bx),'XColor','black')
                    set(subH(bx),'YColor','black')
                end
            end
        end
        
        % One trace
        reDoThisTrace = 'Yes';
        while strcmp(reDoThisTrace,'Yes') && ~isempty(selectedPlot)
            if discardModeFlag == true % discard the trace by deleting contant
                spikeIdx{selectedPlot+(subFigPage *subFigSize(1)* subFigSize(2))} = [];
                EPRaster(selectedPlot+(subFigPage *subFigSize(1)* subFigSize(2))).spikeTime = [];
                subplot(subH(selectedPlot))
                hold on
                plot(xlimits,ylimits,'r')
                discardList = [discardList selectedPlot];
                reDoThisTrace = 'No';
            else
                clf(hF1)
                traceNo = selectedPlot + subFigPage *subFigSize(1)* subFigSize(2);
                plot(   data.time(:, traceNo),...
                    data.values(:,traceNo)  )
                hold on
                plot(EPRaster(traceNo).spikeTime,...
                    data.values(spikeIdx{traceNo},traceNo), '*r')
                
                [~,newThreshold] = ginput(1);
                [traceSpikeTime, traceSpikeIndices] = peakTimeDetection(data.values(:,traceNo),data.time(:, traceNo), newThreshold, recordingType);
                
                plot(   data.time(traceSpikeIndices{1},traceNo),...
                    data.values(traceSpikeIndices{1},traceNo),'*g' )
                
                spikeIdx{traceNo} = traceSpikeIndices{1};
                EPRaster(traceNo). spikeTime = traceSpikeTime{1};
                if ismember(selectedPlot,discardList)
                    discardList(discardList == selectedPlot) = [];% the previously discarded trace removed from the discardList when it is anaylzed again
                end
            end
             
            if discardModeFlag == false
                opts.Default = 'No';
                reDoThisTrace = questdlg('\fontsize{15} Do you want to redo this Trace?',...
                    'Redo?',...
                    'Yes','No',opts);
                if ~(     strcmp(reDoThisTrace,'Yes') || strcmp(reDoThisTrace,'No')     )
                    error('Wrong selection')
                end
            end
        end
        
        
    end% while doneFlag == false
end% if strcmp(traceBytrace, 'Yes')

close(hF1)
clear fh1

end% spikeVisualize
%% Callback Functions
function doneButtonCallback(src,evnt)
src.Value = 1;
global doneFlag selectionFlag selectedPlot
doneFlag = true; selectionFlag = true;
selectedPlot = [];%prevents from drawing the last shown trace
pause(0.1)
src.Value = 0;
end

function nextButtonCallback(src,evnt)
src.Value = 1;
global nextFlag selectionFlag selectedPlot
nextFlag = true; selectionFlag = true;
selectedPlot = [];%prevents from drawing the last shown trace
pause(0.1)
src.Value = 0;
end

function subFigSelectionCallback(src,evnt)
global selectedPlot selectionFlag
selectedPlot = str2num(get(src,'tag'));
selectionFlag = true;
end


(filtDataValue);
    
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