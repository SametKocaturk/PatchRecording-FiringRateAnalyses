function rasterData = otherInfos(rasterData)

defStimTime = 1000;
defstimDur = 1000;

stimTime = nan;
while stimTime > rasterData.traceDur || isnan(stimTime) || (stimTime < 0)
    stimTime = str2double( inputdlg('Enter stimulation time in ms', 'Stimulation Time', [1 50],{num2str(defStimTime)}) )/1000;
    if isempty(stimTime)
        error('Terminated by user')
    end
end
stimDur = nan;
while (stimDur/1000 + stimTime >  rasterData.traceDur) || isnan(stimDur) || (stimDur < 0)% stimDur + stimTime cannot be bigger than recording duration
    stimDur = str2double( inputdlg('Enter stimulation duration in ms', 'Stimulation Duration', [1 50],{num2str(defstimDur)}) );
    if isempty(stimDur)
        error('Terminated by user')
    end
end

rasterData(end).stimTime = stimTime;
rasterData(end).stimDur = stimDur;
