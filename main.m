close all;clear all; clc
% This script searches for 'yymmdddata*' formated .mat electrophysiology
% recording in the current directory. The file is struct array with fields 
% interval : sampling period and 
% values : rows(data)xcolumns(channels)xpages(traces).
%
% The output 'rasterData.mat' file contains a struct array named
% 'rasterData'. Each element of the struct contains an episode. The fields:
% DataName: original recording's nameclose all;clear all; clc
% This script searches for 'yymmdddata*' formated .mat electrophysiology
% recording in the current directory. The file is struct array with fields 
% interval : sampling period and 
% values : rows(data)xcolumns(channels)xpages(traces).
%
% The output 'rasterData.mat' file contains a struct array named
% 'rasterData'. Each element of the struct contains an episode. The fields:
% DataName: original recording's name
% EpRasterData: Is a struct. Each element is a trace of the episode.
%    Contains field 'spikeTime': coloumn vector of spike times.
% recordingType: string of recording mode(Attach Mode or Whole Cell)
% stimTime: stimulation time in sec
% stimDur : stimulation duration in ms
% traceDur: trace duration in sec

if exist('rasterData.mat','file')
    load rasterData
else
    rasterData = struct('DataName',[], 'EpRasterData',[], 'recordingType',[], 'stimTime',[], 'stimDur',[], 'traceDur',[]);
    save('rasterData.mat','rasterData')
end

recordingList = dir('*data*.mat');% search for exported data with ‘.mat’ extention
validDataNameIndxs = [];
for i = 1:size(recordingList,1)
    if ~isempty( regexp(recordingList(i).name,'\d{4}data\d', 'once') )% Exported data name format: ['yymmdddata' description]
        validDataNameIndxs = [validDataNameIndxs i];
    end
end
recordingList = recordingList(validDataNameIndxs);

opts.Interpreter = 'tex';

for i =1:length(recordingList)
    redo = 'Redo';
    if any(  arrayfun( @(x) strcmp( x.DataName , recordingList(i).name ), rasterData  )     )
        opts.Default = 'Skip';
        redo = questdlg(['\fontsize{15}\bf ' recordingList(i).name ' \rmwas analyised. Do you want to redo?'],...
            'Respike detection?',...
            'Skip','Redo',opts);
        if ~(     strcmp(redo,'Skip') || strcmp(redo,'Redo')     )
            error('Wrong selection')
        end
        rasterDataNo = find(  arrayfun(@(x) strcmp(x.DataName, recordingList(i).name), rasterData), 1,'first'    );% if data was previously analysed, its location is detected
    else
        rasterDataNo = length(rasterData) + 1;% if the data is new, it will be added at the end
    end
    
    if strcmp(redo,'Redo') %redo = 'Redo' by default so that this part evaluated if the data hasn't analyzed or will be reanalyzed
        fileFullName = fullfile(recordingList(i).folder, recordingList(i).name);
        tempdata = load(fileFullName);
        dataName = fieldnames(tempdata); dataName = dataName{1};
        data = eval([ 'tempdata.' dataName]);% tempdata.Data1_wave_data.(values and other info)
        
        [r,~,p] = size(data.values);
        data.values = reshape (data.values,r,p);
        t = (double(0:data.points-1) * data.interval)';
        t = repmat(t,1,p);
        data.time = t;
        
        [EpRasterData,recordingType] = spikeVisualize(data);
        
        
        %% Fill rasterData table
        rasterData(rasterDataNo).DataName = recordingList(i).name;
        rasterData(rasterDataNo).EpRasterData = EpRasterData;
        rasterData(rasterDataNo).recordingType = recordingType;
        rasterData(rasterDataNo).traceDur = data.time(end);
        rasterData(rasterDataNo) = otherInfos(rasterData(rasterDataNo));
        
        save('rasterData.mat','rasterData')        
    end
end
close all
opts.Default = 'Yes';
raster = questdlg('\fontsize{15} Do you want to see raster plot?',...
    'Raster Plot?',...
    'Yes','No',opts);
if strcmp(raster,'Yes')
    rasterPlot;
end


% EpRasterData: Is a struct. Each element is a trace of the episode.
%    Contains field 'spikeTime': coloumn vector of spike times.
% recordingType: string of recording mode(Attach Mode or Whole Cell)
% stimTime: stimulation time in sec
% stimDur : stimulation duration in ms
% traceDur: trace duration in sec

if exist('rasterData.mat','file')
    load rasterData
else
    rasterData = struct('DataName',[], 'EpRasterData',[], 'recordingType',[], 'stimTime',[], 'stimDur',[], 'traceDur',[]);
    save('rasterData.mat','rasterData')
end

recordingList = dir('*data*.mat');% search for exported data with ‘.mat’ extention
validDataNameIndxs = [];
for i = 1:size(recordingList,1)
    if ~isempty( regexp(recordingList(i).name,'\d{4}data\d', 'once') )% Exported data name format: ['yymmdddata' description]
        validDataNameIndxs = [validDataNameIndxs i];
    end
end
recordingList = recordingList(validDataNameIndxs);

opts.Interpreter = 'tex';

for i =1:length(recordingList)
    redo = 'Redo';
    if any(  arrayfun( @(x) strcmp( x.DataName , recordingList(i).name ), rasterData  )     )
        opts.Default = 'Skip';
        redo = questdlg(['\fontsize{15}\bf ' recordingList(i).name ' \rmwas analyised. Do you want to redo?'],...
            'Respike detection?',...
            'Skip','Redo',opts);
        if ~(     strcmp(redo,'Skip') || strcmp(redo,'Redo')     )
            error('Wrong selection')
        end
        rasterDataNo = find(  arrayfun(@(x) strcmp(x.DataName, recordingList(i).name), rasterData), 1,'first'    );% if data was previously analysed, its location is detected
    else
        rasterDataNo = length(rasterData) + 1;% if the data is new, it will be added at the end
    end
    
    if strcmp(redo,'Redo') %redo = 'Redo' by default so that this part evaluated if the data hasn't analyzed or will be reanalyzed
        fileFullName = fullfile(recordingList(i).folder, recordingList(i).name);
        tempdata = load(fileFullName);
        dataName = fieldnames(tempdata); dataName = dataName{1};
        data = eval([ 'tempdata.' dataName]);% tempdata.Data1_wave_data.(values and other info)
        
        [r,~,p] = size(data.values);
        data.values = reshape (data.values,r,p);
        t = (double(0:data.points-1) * data.interval)';
        t = repmat(t,1,p);
        data.time = t;
        
        [EpRasterData,recordingType] = spikeVisualize(data);
        
        
        %% Fill rasterData table
        rasterData(rasterDataNo).DataName = recordingList(i).name;
        rasterData(rasterDataNo).EpRasterData = EpRasterData;
        rasterData(rasterDataNo).recordingType = recordingType;
        rasterData(rasterDataNo).traceDur = data.time(end);
        rasterData(rasterDataNo) = otherInfos(rasterData(rasterDataNo));
        
        save('rasterData.mat','rasterData')        
    end
end
close all
opts.Default = 'Yes';
raster = questdlg('\fontsize{15} Do you want to see raster plot?',...
    'Raster Plot?',...
    'Yes','No',opts);
if strcmp(raster,'Yes')
    rasterPlot;
end

