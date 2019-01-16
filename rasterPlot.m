function rasterPlot(rasterDataName,varargin)
% function rasterPlot(rasterDataName,varargin)
% This function uses rasterData created by main.m and plots a raster plot
% centered at the begining of the stimulation. When the optinal inputs are
% given, the output figure is formated accordingly. Otherwise, it is
% formated by MATLAB's default format. Each neuron is color coded. Attach
% mode recordings are plotted first.
% INPUT
%  rasterDataName = When input is rasterData loads data from the workspace.
%      When input 'rasterData.mat' loads data from saved file
% OPTIONAL INPUTS
%  'orderOfEpisodes' = When it is not given, all of the episodes are drawn 
%       starting from attach mode recordings. When a vector of episode 
%       numbers are given, the figure is drawn by this order.
%  'journalName' = the output figure is formated according to the given
%       journal format
%  'width' = width of the figure in cm. Overwrites the journal format
%  'height' = height of the figure in cm. Overwrites the journal format
%  'fontName' = defines the font name used for axes and texts on the
%       figure. Overwrites the journal format
%  'fontSize' = defines the font size used for axes and texts on the
%       figure. Overwrites the journal format

expectedJournals = {'JNeurophys'};
defaultJournal = 'JNeurophys';
defaultFontName = 'Helvetica';
defaultFontSize = 11;

isValidScalarNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

p = inputParser;
addParameter( p,'journalName',defaultJournal,@(x) any(validatestring(x,expectedJournals)) );
addParameter( p,'width',17.6, isValidScalarNum );
addParameter( p,'height',7.4976, isValidScalarNum );
addParameter( p,'fontName',defaultFontName,@(x) any(validatestring(x,listfonts)) );
addParameter( p,'fontSize',defaultFontSize, isValidScalarNum );
if strcmpi (rasterDataName,'rasterData.mat')% Loads rasterData file if the function is called for rasterData.mat. Otherwise, it can be called for a variable named rasterData
    load('rasterData.mat','rasterData')
    rasterData = rasterData(2:end);
else
    rasterData = rasterDataName(2:end);
end
addParameter( p,'orderOfEpisodes',0,@(x) isnumeric(x) && all(round(x)==x)  && length(unique(x)) <= length(rasterData) && length(x) == length(unique(x)))
parse(p,varargin{:});

if nargin > 1 && strcmpi( p.Results.journalName,'JNeurophys' )
    fH = figure;
    fH.PaperUnits = 'centimeters';
    fH.PaperOrientation ='portrait';
    fH.PaperPosition = [0 0 17.6 7.4976];
    fH.PaperPositionMode = 'manual';
    fontName = 'Arial';
    fontSize = 10;
elseif nargin > 1
    fH = figure;
    fH.PaperUnits = 'centimeters';
    fH.PaperOrientation ='portrait';
    fH.PaperPosition = [0 0 p.Results.width p.Results.height];
    fH.PaperPositionMode = 'manual';
    fontName = p.Results.fontName;
    fontSize = p.Results.fontSize;
else
    figure;% creates an unformated figure when the function is called without additional inputs
    fontName = p.Results.fontName;
    fontSize = p.Results.fontSize;
end

%%

randColorIdx = randperm(length(rasterData));
randColors = hsv(length(rasterData));
randColors = randColors(randColorIdx,:);% hsv colors are shuffled to increase the contrast between neuron colors

attachIdx = arrayfun(@(x) strcmp(x.recordingType,'Attach Mode'),rasterData);
attachIdx = find(attachIdx == 1);
wholeIdx = arrayfun(@(x) strcmp(x.recordingType,'Whole Cell'),rasterData);
wholeIdx = find(wholeIdx == 1);

hold on
yAxNo = 0;
if p.Results.orderOfEpisodes == 0
    o = [attachIdx wholeIdx];%attach mode recordings are drawn first
else
    o = p.Results.orderOfepisodes;
end
for i = o
    for ii = 1: length(rasterData(i).EpRasterData)
        if ~isempty(rasterData(i).EpRasterData(ii).spikeTime)
            for iii = 1:length(rasterData(i).EpRasterData(ii).spikeTime)
                periSpike = rasterData(i).EpRasterData(ii).spikeTime(iii) - rasterData(i).stimTime;
                plot([periSpike periSpike],[yAxNo (yAxNo +1)],...% Each spike is ploted as a vertical line
                    'Color', randColors(i,:));
            end
            yAxNo = yAxNo +1;
        end
        
    end
    
end
ylabel('Trace Number','FontName',fontName,'FontSize',fontSize)
xlabel('sec','FontName',fontName,'FontSize',fontSize)
title('Raster Plot','FontName',fontName,'FontSize',fontSize)
set(gca,'FontName',fontName,'fontsize',fontSize)
axis tight