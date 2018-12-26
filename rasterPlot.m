% Plots a raster plot from rasterData.mat. Each neuron is color coded.
close all;clear all;clc

JNeurophys = figure;
width = 17.6;%cm
height = width*0.4260;
JNeurophys.PaperUnits = 'centimeters';
JNeurophys.PaperOrientation ='portrait';
JNeurophys.PaperPosition = [0 0 width height];
JNeurophys.PaperPositionMode = 'manual';

fontName = 'Arial';
fontSize = 10;

%%
load rasterData
rasterData = rasterData(2:end);
randColorIdx = randperm(length(rasterData));
randColors = hsv(length(rasterData));
randColors = randColors(randColorIdx,:);% hsv colors are shuffled to increase the contrast between neuron colors

attachIdx = arrayfun(@(x) strcmp(x.recordingType,'Attach Mode'),rasterData);
attachIdx = find(attachIdx == 1);
wholeIdx = arrayfun(@(x) strcmp(x.recordingType,'Whole Cell'),rasterData);
wholeIdx = find(wholeIdx == 1);

hold on
yAxNo = 0;
for i = [attachIdx wholeIdx]%attach mode recordings are drawn first
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