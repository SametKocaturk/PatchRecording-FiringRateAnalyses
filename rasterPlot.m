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


load rasterData
rasterData = rasterData(2:end);
randColors = jet(length(rasterData));

attachIdx = arrayfun(@(x) strcmp(x.recordingType,'Attach Mode'),rasterData);
attachIdx = find(attachIdx == 1);
wholeIdx = arrayfun(@(x) strcmp(x.recordingType,'Whole Cell'),rasterData);
wholeIdx = find(wholeIdx == 1);

hold on
yAxNo = 0;
for i = [attachIdx wholeIdx]
    for ii = 1: length(rasterData(i).EpRasterData)
        if ~isempty(rasterData(i).EpRasterData(ii).spikeTime)
            for iii = 1:length(rasterData(i).EpRasterData(ii).spikeTime)
                periSpike = rasterData(i).EpRasterData(ii).spikeTime(iii) - rasterData(i).stimTime;
                plot([periSpike periSpike],[yAxNo (yAxNo +1)],...% Each spike is ploted as a vertical line(first trace at [0 1]
                    'Color', randColors(i,:));
            end
            yAxNo = yAxNo +1;
        end
        
    end
    
end
ylabel('Trace Number','FontName',fontName,'FontSize',fontSize)
xlabel('sec','FontName',fontName,'FontSize',fontSize)
title('Raster Plot','FontName',fontName,'FontSize',fontSize)

xAxisLabels = get(gca,'XTickLabel');
set(gca,'XTickLabel',xAxisLabels,'FontName',fontName,'fontsize',fontSize)
yAxisLabels = get(gca,'YTickLabel');
set(gca,'YTickLabel',yAxisLabels,'FontName',fontName,'fontsize',fontSize)