function [] = makeHorzBarChart(rawdata, labels_cell, toAnnotate, RFCutoff, fileName, fileFolder)

    rawdata = flip(rawdata,2);
    labels_cell = flip(labels_cell);
% Inputs: 
% - data should be normalized so that columns add up to 1
% - labels should be a cell array with column labels
% - file names should be strings with location to save files
    nBCs = size(rawdata,1);
    map = distinguishable_colors(nBCs);
    data = zeros(size(rawdata));
    countCutoff = 5;

    for i = 1:size(rawdata, 1)
        for j = 1:size(rawdata, 2)
            if rawdata(i,j) > countCutoff
                data(i,j) = rawdata(i,j);
            end
        end
    end

    data = data ./ sum(data);


    f = figure();
    nBars = length(labels_cell);
    labels = categorical(labels_cell); labels = reordercats(labels, cellstr(labels));
    bars1 = barh(labels,data','stacked','BarWidth',0.8,'FaceColor','flat');
    xlim([0 1]); xlabel('Barcode Proportions'); %ylabel('Hamster Number');
    title(fileName);

    for i = 1:nBCs
        colors = zeros(nBars,3);
        colors(:,1) = colors(:,1) + map(i,1);
        colors(:,2) = colors(:,2) + map(i,2);
        colors(:,3) = colors(:,3) + map(i,3);
        bars1(i).CData = colors;
    end
    
    BCCounts = countBCs(rawdata, RFCutoff, countCutoff);
    annotateHorzBarChart(f, BCCounts, toAnnotate);

    figHeight = 75*nBars+165;
    f.Position = [40 40 1440 figHeight];

    fontsize(gcf, scale=1.5)

    saveas(f, [fileFolder, '/', fileName, '.fig']);
    saveas(f, [fileFolder, '/', fileName, '.png']);
end

function [] = annotateHorzBarChart(f, BCCounts, totalAnnotate)
    
    for i = 1:length(BCCounts)
        text(1.019,i,num2str(BCCounts(i)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'Rotation',270);
    end

    % if ~isempty(totalAnnotate)
    %     text(1.053,(nBars+1)/2,[num2str(totalAnnotate),' Unique BCs'], ...
    %         'HorizontalAlignment','center', ...
    %         'VerticalAlignment','middle', ...
    %         'Rotation',270, ...
    %         'FontWeight','bold');
    % end
end

% function [] = annotateHorzBarChart(f, BCCounts, toAnnotate)
%     for i = 1:length(BCCounts)
%         if toAnnotate(i) > 0
%             annot = annotation('textbox',[.3 .8 .01 .01],'String',num2str(BCCounts(i)));
%             annot.Parent = f.CurrentAxes;
%             annot.Margin = 4;
%             annot.BackgroundColor = [1 1 1];
%             annot.HorizontalAlignment = 'center';
%             annot.VerticalAlignment = 'middle';
%             annot.Position = [0.93 i-0.2 0.05 0.4];
%         end
%     end
% end


function [counts] = countBCs(rawdata, RFCutoff, countCutoff)    
    dataNorm = rawdata ./ sum(rawdata);
    
    counts = zeros(1,size(rawdata,2));

    for j = 1:size(rawdata,2)
        counts(j) = sum((rawdata(:,j)>countCutoff) & (dataNorm(:,j)>RFCutoff));
    end
end