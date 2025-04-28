%% Trim raw BC Data files
clc; close all; clear all;

%% Set reference file
fileName = input('Input raw file name: ', 's');
BC_list = 'Master_BC_List.xlsx';

%% Import Data
initBCs = readtable(fileName, 'ReadRowNames', true, 'PreserveVariableNames', true);
Master_BC_List = readtable(BC_list);

[inInitBCs, ~] = ismember(Master_BC_List.BC_Sequence, initBCs.Properties.RowNames);
emptyRow = cell(1,size(initBCs,2));
emptyRow(1,:) = {0};
finalBCs = initBCs;

for BCpres = 1:length(inInitBCs)
    if ~inInitBCs(BCpres)
        rowToAdd = emptyRow;
        rowToAdd = cell2table(rowToAdd, ...
            'VariableNames',initBCs.Properties.VariableNames, ...
            'RowNames', Master_BC_List.BC_Sequence(BCpres));
        finalBCs = [finalBCs;rowToAdd];
    end
end

outputData = finalBCs(Master_BC_List.BC_Sequence,:);

%% STATS ABOUT BCS REMOVED
initCounts = sum(initBCs);
finalCounts = sum(outputData);

percentTrimmed = 1 - finalCounts./initCounts;

% other parameters to calculate:
% - how many BCs were removed
% - Distance to nearest BC that was kept











%% Write merged file to Excel
if strcmp(fileName(end-8:end-5), '_raw')
    outFile = [fileName(1:end-9), '_trim.xlsx'];
else
    outFile = [fileName(1:end-5), '_trim.xlsx'];
end

slashPres = strfind(fileName, '/');
if ~isempty(slashPres)
    outFile = [outFile((slashPres(end)+1):end)];
end

writetable(outputData, outFile,'WriteRowNames',true);