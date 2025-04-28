clc; clear all; close all;

%% Load data
% Read in BC counts
file = 'Figure S3 - MvF/FigS3-MvF.xlsx';
data_raw = readtable(file, 'ReadRowNames', true, 'PreserveVariableNames', true);
fprintf('Data read in\n')
data_norm = data_raw ./ sum(data_raw);

MADCutoff = 0.00056;
countCutoff = 10;

nBCs = size(data_raw,1);

% data_norm(or(data_norm<MADCutoff, data_raw<countCutoff)) = 0;
% data_norm = data_norm ./ sum(data_norm);

%% Organize data into structs
baseHamster.DonCon = ''; % donor or contact
baseHamster.sex = ''; % male or female
hamster.dur = '8hr'; % duration of exposure for contacts
hamster.route = 'AB'; % route of exposure for contacts
% baseHamster.time = ''; % For donors, hpi of collection. For contacts, dpe of harvest
baseHamster.Exp = 'T42';

baseHamster.NT.raw = zeros(nBCs,1); % raw reads for each BC in nasal turbinate
baseHamster.Trach.raw = zeros(nBCs,1); % raw reads for each BC in trachea
baseHamster.WL.raw = zeros(nBCs,1); % raw reads for each BC in lungs
baseHamster.LL.raw = zeros(nBCs,1); % raw reads for each BC in left lung lobe
baseHamster.RL1.raw = zeros(nBCs,1); % raw reads for each BC in right apical lung lobe
baseHamster.RL2.raw = zeros(nBCs,1); % raw reads for each BC in right middle lung lobe
baseHamster.RL3.raw = zeros(nBCs,1); % raw reads for each BC in right caudal lung lobe
baseHamster.RL4.raw = zeros(nBCs,1); % raw reads for each BC in infracardial lung lobe

T42 = {baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, ...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, ...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, ...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster};
% First 15 cells are spacers so that index = hamster number

samples = data_raw.Properties.VariableNames;

for i = 1:length(samples)
    % Determine hamster within experiment
    hamsterNumIndex = strfind(samples{i},'T42_H');
    hamsterNum = str2num(samples{i}(hamsterNumIndex+5:hamsterNumIndex+6));
    T42{hamsterNum}.animalNum = hamsterNum;

    % discern collection time
    %  note to self: make collection time read in more robust by
    %  checking if a value has already been assigned to it
    if strfind(samples{i},'Don')
        T42{hamsterNum}.DonCon = 'Donor';
        T42{hamsterNum}.time = '32h';
    elseif strfind(samples{i},'Con')
        T42{hamsterNum}.DonCon = 'Contact';
        T42{hamsterNum}.time = '3dpe';
    else
        error('Donor or contact not properly specified');
    end
    
    if strfind(samples{i},'Male')
        T42{hamsterNum}.sex = 'M';
    elseif strfind(samples{i},'Fem')
        T42{hamsterNum}.sex = 'Fem';
    else
        error('Sex not properly specified');
    end

    % Discern tissue
    if strfind(samples{i},'NT')
        T42{hamsterNum}.NT.raw = T42{hamsterNum}.NT.raw + data_raw{:,i};
    elseif strfind(samples{i},'Trach')
        T42{hamsterNum}.Trach.raw = T42{hamsterNum}.Trach.raw + data_raw{:,i};
    elseif strfind(samples{i},'WL')
        T42{hamsterNum}.WL.raw = T42{hamsterNum}.WL.raw + data_raw{:,i};
    elseif strfind(samples{i},'LL')
        T42{hamsterNum}.LL.raw = T42{hamsterNum}.LL.raw + data_raw{:,i};
    elseif strfind(samples{i},'RL1')
        T42{hamsterNum}.RL1.raw = T42{hamsterNum}.RL1.raw + data_raw{:,i};
    elseif strfind(samples{i},'RL2')
        T42{hamsterNum}.RL2.raw = T42{hamsterNum}.RL2.raw + data_raw{:,i};
    elseif strfind(samples{i},'RL3')
        T42{hamsterNum}.RL3.raw = T42{hamsterNum}.RL3.raw + data_raw{:,i};
    elseif strfind(samples{i},'RL4')
        T42{hamsterNum}.RL4.raw = T42{hamsterNum}.RL4.raw + data_raw{:,i};
    else
        error('Tissue of collection not properly specified');
    end
    
end

% for i = 1:length(T42)
%     if sum(T42{i}.WL.raw) == 0
%         T42{i}.WL.raw = T42{i}.LL.raw + T42{i}.RL1.raw + ...
%             T42{i}.RL2.raw + T42{i}.RL3.raw + T42{i}.RL4.raw;
%     end
% end

fprintf('Raw data organized into structs\n')

%% Add info about which BCs are transmitted, how many BCs are transmitted
for i = 1:length(T42)
    T42{i}.NT.norm = T42{i}.NT.raw ./ sum(T42{i}.NT.raw);
    T42{i}.NT.trans = and((T42{i}.NT.raw>countCutoff), (T42{i}.NT.norm>MADCutoff));
    T42{i}.NT.raw = T42{i}.NT.raw .* T42{i}.NT.trans;
    T42{i}.NT.norm = T42{i}.NT.raw ./ sum(T42{i}.NT.raw);
    T42{i}.NT.nBCs = sum(T42{i}.NT.trans);

    T42{i}.Trach.norm = T42{i}.Trach.raw ./ sum(T42{i}.Trach.raw);
    T42{i}.Trach.trans = and((T42{i}.Trach.raw>countCutoff), (T42{i}.Trach.norm>MADCutoff));
    T42{i}.Trach.raw = T42{i}.Trach.raw .* T42{i}.Trach.trans;
    T42{i}.Trach.norm = T42{i}.Trach.raw ./ sum(T42{i}.Trach.raw);
    T42{i}.Trach.nBCs = sum(T42{i}.Trach.trans);

    T42{i}.WL.norm = T42{i}.WL.raw ./ sum(T42{i}.WL.raw);
    T42{i}.WL.trans = and((T42{i}.WL.raw>countCutoff), (T42{i}.WL.norm>MADCutoff));
    T42{i}.WL.raw = T42{i}.WL.raw .* T42{i}.WL.trans;
    T42{i}.WL.norm = T42{i}.WL.raw ./ sum(T42{i}.WL.raw);
    T42{i}.WL.nBCs = sum(T42{i}.WL.trans);

    T42{i}.totalBCs = sum(max([T42{i}.NT.trans, T42{i}.Trach.trans, T42{i}.WL.trans],[],2));
end

fprintf('Total BCs present determined\n')

%% Make stacked bar charts
% Set output folder
folder = 'Figure S3 - MvF/figs';
colmap = distinguishable_colors(nBCs);

makeFigs = input('Re-make stacked bar charts? ', 's');
labels = {'NT','Trachea','Lungs'};

if upper(makeFigs) == 'Y'
    for i = [7,19]% 1:length(T42)
        ylabelstring = [T42{i}.sex, 'ale ', T42{i}.DonCon];
        
        makeHorzBarChartVenn(T42{i}, labels, T42{i}.totalBCs, ylabelstring, folder);
    end
end