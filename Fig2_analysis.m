clc; close all; %clear all;

%% Load data
% Read in BC counts
file = 'Figure 2, S2 - T43 contacts/Fig2_T43_Contacts.xlsx';
data_raw = readtable(file, 'ReadRowNames', true, 'PreserveVariableNames', true);
fprintf('Data read in\n')
data_norm = data_raw ./ sum(data_raw);

MADCutoff = 0.00056;
countCutoff = 10;

nBCs = size(data_raw,1);

% data_norm(or(data_norm<MADCutoff, data_raw<countCutoff)) = 0;
% data_norm = data_norm ./ sum(data_norm);

%% Organize data into structs
baseHamster.DonCon = 'Contact'; % donor or contact
baseHamster.sex = 'M'; % male or female
hamster.dur = '8h'; % duration of exposure for contacts
hamster.route = 'AB'; % route of exposure for contacts
baseHamster.time = ''; % For donors, hpi of collection. For contacts, dpe of harvest
baseHamster.Exp = 'T43';

baseHamster.NT.raw = zeros(nBCs,1); % raw reads for each BC in nasal turbinate
baseHamster.Trach.raw = zeros(nBCs,1); % raw reads for each BC in trachea
baseHamster.WL.raw = zeros(nBCs,1); % raw reads for each BC in lungs
baseHamster.LL.raw = zeros(nBCs,1); % raw reads for each BC in left lung lobe
baseHamster.RL1.raw = zeros(nBCs,1); % raw reads for each BC in right apical lung lobe
baseHamster.RL2.raw = zeros(nBCs,1); % raw reads for each BC in right middle lung lobe
baseHamster.RL3.raw = zeros(nBCs,1); % raw reads for each BC in right caudal lung lobe
baseHamster.RL4.raw = zeros(nBCs,1); % raw reads for each BC in infracardial lung lobe

T43_con = {'','','','','','','','','','','','','','','',...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, ...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, ...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster};
% First 15 cells are spacers so that index = hamster number

samples = data_raw.Properties.VariableNames;

for i = 1:length(samples)
    % Determine hamster within experiment
    hamsterNumIndex = strfind(samples{i},'T43_H');
    hamsterNum = str2num(samples{i}(hamsterNumIndex+5:hamsterNumIndex+6));
    T43_con{hamsterNum}.animalNum = hamsterNum;

    % discern collection time
    %  note to self: make collection time read in more robust by
    %  checking if a value has already been assigned to it
    if strfind(samples{i},'1dpe')
        T43_con{hamsterNum}.time = '1dpe';
    elseif strfind(samples{i},'2dpe')
        T43_con{hamsterNum}.time = '2dpe';
    elseif strfind(samples{i},'3dpe')
        T43_con{hamsterNum}.time = '3dpe';
    else
        error('Time of collection not properly specified');
    end

    % Discern tissue
    if strfind(samples{i},'NT')
        T43_con{hamsterNum}.NT.raw = T43_con{hamsterNum}.NT.raw + data_raw{:,i};
    elseif strfind(samples{i},'Trach')
        T43_con{hamsterNum}.Trach.raw = T43_con{hamsterNum}.Trach.raw + data_raw{:,i};
    elseif strfind(samples{i},'WL')
        T43_con{hamsterNum}.WL.raw = T43_con{hamsterNum}.WL.raw + data_raw{:,i};
    elseif strfind(samples{i},'LL')
        T43_con{hamsterNum}.LL.raw = T43_con{hamsterNum}.LL.raw + data_raw{:,i};
    elseif strfind(samples{i},'RL1')
        T43_con{hamsterNum}.RL1.raw = T43_con{hamsterNum}.RL1.raw + data_raw{:,i};
    elseif strfind(samples{i},'RL2')
        T43_con{hamsterNum}.RL2.raw = T43_con{hamsterNum}.RL2.raw + data_raw{:,i};
    elseif strfind(samples{i},'RL3')
        T43_con{hamsterNum}.RL3.raw = T43_con{hamsterNum}.RL3.raw + data_raw{:,i};
    elseif strfind(samples{i},'RL4')
        T43_con{hamsterNum}.RL4.raw = T43_con{hamsterNum}.RL4.raw + data_raw{:,i};
    else
        error('Tissue of collection not properly specified');
    end
    
end

for i = 16:length(T43_con)
    if sum(T43_con{i}.WL.raw) == 0
        T43_con{i}.WL.raw = T43_con{i}.LL.raw + T43_con{i}.RL1.raw + ...
            T43_con{i}.RL2.raw + T43_con{i}.RL3.raw + T43_con{i}.RL4.raw;
    end
end

fprintf('Raw data organized into structs\n')

%% Add info about which BCs are transmitted, how many BCs are transmitted
for i = 16:length(T43_con)
    T43_con{i}.NT.norm = T43_con{i}.NT.raw ./ sum(T43_con{i}.NT.raw);
    T43_con{i}.NT.trans = and((T43_con{i}.NT.raw>countCutoff), (T43_con{i}.NT.norm>MADCutoff));
    T43_con{i}.NT.raw = T43_con{i}.NT.raw .* T43_con{i}.NT.trans;
    T43_con{i}.NT.norm = T43_con{i}.NT.raw ./ sum(T43_con{i}.NT.raw);
    T43_con{i}.NT.nBCs = sum(T43_con{i}.NT.trans);

    T43_con{i}.Trach.norm = T43_con{i}.Trach.raw ./ sum(T43_con{i}.Trach.raw);
    T43_con{i}.Trach.trans = and((T43_con{i}.Trach.raw>countCutoff), (T43_con{i}.Trach.norm>MADCutoff));
    T43_con{i}.Trach.raw = T43_con{i}.Trach.raw .* T43_con{i}.Trach.trans;
    T43_con{i}.Trach.norm = T43_con{i}.Trach.raw ./ sum(T43_con{i}.Trach.raw);
    T43_con{i}.Trach.nBCs = sum(T43_con{i}.Trach.trans);

    T43_con{i}.WL.norm = T43_con{i}.WL.raw ./ sum(T43_con{i}.WL.raw);
    T43_con{i}.WL.trans = and((T43_con{i}.WL.raw>countCutoff), (T43_con{i}.WL.norm>MADCutoff));
    T43_con{i}.WL.raw = T43_con{i}.WL.raw .* T43_con{i}.WL.trans;
    T43_con{i}.WL.norm = T43_con{i}.WL.raw ./ sum(T43_con{i}.WL.raw);
    T43_con{i}.WL.nBCs = sum(T43_con{i}.WL.trans);

    T43_con{i}.totalBCs = sum(max([T43_con{i}.NT.trans, T43_con{i}.Trach.trans, T43_con{i}.WL.trans],[],2));
end

fprintf('Total BCs present determined\n')

%% Make stacked bar charts
% Set output folder
folder = 'Figure 2, S2 - T43 contacts/figs';
colmap = distinguishable_colors(nBCs);

makeFigs = input('Re-make stacked bar charts? ', 's');
labels = {'NT','Trachea','Lungs'};

if upper(makeFigs) == 'Y'
    for i = [17,23]%16:length(T43_con)
        if i < 26
            contactNum = i-10;
        else
            contactNum = i-25;
        end

        ylabelString = ['Contact #',int2str(contactNum)];

        if i == 17
            ylabelString = 'Contact, 24 hpe';
        else
            ylabelString = 'Contact, 48 hpe';
        end

        makeHorzBarChartVenn(T43_con{i}, labels, T43_con{i}.totalBCs, ylabelString, folder);
    end
end