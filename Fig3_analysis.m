clc; close all; %clear all;

%% Load data
% Read in BC counts
file = 'Figure 3 - time of exposure/Fig3_Time_Contacts.xlsx';
data_raw = readtable(file, 'ReadRowNames', true, 'PreserveVariableNames', true);
fprintf('Data read in\n')
data_norm = data_raw ./ sum(data_raw);

MADCutoff = 0.00056;
countCutoff = 10;

nBCs = size(data_raw,1);

% data_norm(or(data_norm<MADCutoff, data_raw<countCutoff)) = 0;
% data_norm = data_norm ./ sum(data_norm);

%% Organize data into structs
baseHamster.DonCon = 'Con'; % donor or contact
baseHamster.sex = 'M'; % male or female
% hamster.dur = '8hr'; % duration of exposure for contacts
hamster.route = 'AB'; % route of exposure for contacts
baseHamster.time = '3dpe'; % For donors, hpi of collection. For contacts, dpe of harvest
% baseHamster.Exp = 'T39';

baseHamster.NT.raw = zeros(nBCs,1); % raw reads for each BC in nasal turbinate
baseHamster.Trach.raw = zeros(nBCs,1); % raw reads for each BC in trachea
baseHamster.WL.raw = zeros(nBCs,1); % raw reads for each BC in lungs
baseHamster.LL.raw = zeros(nBCs,1); % raw reads for each BC in left lung lobe
baseHamster.RL1.raw = zeros(nBCs,1); % raw reads for each BC in right apical lung lobe
baseHamster.RL2.raw = zeros(nBCs,1); % raw reads for each BC in right middle lung lobe
baseHamster.RL3.raw = zeros(nBCs,1); % raw reads for each BC in right caudal lung lobe
baseHamster.RL4.raw = zeros(nBCs,1); % raw reads for each BC in infracardial lung lobe

T39 = {'','','','','','','','', ...
    baseHamster, baseHamster, ...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, ...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster};
% First 8 cells are spacers so that index = hamster number
T46 = {'','','','','','','','', ...
    baseHamster, baseHamster, ...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, ...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster};
% First 8 cells are spacers so that index = hamster number

samples = data_raw.Properties.VariableNames;

for i = 1:length(samples)
    if strfind(samples{i},'T39')
        % Determine hamster within experiment
        hamsterNumIndex = strfind(samples{i},'T39_H');
        hamsterNum = str2num(samples{i}(hamsterNumIndex+5:hamsterNumIndex+6));
        T39{hamsterNum}.animalNum = hamsterNum;
        T39{hamsterNum}.Exp = 'T39';
    
        % discern collection time
        %  note to self: make collection time read in more robust by
        %  checking if a value has already been assigned to it
        if strfind(samples{i},'1hr')
            T39{hamsterNum}.dur = '1hr';
        elseif strfind(samples{i},'4hr')
            T39{hamsterNum}.DonCon = '4hr';
        elseif strfind(samples{i},'8hr')
            T39{hamsterNum}.DonCon = '8hr';
        else
            error('Time of exposure not properly specified');
        end
        
    
        % Discern tissue
        if strfind(samples{i},'NT')
            T39{hamsterNum}.NT.raw = T39{hamsterNum}.NT.raw + data_raw{:,i};
        elseif strfind(samples{i},'Trach')
            T39{hamsterNum}.Trach.raw = T39{hamsterNum}.Trach.raw + data_raw{:,i};
        elseif strfind(samples{i},'WL')
            T39{hamsterNum}.WL.raw = T39{hamsterNum}.WL.raw + data_raw{:,i};
        elseif strfind(samples{i},'LL')
            T39{hamsterNum}.LL.raw = T39{hamsterNum}.LL.raw + data_raw{:,i};
        elseif strfind(samples{i},'RL1')
            T39{hamsterNum}.RL1.raw = T39{hamsterNum}.RL1.raw + data_raw{:,i};
        elseif strfind(samples{i},'RL2')
            T39{hamsterNum}.RL2.raw = T39{hamsterNum}.RL2.raw + data_raw{:,i};
        elseif strfind(samples{i},'RL3')
            T39{hamsterNum}.RL3.raw = T39{hamsterNum}.RL3.raw + data_raw{:,i};
        elseif strfind(samples{i},'RL4')
            T39{hamsterNum}.RL4.raw = T39{hamsterNum}.RL4.raw + data_raw{:,i};
        else
            error('Tissue of collection not properly specified');
        end

    elseif strfind(samples{i},'T46')
        % Determine hamster within experiment
        hamsterNumIndex = strfind(samples{i},'T46_H');
        hamsterNum = str2num(samples{i}(hamsterNumIndex+5:hamsterNumIndex+6));
        T46{hamsterNum}.animalNum = hamsterNum;
        T46{hamsterNum}.Exp = 'T46';
    
        % discern collection time
        %  note to self: make collection time read in more robust by
        %  checking if a value has already been assigned to it
        if strfind(samples{i},'1hr')
            T46{hamsterNum}.dur = '1hr';
        elseif strfind(samples{i},'4hr')
            T46{hamsterNum}.DonCon = '4hr';
        elseif strfind(samples{i},'8hr')
            T46{hamsterNum}.DonCon = '8hr';
        else
            error('Time of exposure not properly specified');
        end
        
    
        % Discern tissue
        if strfind(samples{i},'NT')
            T46{hamsterNum}.NT.raw = T46{hamsterNum}.NT.raw + data_raw{:,i};
        elseif strfind(samples{i},'Trach')
            T46{hamsterNum}.Trach.raw = T46{hamsterNum}.Trach.raw + data_raw{:,i};
        elseif strfind(samples{i},'WL')
            T46{hamsterNum}.WL.raw = T46{hamsterNum}.WL.raw + data_raw{:,i};
        elseif strfind(samples{i},'LL')
            T46{hamsterNum}.LL.raw = T46{hamsterNum}.LL.raw + data_raw{:,i};
        elseif strfind(samples{i},'RL1')
            T46{hamsterNum}.RL1.raw = T46{hamsterNum}.RL1.raw + data_raw{:,i};
        elseif strfind(samples{i},'RL2')
            T46{hamsterNum}.RL2.raw = T46{hamsterNum}.RL2.raw + data_raw{:,i};
        elseif strfind(samples{i},'RL3')
            T46{hamsterNum}.RL3.raw = T46{hamsterNum}.RL3.raw + data_raw{:,i};
        elseif strfind(samples{i},'RL4')
            T46{hamsterNum}.RL4.raw = T46{hamsterNum}.RL4.raw + data_raw{:,i};
        else
            error('Tissue of collection not properly specified');
        end
    end
end

% for i = 9:length(T39)
%     if sum(T39{i}.WL.raw) == 0
%         T39{i}.WL.raw = T39{i}.LL.raw + T39{i}.RL1.raw + ...
%             T39{i}.RL2.raw + T39{i}.RL3.raw + T39{i}.RL4.raw;
%     end
% end
% 
% for i = 9:length(T46)
%     if sum(T46{i}.WL.raw) == 0
%         T46{i}.WL.raw = T46{i}.LL.raw + T46{i}.RL1.raw + ...
%             T46{i}.RL2.raw + T46{i}.RL3.raw + T46{i}.RL4.raw;
%     end
% end

fprintf('Raw data organized into structs\n')

%% Add info about which BCs are transmitted, how many BCs are transmitted
for i = 9:length(T39)
    T39{i}.NT.norm = T39{i}.NT.raw ./ sum(T39{i}.NT.raw);
    T39{i}.NT.trans = and((T39{i}.NT.raw>countCutoff), (T39{i}.NT.norm>MADCutoff));
    T39{i}.NT.raw = T39{i}.NT.raw .* T39{i}.NT.trans;
    T39{i}.NT.norm = T39{i}.NT.raw ./ sum(T39{i}.NT.raw);
    T39{i}.NT.nBCs = sum(T39{i}.NT.trans);

    T39{i}.Trach.norm = T39{i}.Trach.raw ./ sum(T39{i}.Trach.raw);
    T39{i}.Trach.trans = and((T39{i}.Trach.raw>countCutoff), (T39{i}.Trach.norm>MADCutoff));
    T39{i}.Trach.raw = T39{i}.Trach.raw .* T39{i}.Trach.trans;
    T39{i}.Trach.norm = T39{i}.Trach.raw ./ sum(T39{i}.Trach.raw);
    T39{i}.Trach.nBCs = sum(T39{i}.Trach.trans);

    T39{i}.WL.norm = T39{i}.WL.raw ./ sum(T39{i}.WL.raw);
    T39{i}.WL.trans = and((T39{i}.WL.raw>countCutoff), (T39{i}.WL.norm>MADCutoff));
    T39{i}.WL.raw = T39{i}.WL.raw .* T39{i}.WL.trans;
    T39{i}.WL.norm = T39{i}.WL.raw ./ sum(T39{i}.WL.raw);
    T39{i}.WL.nBCs = sum(T39{i}.WL.trans);

    T39{i}.totalBCs = sum(max([T39{i}.NT.trans, T39{i}.Trach.trans, T39{i}.WL.trans],[],2));
end

for i = 9:length(T46)
    T46{i}.NT.norm = T46{i}.NT.raw ./ sum(T46{i}.NT.raw);
    T46{i}.NT.trans = and((T46{i}.NT.raw>countCutoff), (T46{i}.NT.norm>MADCutoff));
    T46{i}.NT.raw = T46{i}.NT.raw .* T46{i}.NT.trans;
    T46{i}.NT.norm = T46{i}.NT.raw ./ sum(T46{i}.NT.raw);
    T46{i}.NT.nBCs = sum(T46{i}.NT.trans);

    T46{i}.Trach.norm = T46{i}.Trach.raw ./ sum(T46{i}.Trach.raw);
    T46{i}.Trach.trans = and((T46{i}.Trach.raw>countCutoff), (T46{i}.Trach.norm>MADCutoff));
    T46{i}.Trach.raw = T46{i}.Trach.raw .* T46{i}.Trach.trans;
    T46{i}.Trach.norm = T46{i}.Trach.raw ./ sum(T46{i}.Trach.raw);
    T46{i}.Trach.nBCs = sum(T46{i}.Trach.trans);

    T46{i}.WL.norm = T46{i}.WL.raw ./ sum(T46{i}.WL.raw);
    T46{i}.WL.trans = and((T46{i}.WL.raw>countCutoff), (T46{i}.WL.norm>MADCutoff));
    T46{i}.WL.raw = T46{i}.WL.raw .* T46{i}.WL.trans;
    T46{i}.WL.norm = T46{i}.WL.raw ./ sum(T46{i}.WL.raw);
    T46{i}.WL.nBCs = sum(T46{i}.WL.trans);

    T46{i}.totalBCs = sum(max([T46{i}.NT.trans, T46{i}.Trach.trans, T46{i}.WL.trans],[],2));
end

fprintf('Total BCs present determined\n')

%% Make stacked bar charts
% Set output folder
folder = 'Figure 3 - time of exposure/figs';
colmap = distinguishable_colors(nBCs);

makeFigs = input('Re-make stacked bar charts? ', 's');
labels = {'NT','Trachea','Lungs'};

if upper(makeFigs) == 'Y'
    for i = [9,10,11,12,13,14,15,16,18,19]
        makeHorzBarChartSimVenn(T39{i}, labels, T39{i}.totalBCs, folder);
    end
end

if upper(makeFigs) == 'Y'
    for i = [9,10,11,12,13,14,15,16,18,20]
        makeHorzBarChartSimVenn(T46{i}, labels, T46{i}.totalBCs, folder);
    end
end