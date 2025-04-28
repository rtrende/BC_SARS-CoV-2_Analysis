clc; close all; %clear all; 

%% Load data
% Read in BC counts
file = 'Figure 4 - route of exposure/Fig4_Route_Contacts.xlsx';
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
hamster.dur = '8hr'; % duration of exposure for contacts
hamster.route = ''; % route of exposure for contacts
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

T41 = {'','','','','','','','','','', ...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, ...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster};
% First 10 cells are spacers so that index = hamster number

T46 = {'','','','','','','','','','','','','','','','','','','','', ...
    baseHamster, baseHamster, baseHamster, baseHamster};
% first 20 cells are spacers so that index = hamster number

T47 = {'','','','','','','','','','', ...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, ...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, ...
    baseHamster, baseHamster, baseHamster, baseHamster};
% First 10 cells are spacers so that index = hamster number

samples = data_raw.Properties.VariableNames;

for i = 1:length(samples)
    if strfind(samples{i},'T41')
        % Determine hamster within experiment
        hamsterNumIndex = strfind(samples{i},'T41_H');
        hamsterNum = str2num(samples{i}(hamsterNumIndex+5:hamsterNumIndex+6));
        T41{hamsterNum}.animalNum = hamsterNum;
        T41{hamsterNum}.Exp = 'T41';
    
        % discern collection time
        %  note to self: make collection time read in more robust by
        %  checking if a value has already been assigned to it
        if strfind(samples{i},'AB')
            T41{hamsterNum}.route = 'AB';
        elseif strfind(samples{i},'DC')
            T41{hamsterNum}.route = 'DC';
        elseif strfind(samples{i},'Fom')
            T41{hamsterNum}.route = 'Fom';
        else
            error('Route of exposure not properly specified');
        end
        
    
        % Discern tissue
        if strfind(samples{i},'NT')
            T41{hamsterNum}.NT.raw = T41{hamsterNum}.NT.raw + data_raw{:,i};
        elseif strfind(samples{i},'Trach')
            T41{hamsterNum}.Trach.raw = T41{hamsterNum}.Trach.raw + data_raw{:,i};
        elseif strfind(samples{i},'WL')
            T41{hamsterNum}.WL.raw = T41{hamsterNum}.WL.raw + data_raw{:,i};
        elseif strfind(samples{i},'LL')
            T41{hamsterNum}.LL.raw = T41{hamsterNum}.LL.raw + data_raw{:,i};
        elseif strfind(samples{i},'RL1')
            T41{hamsterNum}.RL1.raw = T41{hamsterNum}.RL1.raw + data_raw{:,i};
        elseif strfind(samples{i},'RL2')
            T41{hamsterNum}.RL2.raw = T41{hamsterNum}.RL2.raw + data_raw{:,i};
        elseif strfind(samples{i},'RL3')
            T41{hamsterNum}.RL3.raw = T41{hamsterNum}.RL3.raw + data_raw{:,i};
        elseif strfind(samples{i},'RL4')
            T41{hamsterNum}.RL4.raw = T41{hamsterNum}.RL4.raw + data_raw{:,i};
        else
            error('Tissue of collection not properly specified');
        end

    elseif strfind(samples{i},'T47')
        % Determine hamster within experiment
        hamsterNumIndex = strfind(samples{i},'T47_H');
        hamsterNum = str2num(samples{i}(hamsterNumIndex+5:hamsterNumIndex+6));
        T47{hamsterNum}.animalNum = hamsterNum;
        T47{hamsterNum}.Exp = 'T47';
    
        % discern collection time
        %  note to self: make collection time read in more robust by
        %  checking if a value has already been assigned to it
        if strfind(samples{i},'AB')
            T47{hamsterNum}.route = 'AB';
        elseif strfind(samples{i},'DC')
            T47{hamsterNum}.route = 'DC';
        elseif strfind(samples{i},'Fom')
            T47{hamsterNum}.route = 'Fom';
        else
            error('Route of exposure not properly specified');
        end
        
    
        % Discern tissue
        if strfind(samples{i},'NT')
            T47{hamsterNum}.NT.raw = T47{hamsterNum}.NT.raw + data_raw{:,i};
        elseif strfind(samples{i},'Trach')
            T47{hamsterNum}.Trach.raw = T47{hamsterNum}.Trach.raw + data_raw{:,i};
        elseif strfind(samples{i},'WL')
            T47{hamsterNum}.WL.raw = T47{hamsterNum}.WL.raw + data_raw{:,i};
        elseif strfind(samples{i},'LL')
            T47{hamsterNum}.LL.raw = T47{hamsterNum}.LL.raw + data_raw{:,i};
        elseif strfind(samples{i},'RL1')
            T47{hamsterNum}.RL1.raw = T47{hamsterNum}.RL1.raw + data_raw{:,i};
        elseif strfind(samples{i},'RL2')
            T47{hamsterNum}.RL2.raw = T47{hamsterNum}.RL2.raw + data_raw{:,i};
        elseif strfind(samples{i},'RL3')
            T47{hamsterNum}.RL3.raw = T47{hamsterNum}.RL3.raw + data_raw{:,i};
        elseif strfind(samples{i},'RL4')
            T47{hamsterNum}.RL4.raw = T47{hamsterNum}.RL4.raw + data_raw{:,i};
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
        if strfind(samples{i},'AB')
            T46{hamsterNum}.route = 'AB';
        elseif strfind(samples{i},'DC')
            T46{hamsterNum}.route = 'DC';
        elseif strfind(samples{i},'Fom')
            T46{hamsterNum}.route = 'Fom';
        else
            error('Route of exposure not properly specified');
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

% for i = 11:length(T41)
%     if sum(T41{i}.WL.raw) == 0
%         T41{i}.WL.raw = T41{i}.LL.raw + T41{i}.RL1.raw + ...
%             T41{i}.RL2.raw + T41{i}.RL3.raw + T41{i}.RL4.raw;
%     end
% end
% 
% for i = 11:length(T47)
%     if sum(T47{i}.WL.raw) == 0
%         T47{i}.WL.raw = T47{i}.LL.raw + T47{i}.RL1.raw + ...
%             T47{i}.RL2.raw + T47{i}.RL3.raw + T47{i}.RL4.raw;
%     end
% end
% 
% for i = 21:length(T46)
%     if sum(T46{i}.WL.raw) == 0
%         T46{i}.WL.raw = T46{i}.LL.raw + T46{i}.RL1.raw + ...
%             T46{i}.RL2.raw + T46{i}.RL3.raw + T46{i}.RL4.raw;
%     end
% end

fprintf('Raw data organized into structs\n')

%% Add info about which BCs are transmitted, how many BCs are transmitted
for i = 11:length(T41)
    T41{i}.NT.norm = T41{i}.NT.raw ./ sum(T41{i}.NT.raw);
    T41{i}.NT.trans = and((T41{i}.NT.raw>countCutoff), (T41{i}.NT.norm>MADCutoff));
    T41{i}.NT.raw = T41{i}.NT.raw .* T41{i}.NT.trans;
    T41{i}.NT.norm = T41{i}.NT.raw ./ sum(T41{i}.NT.raw);
    T41{i}.NT.nBCs = sum(T41{i}.NT.trans);

    T41{i}.Trach.norm = T41{i}.Trach.raw ./ sum(T41{i}.Trach.raw);
    T41{i}.Trach.trans = and((T41{i}.Trach.raw>countCutoff), (T41{i}.Trach.norm>MADCutoff));
    T41{i}.Trach.raw = T41{i}.Trach.raw .* T41{i}.Trach.trans;
    T41{i}.Trach.norm = T41{i}.Trach.raw ./ sum(T41{i}.Trach.raw);
    T41{i}.Trach.nBCs = sum(T41{i}.Trach.trans);

    T41{i}.WL.norm = T41{i}.WL.raw ./ sum(T41{i}.WL.raw);
    T41{i}.WL.trans = and((T41{i}.WL.raw>countCutoff), (T41{i}.WL.norm>MADCutoff));
    T41{i}.WL.raw = T41{i}.WL.raw .* T41{i}.WL.trans;
    T41{i}.WL.norm = T41{i}.WL.raw ./ sum(T41{i}.WL.raw);
    T41{i}.WL.nBCs = sum(T41{i}.WL.trans);

    T41{i}.totalBCs = sum(max([T41{i}.NT.trans, T41{i}.Trach.trans, T41{i}.WL.trans],[],2));
end

for i = 11:length(T47)
    T47{i}.NT.norm = T47{i}.NT.raw ./ sum(T47{i}.NT.raw);
    T47{i}.NT.trans = and((T47{i}.NT.raw>countCutoff), (T47{i}.NT.norm>MADCutoff));
    T47{i}.NT.raw = T47{i}.NT.raw .* T47{i}.NT.trans;
    T47{i}.NT.norm = T47{i}.NT.raw ./ sum(T47{i}.NT.raw);
    T47{i}.NT.nBCs = sum(T47{i}.NT.trans);

    T47{i}.Trach.norm = T47{i}.Trach.raw ./ sum(T47{i}.Trach.raw);
    T47{i}.Trach.trans = and((T47{i}.Trach.raw>countCutoff), (T47{i}.Trach.norm>MADCutoff));
    T47{i}.Trach.raw = T47{i}.Trach.raw .* T47{i}.Trach.trans;
    T47{i}.Trach.norm = T47{i}.Trach.raw ./ sum(T47{i}.Trach.raw);
    T47{i}.Trach.nBCs = sum(T47{i}.Trach.trans);

    T47{i}.WL.norm = T47{i}.WL.raw ./ sum(T47{i}.WL.raw);
    T47{i}.WL.trans = and((T47{i}.WL.raw>countCutoff), (T47{i}.WL.norm>MADCutoff));
    T47{i}.WL.raw = T47{i}.WL.raw .* T47{i}.WL.trans;
    T47{i}.WL.norm = T47{i}.WL.raw ./ sum(T47{i}.WL.raw);
    T47{i}.WL.nBCs = sum(T47{i}.WL.trans);

    T47{i}.totalBCs = sum(max([T47{i}.NT.trans, T47{i}.Trach.trans, T47{i}.WL.trans],[],2));
end

for i = 21:length(T46)
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
folder = 'Figure 4 - route of exposure/figs';
colmap = distinguishable_colors(nBCs);

makeFigs = input('Re-make stacked bar charts? ', 's');
labels = {'NT','Trachea','Lungs'};

if upper(makeFigs) == 'Y'
    for i = [14 19]%11:length(T41)
        if i < 16
            ylabelstring = 'Airborne';
        elseif i < 21
            ylabelstring = 'Direct Contact';
        else
            ylabelstring = 'Fomite';
        end
        makeHorzBarChartVenn(T41{i}, labels, T41{i}.totalBCs, ylabelstring, folder);
    end

    % for i = 21:length(T46)
    %     makeHorzBarChartVenn(T46{i}, labels, T46{i}.totalBCs, 'Fomites', folder);
    % end

    for i = 23% 11:length(T47)
        if i < 16
            ylabelstring = 'Airborne';
        elseif i < 21
            ylabelstring = 'Direct Contact';
        else
            ylabelstring = 'Fomite';
        end
        makeHorzBarChartVenn(T47{i}, labels, T47{i}.totalBCs, ylabelstring, folder);
    end
end
% makeHorzBarChartVenn(data, labels_cell, totalBCs, ylabelString, outputFolder)