clc; close all; %clear all

%% Load data
% Read in BC counts
file = 'Figure 5 - Onwards/Fig5_onwards_contacts.xlsx';
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
baseHamster.dur = '8hr'; % duration of exposure for contacts
baseHamster.contact = 'AB'; % route of exposure for contacts
baseHamster.Exp = 'T48';
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

T48 = {'','','','','','','','', ...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, ...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, ...
    baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, ...
    baseHamster};
% First 8 cells are spacers so that index = hamster number

samples = data_raw.Properties.VariableNames;

for i = 1:length(samples)
    % Determine hamster within experiment
    hamsterNumIndex = strfind(samples{i},'T48_H');
    hamsterNum = str2num(samples{i}(hamsterNumIndex+5:hamsterNumIndex+6));
    T48{hamsterNum}.animalNum = hamsterNum;

    % discern collection time
    %  note to self: make collection time read in more robust by
    %  checking if a value has already been assigned to it
    if strfind(samples{i},'Con')
        T48{hamsterNum}.contact = 'Con1';
    elseif strfind(samples{i},'C2')
        T48{hamsterNum}.contact = 'Con2';
    else
        error('Primary or secondary contact not properly specified');
    end
    

    % Discern tissue
    if strfind(samples{i},'NT')
        T48{hamsterNum}.NT.raw = T48{hamsterNum}.NT.raw + data_raw{:,i};
    elseif strfind(samples{i},'Trach')
        T48{hamsterNum}.Trach.raw = T48{hamsterNum}.Trach.raw + data_raw{:,i};
    elseif strfind(samples{i},'WL')
        T48{hamsterNum}.WL.raw = T48{hamsterNum}.WL.raw + data_raw{:,i};
    elseif strfind(samples{i},'LL')
        T48{hamsterNum}.LL.raw = T48{hamsterNum}.LL.raw + data_raw{:,i};
    elseif strfind(samples{i},'RL1')
        T48{hamsterNum}.RL1.raw = T48{hamsterNum}.RL1.raw + data_raw{:,i};
    elseif strfind(samples{i},'RL2')
        T48{hamsterNum}.RL2.raw = T48{hamsterNum}.RL2.raw + data_raw{:,i};
    elseif strfind(samples{i},'RL3')
        T48{hamsterNum}.RL3.raw = T48{hamsterNum}.RL3.raw + data_raw{:,i};
    elseif strfind(samples{i},'RL4')
        T48{hamsterNum}.RL4.raw = T48{hamsterNum}.RL4.raw + data_raw{:,i};
    else
        error('Tissue of collection not properly specified');
    end
end

% for i = 9:length(T48)
%     if sum(T48{i}.WL.raw) == 0
%         T48{i}.WL.raw = T48{i}.LL.raw + T48{i}.RL1.raw + ...
%             T48{i}.RL2.raw + T48{i}.RL3.raw + T48{i}.RL4.raw;
%     end
% end

fprintf('Raw data organized into structs\n')

%% Add info about which BCs are transmitted, how many BCs are transmitted
for i = 9:length(T48)
    T48{i}.NT.norm = T48{i}.NT.raw ./ sum(T48{i}.NT.raw);
    T48{i}.NT.trans = and((T48{i}.NT.raw>countCutoff), (T48{i}.NT.norm>MADCutoff));
    T48{i}.NT.raw = T48{i}.NT.raw .* T48{i}.NT.trans;
    T48{i}.NT.norm = T48{i}.NT.raw ./ sum(T48{i}.NT.raw);
    T48{i}.NT.nBCs = sum(T48{i}.NT.trans);

    T48{i}.Trach.norm = T48{i}.Trach.raw ./ sum(T48{i}.Trach.raw);
    T48{i}.Trach.trans = and((T48{i}.Trach.raw>countCutoff), (T48{i}.Trach.norm>MADCutoff));
    T48{i}.Trach.raw = T48{i}.Trach.raw .* T48{i}.Trach.trans;
    T48{i}.Trach.norm = T48{i}.Trach.raw ./ sum(T48{i}.Trach.raw);
    T48{i}.Trach.nBCs = sum(T48{i}.Trach.trans);

    T48{i}.WL.norm = T48{i}.WL.raw ./ sum(T48{i}.WL.raw);
    T48{i}.WL.trans = and((T48{i}.WL.raw>countCutoff), (T48{i}.WL.norm>MADCutoff));
    T48{i}.WL.raw = T48{i}.WL.raw .* T48{i}.WL.trans;
    T48{i}.WL.norm = T48{i}.WL.raw ./ sum(T48{i}.WL.raw);
    T48{i}.WL.nBCs = sum(T48{i}.WL.trans);

    T48{i}.totalTrans = max([T48{i}.NT.trans, T48{i}.Trach.trans, T48{i}.WL.trans],[],2);
    T48{i}.totalBCs = sum(T48{i}.totalTrans);

    T48{i}.NT.spec = and(T48{i}.NT.trans, ~max(T48{i}.Trach.trans, T48{i}.WL.trans));
    T48{i}.NT.nspec = sum(T48{i}.NT.spec);
    T48{i}.Trach.spec = and(T48{i}.Trach.trans, ~max(T48{i}.NT.trans, T48{i}.WL.trans));
    T48{i}.Trach.nspec = sum(T48{i}.Trach.spec);
    T48{i}.WL.spec = and(T48{i}.WL.trans, ~max(T48{i}.NT.trans, T48{i}.Trach.trans));
    T48{i}.WL.nspec = sum(T48{i}.WL.spec);
end

fprintf('Total BCs present determined\n')

%% Make stacked bar charts
% Set output folder
folder = 'Figure 5 - Onwards/figs';
colmap = distinguishable_colors(nBCs);

makeFigs = input('Re-make stacked bar charts? ', 's');
labels = {'NT','Trachea','Lungs'};

if upper(makeFigs) == 'Y'
    for i = [16 24]%[9 10 11 12 13 14 15 16 19 20 22 23 24]
        if i < 17
            ylabelstring = 'Primary Contact';
        else
            ylabelstring = 'Paired Secondary';
        end
        makeHorzBarChartSimVenn(T48{i}, labels, T48{i}.totalBCs, ylabelstring, folder);
    end
end

%% Compare C1s to C2s
% look for NT-specific BCs
NTfreqs = cell(1,24); allNTfreqs = []; 
TransNTfreqs = cell(1,24); allTransNTfreqs = [];
UntransNTfreqs = cell(1,24); allUntransNTfreqs = [];

Trachfreqs = cell(1,24); allTrachfreqs = []; 
TransTrachfreqs = cell(1,24); allTransTrachfreqs = [];
UntransTrachfreqs = cell(1,24); allUntransTrachfreqs = [];

WLfreqs = cell(1,24); allWLfreqs = [];
TransWLfreqs = cell(1,24); allTransWLfreqs = [];
UntransWLfreqs = cell(1,24); allUntransWLfreqs = [];
index = 1; allHamsterNums = [];
for i = 9:16
    % Identify all BCs present in C1 hamsters and their freqs in NT
    NTfreqs{i} = T48{i}.NT.norm(T48{i}.totalTrans);
    allNTfreqs = [allNTfreqs; NTfreqs{i}];

    % Identify all BCs present in C1 hamsters and their freqs in Trach
    Trachfreqs{i} = T48{i}.Trach.norm(T48{i}.totalTrans);
    allTrachfreqs = [allTrachfreqs; Trachfreqs{i}];

    % Identify all BCs present in C1 hamsters and their freqs in WL
    WLfreqs{i} = T48{i}.WL.norm(T48{i}.totalTrans);
    allWLfreqs = [allWLfreqs; WLfreqs{i}];

    allHamsterNums = [allHamsterNums, ones(1,T48{9}.totalBCs)*i];

    if any(i == [11 12 14 15 16])
        % Identify all BCs transmitted to C2 hamsters and their freqs in NT
        TransNTfreqs{i+8} = T48{i}.NT.norm(T48{i+8}.totalTrans);
        allTransNTfreqs = [allTransNTfreqs; TransNTfreqs{i+8}];
        % Identify all BCs not transmitted to C2 hamsters and their freqs in NT
        UntransNTfreqs{i+8} = T48{i}.NT.norm(and(T48{i}.totalTrans, ~T48{i+8}.totalTrans));
        allUntransNTfreqs = [allUntransNTfreqs; UntransNTfreqs{i+8}];
        % Identify all BCs transmitted to C2 hamsters and their freqs in Trach
        TransTrachfreqs{i+8} = T48{i}.Trach.norm(T48{i+8}.totalTrans);
        allTransTrachfreqs = [allTransTrachfreqs; TransTrachfreqs{i+8}];
        % Identify all BCs not transmitted to C2 hamsters and their freqs in Trach
        UntransTrachfreqs{i+8} = T48{i}.Trach.norm(and(T48{i}.totalTrans, ~T48{i+8}.totalTrans));
        allUntransTrachfreqs = [allUntransTrachfreqs; UntransTrachfreqs{i+8}];
        % Identify all BCs transmitted to C2 hamsters and their freqs in WL
        TransWLfreqs{i+8} = T48{i}.WL.norm(T48{i+8}.totalTrans);
        allTransWLfreqs = [allTransWLfreqs; TransWLfreqs{i+8}];
        % Identify all BCs not transmitted to C2 hamsters and their freqs in WL
        UntransWLfreqs{i+8} = T48{i}.WL.norm(and(T48{i}.totalTrans, ~T48{i+8}.totalTrans));
        allUntransWLfreqs = [allUntransWLfreqs; UntransWLfreqs{i+8}];
    end
end

nTrans = length(allTransNTfreqs);
nUntrans = length(allUntransNTfreqs);
resp = logical([ones(nTrans,1); zeros(nUntrans,1)]);

% similarity of C2 tissues to each C1 tissue
allSims = zeros(3,3,5);
index = 1;

for i = [19 20 22 23 24]
    simMatrix = CSCsim([T48{i-8}.NT.norm, T48{i-8}.Trach.norm, T48{i-8}.WL.norm, T48{i}.NT.norm, T48{i}.Trach.norm, T48{i}.WL.norm]);
    allSims(:,:,index) = simMatrix(1:3,4:6);
    index = index+1;
end

fprintf('Measured similarity to donor tissues\n');

%% Fit logistic model to each data set
mdlNT = fitglm([allTransNTfreqs; allUntransNTfreqs], ... % freqs in NT
    resp, ... % logical array indicating transmitted or not
    'Distribution', 'binomial', ...
    'Link', 'logit');
score_log_NT = mdlNT.Fitted.Probability;
[Xlog_NT,Ylog_NT,Tlog_NT,True_AUClog_NT] = perfcurve(resp,score_log_NT,'true');

mdlTrach = fitglm([allTransTrachfreqs; allUntransTrachfreqs], ... % freqs in Trach
    resp, ... % logical array indicating transmitted or not
    'Distribution', 'binomial', ...
    'Link', 'logit');
score_log_Trach = mdlTrach.Fitted.Probability;
[Xlog_Trach,Ylog_Trach,Tlog_Trach,True_AUClog_Trach] = perfcurve(resp,score_log_Trach,'true');

mdlWL = fitglm([allTransWLfreqs; allUntransWLfreqs], ... % freqs in WL
    resp, ... % logical array indicating transmitted or not
    'Distribution', 'binomial', ...
    'Link', 'logit');
score_log_WL = mdlWL.Fitted.Probability;
[Xlog_WL,Ylog_WL,Tlog_WL,True_AUClog_WL] = perfcurve(resp,score_log_WL,'true');


mdlAll_log = fitglm([allTransNTfreqs,allTransTrachfreqs,allTransWLfreqs; ...
    allUntransNTfreqs,allUntransTrachfreqs,allUntransWLfreqs], ... % freqs in NT
    resp, ... % logical array indicating transmitted or not
    'Distribution', 'binomial', ...
    'Link', 'logit');
score_log_All = mdlAll_log.Fitted.Probability;
[Xlog_All,Ylog_All,Tlog_All,True_AUClog_All] = perfcurve(resp,score_log_All,'true');

%% Simulate tissue-specific transmission events
nSim = 10000; percentDone = 0;
nBCTrans = zeros(nSim,1);
BCTrans = zeros(nSim,nTrans);
rand2BC = zeros(size(allNTfreqs));
for j = 2:length(allNTfreqs)
    rand2BC(j) = rand2BC(j-1) + allNTfreqs(j-1);
end
AUCs = zeros(nSim,3);

for i = 1:nSim
    hamsterNums = [];
    while nBCTrans(i) < nTrans
        randNum = 8*rand();
        hamsterNum = ceil(randNum)+8;
        BC = sum(rand2BC<randNum);

        if ~any(BCTrans(i,:) == BC)
            nBCTrans(i) = nBCTrans(i)+1;
            BCTrans(i,nBCTrans(i)) = BC;

            if ~any(hamsterNum == hamsterNums)
                hamsterNums = [hamsterNums,hamsterNum];
            end
        end
    end

    BCsubset = ismember(allHamsterNums, hamsterNums);

    % NOTE: Do not account for hamsters that did not transmit
    resp = zeros(size(allNTfreqs));
    resp(BCTrans(i,:)) = ones(1,nTrans);
    resp = logical(resp);

    mdlNT = fitglm(allNTfreqs(BCsubset), ... % freqs in NT
        resp(BCsubset), ... % logical array indicating transmitted or not
        'Distribution', 'binomial', ...
        'Link', 'logit');
    score_log_NT = mdlNT.Fitted.Probability;
    [Xlog_NT,Ylog_NT,Tlog_NT,AUClog_NT] = perfcurve(resp(BCsubset),score_log_NT,'true');
    
    mdlTrach = fitglm(allTrachfreqs(BCsubset), ... % freqs in Trach
        resp(BCsubset), ... % logical array indicating transmitted or not
        'Distribution', 'binomial', ...
        'Link', 'logit');
    score_log_Trach = mdlTrach.Fitted.Probability;
    [Xlog_Trach,Ylog_Trach,Tlog_Trach,AUClog_Trach] = perfcurve(resp(BCsubset),score_log_Trach,'true');
    
    mdlWL = fitglm(allWLfreqs(BCsubset), ... % freqs in WL
        resp(BCsubset), ... % logical array indicating transmitted or not
        'Distribution', 'binomial', ...
        'Link', 'logit');
    score_log_WL = mdlWL.Fitted.Probability;
    [Xlog_WL,Ylog_WL,Tlog_WL,AUClog_WL] = perfcurve(resp(BCsubset),score_log_WL,'true');
    
    AUCs(i,:) = [AUClog_NT, AUClog_Trach, AUClog_WL];

    if floor(100*i / nSim) > percentDone
        percentDone = percentDone+10;
        fprintf('%d%% done\n', percentDone);
    end
end

nBins = 100;
histEdges = linspace(0, 1, nBins+1);
cols = rainbowCols(3);

f = figure();
hold off;
histogram(AUCs(:,1),histEdges, 'FaceColor',cols(1,:), 'FaceAlpha',0.4);
hold on;
histogram(AUCs(:,2),histEdges, 'FaceColor',cols(2,:), 'FaceAlpha',0.4);
histogram(AUCs(:,3),histEdges, 'FaceColor',cols(3,:), 'FaceAlpha',0.4);
xlim([0 1]);
legend('NT','Trachea','Lungs', 'Location','best');

NTarrow = annotation('textarrow', [True_AUClog_NT,True_AUClog_NT], [0.1,0],'String','NT');
NTarrow.Parent = f.CurrentAxes;
NTarrow.X = [True_AUClog_NT,True_AUClog_NT]; NTarrow.Y = [nSim*0.025,0];

Tracharrow = annotation('textarrow', [True_AUClog_Trach,True_AUClog_Trach], [0.1,0],'String','Trachea');
Tracharrow.Parent = f.CurrentAxes;
Tracharrow.X = [True_AUClog_Trach,True_AUClog_Trach]; Tracharrow.Y = [nSim*0.025,0];

WLarrow = annotation('textarrow', [True_AUClog_WL,True_AUClog_WL], [0.1,0],'String','Lungs');
WLarrow.Parent = f.CurrentAxes;
WLarrow.X = [True_AUClog_WL,True_AUClog_WL]; WLarrow.Y = [nSim*0.025,0];


%% Functions
