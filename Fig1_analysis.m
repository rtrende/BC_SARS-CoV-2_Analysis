clc; close all; % clear all;

%% Load data
% Read in BC counts
file = 'Figure 1 - Representative BC Dists/Fig1_T41-T47_Donors.xlsx';
file2 = 'Figure 1 - Representative BC Dists/Fig1_Inoculum_BC_seqs.xlsx';
data_raw = readtable(file, 'ReadRowNames', true, 'PreserveVariableNames', true);
inoc_data_raw = readtable(file2, 'ReadRowNames', true, 'PreserveVariableNames', true);
fprintf('Data read in\n')
data_norm = data_raw ./ sum(data_raw);
inoc_data_norm = inoc_data_raw ./ sum(inoc_data_raw);

MADCutoff = 0.00056;
countCutoff = 10;

nBCs = size(data_raw,1);

% data_norm(or(data_norm<MADCutoff, data_raw<countCutoff)) = 0;
% data_norm = data_norm ./ sum(data_norm);

%% Organize data into structs
baseHamster.DonCon = 'Donor'; % donor or contact
baseHamster.sex = 'M'; % male or female
% hamster.dur = 'NA'; % duration of exposure for contacts
% hamster.route = 'NA'; % route of exposure for contacts
baseHamster.time = ''; % For donors, hpi of collection

baseHamster.NT.raw = zeros(nBCs,1); % raw reads for each BC in nasal turbinate
baseHamster.Trach.raw = zeros(nBCs,1); % raw reads for each BC in trachea
baseHamster.WL.raw = zeros(nBCs,1); % raw reads for each BC in lungs

T41_donors = {baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, baseHamster};
T47_donors = {baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, baseHamster, baseHamster};

samples = data_raw.Properties.VariableNames;

for i = 1:length(samples)
    if strfind(samples{i}, 'T41') % Determine experiment
        % Determine hamster within experimen
        hamsterNumIndex = strfind(samples{i},'T41_H');
        hamsterNum = str2num(samples{i}(hamsterNumIndex+5:hamsterNumIndex+6));
        T41_donors{hamsterNum}.animalNum = hamsterNum;
        T41_donors{hamsterNum}.Exp = 'T41';

        % discern collection time
        %  note to self: make collection time read in more robust by
        %  checking if a value has already been assigned to it
        if strfind(samples{i},'32h')
            T41_donors{hamsterNum}.time = '32h';
        elseif strfind(samples{i},'72h')
            T41_donors{hamsterNum}.time = '72h';
        else
            error('Time of collection not properly specified');
        end

        % Discern tissue
        if strfind(samples{i},'NT')
            T41_donors{hamsterNum}.NT.raw = T41_donors{hamsterNum}.NT.raw + data_raw{:,i};
        elseif strfind(samples{i},'Trach')
            T41_donors{hamsterNum}.Trach.raw = T41_donors{hamsterNum}.Trach.raw + data_raw{:,i};
        elseif strfind(samples{i},'WL')
            T41_donors{hamsterNum}.WL.raw = T41_donors{hamsterNum}.WL.raw + data_raw{:,i};
        else
            error('Tissue of collection not properly specified');
        end
        
    elseif strfind(samples{i}, 'T47') % Determine experiment
        % Determine hamster within experiment
        hamsterNumIndex = strfind(samples{i},'T47_H');
        hamsterNum = str2num(samples{i}(hamsterNumIndex+5:hamsterNumIndex+6));
        T47_donors{hamsterNum}.animalNum = hamsterNum;
        T47_donors{hamsterNum}.Exp = 'T47';

        % discern collection time
        %  note to self: make collection time read in more robust by
        %  checking if a value has already been assigned to it
        if strfind(samples{i},'32h')
            T47_donors{hamsterNum}.time = '32h';
        elseif strfind(samples{i},'72h')
            T47_donors{hamsterNum}.time = '72h';
        else
            error('Time of collection not properly specified');
        end

        % Discern tissue
        if strfind(samples{i},'NT')
            T47_donors{hamsterNum}.NT.raw = T47_donors{hamsterNum}.NT.raw + data_raw{:,i};
        elseif strfind(samples{i},'Trach')
            T47_donors{hamsterNum}.Trach.raw = T47_donors{hamsterNum}.Trach.raw + data_raw{:,i};
        elseif strfind(samples{i},'WL')
            T47_donors{hamsterNum}.WL.raw = T47_donors{hamsterNum}.WL.raw + data_raw{:,i};
        else
            error('Tissue of collection not properly specified');
        end

    else
        error('Experiment number not properly specified')
    end
end

fprintf('Raw data organized into structs\n')

%% Add info about which BCs are transmitted, how many BCs are transmitted
for i = 1:length(T41_donors)
    T41_donors{i}.NT.norm = T41_donors{i}.NT.raw ./ sum(T41_donors{i}.NT.raw);
    T41_donors{i}.NT.trans = and((T41_donors{i}.NT.raw>countCutoff), (T41_donors{i}.NT.norm>MADCutoff));
    T41_donors{i}.NT.raw = T41_donors{i}.NT.raw .* T41_donors{i}.NT.trans;
    T41_donors{i}.NT.norm = T41_donors{i}.NT.raw ./ sum(T41_donors{i}.NT.raw);
    T41_donors{i}.NT.nBCs = sum(T41_donors{i}.NT.trans);

    T41_donors{i}.Trach.norm = T41_donors{i}.Trach.raw ./ sum(T41_donors{i}.Trach.raw);
    T41_donors{i}.Trach.trans = and((T41_donors{i}.Trach.raw>countCutoff), (T41_donors{i}.Trach.norm>MADCutoff));
    T41_donors{i}.Trach.raw = T41_donors{i}.Trach.raw .* T41_donors{i}.Trach.trans;
    T41_donors{i}.Trach.norm = T41_donors{i}.Trach.raw ./ sum(T41_donors{i}.Trach.raw);
    T41_donors{i}.Trach.nBCs = sum(T41_donors{i}.Trach.trans);

    T41_donors{i}.WL.norm = T41_donors{i}.WL.raw ./ sum(T41_donors{i}.WL.raw);
    T41_donors{i}.WL.trans = and((T41_donors{i}.WL.raw>countCutoff), (T41_donors{i}.WL.norm>MADCutoff));
    T41_donors{i}.WL.raw = T41_donors{i}.WL.raw .* T41_donors{i}.WL.trans;
    T41_donors{i}.WL.norm = T41_donors{i}.WL.raw ./ sum(T41_donors{i}.WL.raw);
    T41_donors{i}.WL.nBCs = sum(T41_donors{i}.WL.trans);

    T41_donors{i}.totalBCs = sum(max([T41_donors{i}.NT.trans, T41_donors{i}.Trach.trans, T41_donors{i}.WL.trans],[],2));
end

for i = 1:length(T47_donors)
    T47_donors{i}.NT.norm = T47_donors{i}.NT.raw ./ sum(T47_donors{i}.NT.raw);
    T47_donors{i}.NT.trans = and((T47_donors{i}.NT.raw>countCutoff), (T47_donors{i}.NT.norm>MADCutoff));
    T47_donors{i}.NT.raw = T47_donors{i}.NT.raw .* T47_donors{i}.NT.trans;
    T47_donors{i}.NT.norm = T47_donors{i}.NT.raw ./ sum(T47_donors{i}.NT.raw);
    T47_donors{i}.NT.nBCs = sum(T47_donors{i}.NT.trans);

    T47_donors{i}.Trach.norm = T47_donors{i}.Trach.raw ./ sum(T47_donors{i}.Trach.raw);
    T47_donors{i}.Trach.trans = and((T47_donors{i}.Trach.raw>countCutoff), (T47_donors{i}.Trach.norm>MADCutoff));
    T47_donors{i}.Trach.raw = T47_donors{i}.Trach.raw .* T47_donors{i}.Trach.trans;
    T47_donors{i}.Trach.norm = T47_donors{i}.Trach.raw ./ sum(T47_donors{i}.Trach.raw);
    T47_donors{i}.Trach.nBCs = sum(T47_donors{i}.Trach.trans);

    T47_donors{i}.WL.norm = T47_donors{i}.WL.raw ./ sum(T47_donors{i}.WL.raw);
    T47_donors{i}.WL.trans = and((T47_donors{i}.WL.raw>countCutoff), (T47_donors{i}.WL.norm>MADCutoff));
    T47_donors{i}.WL.raw = T47_donors{i}.WL.raw .* T47_donors{i}.WL.trans;
    T47_donors{i}.WL.norm = T47_donors{i}.WL.raw ./ sum(T47_donors{i}.WL.raw);
    T47_donors{i}.WL.nBCs = sum(T47_donors{i}.WL.trans);

    T47_donors{i}.totalBCs = sum(max([T47_donors{i}.NT.trans, T47_donors{i}.Trach.trans, T47_donors{i}.WL.trans],[],2));
end

fprintf('Total BCs present determined\n')

%% Make stacked bar charts
% Set output folder
folder = 'Figure 1 - Representative BC Dists/figs';
colmap = distinguishable_colors(nBCs);

makeFigs = input('Re-make stacked bar charts? ', 's');
labels = {'NT','Trachea','Lungs'};

if upper(makeFigs) == 'Y'
    for i = 1:length(T41_donors)
        makeHorzBarChartSimVenn(T41_donors{i}, labels, T41_donors{i}.totalBCs, folder);
    end

    for i = 1:length(T47_donors)
        makeHorzBarChartSimVenn(T47_donors{i}, labels, T41_donors{i}.totalBCs, folder);
    end
end

%% Calculate Relative Replicative Fitness
% RRF = (FC change in freq of BC of interest) / (FC in freq of all other BCs)
RRF_NT_raw = (donorNTNormP1.*(1-inocNormP1)) ./ (inocNormP1 .* (1-donorNTNormP1));
RRF_NT_D1_raw = RRF_NT_raw(:,[1:3,6]); % Remove 7 bc of bad sequencing
RRF_NT_D3_raw = RRF_NT_raw(:,[4:5,8:10]); % Remove 7 bc of bad sequencing
RRF_NT_D1_avg = geomean(RRF_NT_D1_raw,2);
RRF_NT_D3_avg = geomean(RRF_NT_D3_raw,2);
RRF_Trach_raw = (donorTrachNormP1.*(1-inocNormP1)) ./ (inocNormP1 .* (1-donorTrachNormP1));
RRF_Trach_D1_avg = geomean(RRF_Trach_raw(:,[1:3,6:7]),2);
RRF_Trach_D3_avg = geomean(RRF_Trach_raw(:,[4:5,8:10]),2);

% Plot RRF values sorted
figure(); hold off;
semilogy(1:nBCs,sort(RRF_NT_avg,'descend'),'LineWidth',2); hold on;
semilogy(1:nBCs,sort(RRF_Trach_avg,'descend'),'LineWidth',2);
semilogy([0 200],[sqrt(10),sqrt(10)],'k--','LineWidth',2);
semilogy([0 200],[sqrt(0.1),sqrt(0.1)],'k--','LineWidth',2);
ylim([0.1 10]); xlim([0 200]);
title('Barcode Relative Replicative Fitnesses');
ylabel('Relative Replicative Fitness');
xlabel('Rank');
legend('Nasal Turbinate','Trachea','+0.5 logs','-0.5 logs');


avgDonorRRF = avgDonor(:,3);
avgDonorRRF(1) = 0;
avgDonorRRF = avgDonorRRF ./ sum(avgDonorRRF);
for i = 2:40 % skip 1 = skip D614G
    % Calculate fold change in individual BC, all other BCs, and the 10th
    % BC (consistent, used as reference)
    FC_BC_old = avgDonorRRF(i) / inocNorm(i);
    FC_other_old = (1-avgDonorRRF(i)) / (1-inocNorm(i));
    FC_10_old = avgDonorRRF(10) / inocNorm(10);
    
    % Same as above, but for new BC mix rather than old
    FC_BC_new = avgDonorRRF(i) / newInocNorm(i);
    FC_other_new = (1-avgDonorRRF(i)) / (1-newInocNorm(i));
    FC_10_new = avgDonorRRF(10) / newInocNorm(10);
    
    % Compare FC in average donor to FC in individual donor animals
    
    RRF_avg(i,1) = FC_BC_old / FC_other_old;
    RRF_avg(i,2) = FC_BC_new / FC_other_new;
    
    RRF_ind(i,1) = FC_BC_old / FC_10_old;
    RRF_ind(i,2) = FC_BC_new / FC_10_new;
end

% Calculate geometric mean
RRF_avg(:,3) = sqrt(prod(RRF_avg(:,1:2)'));
RRF_ind(:,3) = sqrt(prod(RRF_ind(:,1:2)'));

%% Plot RRF

makeRRFFigs = input('Re-make RRF Figures? ', 's');

if upper(makeRRFFigs) == 'Y'
    makeRRFPlot(RRF_avg(2:end,:), 'Relative Replicative Fitness', folder, 0);
end
