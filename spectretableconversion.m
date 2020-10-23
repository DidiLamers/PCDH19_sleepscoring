%% Program for patient EEG analysis Meyer: to convert ZebraSpecter table to table for analysis in excel

% input file should contain 4 columns only: channelnumber, BP1, BP2, and
% BP3, otherwise adapt program

%% 16 March 2020: I added the analysis to count NREM episodes and average delta power in these episodes
% Note that it is a very crude measurement: no artefacts rejected, no EMG
% included. Merely theta/delta ration or beta/delta ration below 0.75 as in
% PRISM scoring and average delta power in NREM scored episodes

% data to fill out
file = load('F:\Data to be analyzed\Meyer EEG patients\sonno controlli melani\Simiuolo_2analyzed.txt');
% wheretosave = 'F:\Data to be analyzed\Meyer EEG patients\dds2table.txt';
nch = 27;


% set parameters to work with
filen = length(file)/nch; % parameter contains the number of epochs and accordingly .dat files to be analyzed
tableout = zeros(nch*4, filen);
file_startidx = 1:nch:length(file); % vector contains the index in file of each new .dat file
for f = 1:filen
    for c = 1:nch
        row = (file_startidx(f)-1)+c;
        col = ((c-1)*4)+1;
        tableout(col, f) = file(row, 1); % fill out the channel number on the first row of each new channel
        tableout(col+1, f) = file(row, 2); %the BP1 value on the second row
        tableout(col+2, f) = file(row, 3); % and the BP2 value on the third row
        tableout(col+3, f) = file(row, 4); % and the BP3 value on the third row
    end
end

% Now add the sleep scoring
perctime_nrem_td = zeros(nch,1);
perctime_nrem_bd = zeros(nch,1);
deltapow_nrem_td = zeros(nch,1);
deltapow_nrem_bd = zeros(nch,1);
thetadelta = zeros(nch, filen);
betadelta = zeros(nch, filen);
for ii = 1:nch
    idx_nrem_td = [];
    idx_nrem_bd = [];
    c_nrem_td = [];
    c_nrem_bd = [];
    rown = ((ii-1)*4)+1;
    rdelta = rown+1;
    rtheta = rown+2;
    rbeta = rown+3;
    td = tableout(rtheta,:)./tableout(rdelta,:);    
    bd = tableout(rbeta,:)./tableout(rdelta,:);    
    idx_nrem_td = find(td <0.75);
    idx_nrem_bd = find(bd <0.75);
    c_nrem_td = length(idx_nrem_td);
    c_nrem_bd = length(idx_nrem_bd);
    perctime_nrem_td(ii) = (c_nrem_td/filen)*100; % percentage of epochs labeled as NREM by theta/delta
    % ratio <0.75
    perctime_nrem_bd(ii) = (c_nrem_bd/filen)*100;  % percentage of epochs labeled as NREM by beta/delta
    % ratio <0.75 
    deltapow_nrem_td(ii) = nanmean(tableout(rdelta,idx_nrem_td));
    deltapow_nrem_bd(ii) = nanmean(tableout(rdelta,idx_nrem_bd));
    thetadelta(ii, :) = td;
    betadelta(ii, :) = bd;
end

output1 = [thetadelta; betadelta];
output2 = [perctime_nrem_td'; perctime_nrem_bd'; deltapow_nrem_td'; deltapow_nrem_bd'];

display('finished!')
