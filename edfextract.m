%% Code that takes any length edf file and creates .dat files of 30 second epochs for viewing in ZebraExplore

filename = 'E:\Data to be analyzed\EEG recordings\2019 09 20 25\2019-09-20 topo 41 treated\export.edf';
folder = 'E:\Data to be analyzed\EEG recordings\2019 09 20 25\2019-09-20 topo 41 treated\splitedf';
name = '\eeg';
[header, data] = edfread(filename);
dtaLen = length(data);
episodes_seconds = 86400;
acq_rate = header.frequency(1);
episodes = acq_rate*episodes_seconds; 
for i = 1:(dtaLen/episodes)
    if i == 1
        st = 1;
        en = 1+episodes;
    else
        st = st+episodes;
        en = en+episodes;
        if en > dtaLen
            break
        end
    end
    datasave = data(:, st:en)';
    % [rows, columns] = size(datasave);
    if i<10
        add = '000';
    elseif i>=10 && i <100
        add = '00';
    elseif i>=100 && i<1000
        add = '0';
    else add = '';
    end
    fileOut = [folder name '_' add num2str(i) '.dat'];
    save(fileOut, 'datasave', '-ascii');
   
end
display('finished!')
