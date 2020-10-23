%% April 2019 by Didi: script to  automatically analyse EEG recordings, sleep scoring and power analysis

% The program is based on the SCOPRISM algorithm (Bastiani et al. 2014). 
% It takes as an input an edf file. Calls the PRISM_scoring function to
% score sleep stages as in Bastianini et al. 2014.

% It then provides the output in terms of fraction of time spent in each
% sleep phase, power during each sleep phase in various frequency bands,
% and some sleep fragmentation information and circadian rhythm of sleep.

% The output is saved in a struct called SleepScore, as some information is
% copied to clipboard for easy transfer to an excel file. 

%% Note: if you want to change the BP values: they are specified within the PRISM_scoring function

%% Indication if less than 10% of REM episodes have a theta/delta ration >2.5 is explained in Bastianini et al:

% "we found that in 6 mice SCOPRISM had REMS sensitivity that was more than 2 standard deviations lower than the population
% mean. We then investigated whether in these outlier mice, the 3D plots displaying the fraction of 4-s epochs with a 
% given combination of EEG theta/delta ratio and EMG rms had unusual features. In particular, in each of these 6 mice, the 
% ridge positioned at low EMG rms values and extending toward high values of the EEG theta/delta ratio was unusually short. 
% This reflected unusually abundant slow EEG activity during REMS. Accordingly, in each of these 6 outlier mice, <10% of 
% automatically-scored REMS epochs had an EEG theta/delta ratio >2.5. This was also the case in 2 additional WT mice whose 
% REMS sensitivity, although low differed from the population mean by less than 2 standard deviations. Based on these results, 
% we elected conservatively to employ this criterion (namely, <10% of automatically-scored REMS epochs with an EEG theta/delta
% ratio >2.5) as an exclusion criterion for the correct application of the
% automatic sleep-scoring algorithm. "

%% Provide necessary input
clc; clear;
filename = 'E:\Data to be analyzed\EEG recordings\2020 06\2020-07-13 6pko bis.edf';
savefigname = 'E:\Data to be analyzed\EEG recordings\2020 06\mouse 6 rep2 figures'; % folder where you want to save the figures of the analysis
nCh = 2; % Fill out the number of eeg channels
emg = 3; % Fill out which channel contains the emg signal
samprate = 400; % Fill out the sampling rate in Hz (should be the same for eeg and emg signal)
plotshow = 1; % 1 if you want to look at data during process and do data validation, 0 if not

%% Read the file and perform automated sleep scoring with adapted SCOPRISM

[header, data] = edfread(filename); % header stores information about edf file, data contains the actual data
clc
emg = data(emg, :)'; % this will be the emg input to the SCOPRISM function, it needs to be a column vector
SleepScore(nCh) = struct(); % Initialize struct that will contain the output
for ch = 1:nCh %automated sleep scoring for both channels
    clc    
    disp(['Processing channel:',num2str(ch)])
    emgpeak2 = 0; 
    eeg = data(ch, :)'; % the eeg input to the SCOPRISM function: again it needs to be a column vector
    MSCO = PRISM_scoring(eeg, emg, samprate, plotshow, emgpeak2); % run the automated sleep scoring function
    prompt = 'EMG scoring okay? Press 0 if okay, 1 to shift to the second and third highest peak';
    emgpeak2 = input(prompt);
    if emgpeak2 == 1
        MSCO = PRISM_scoring(eeg, emg, samprate, plotshow, emgpeak2);
    end
    SleepScore(ch).MSCO = MSCO; % and save the output in the SleepScore struct 
    if plotshow == 1
        mkdir(savefigname, ['Ch' num2str(ch)]); % create folder to save the figures
        newdir = [savefigname '\Ch' num2str(ch)]; % find its name
        movefile('scoprism_fig1.fig', newdir); % and move the saved figures to the specified folder
        movefile('scoprism_fig2.fig', newdir);
        movefile('scoprism_fig3.fig', newdir);
    end
end

clc

%% find fractions spent in wake, rem, and NREM according to EEG/EMG Scoring algorithm from Bastianini et al. 2014

for ch = 1:nCh
    num_wake = length(find(SleepScore(ch).MSCO.S(:,4) == 1)); % wake is scored as '1' in the PRISM_scoring function
    frac_wake = num_wake/length(SleepScore(ch).MSCO.S(:,4));
    num_nrem = length(find(SleepScore(ch).MSCO.S(:,4) == 2)); % NREM sleep is scored as '2' in the PRISM_scoring function
    frac_nrem = num_nrem/length(SleepScore(ch).MSCO.S(:,4));
    num_rem = length(find(SleepScore(ch).MSCO.S(:,4) == 3)); % REM sleep is scored as '3' in the PRISM_scoring function
    frac_rem = num_rem/length(SleepScore(ch).MSCO.S(:,4));
    
    % save in output struct
    SleepScore(ch).fractions.WAKE = frac_wake;
    SleepScore(ch).fractions.NREM = frac_nrem;
    SleepScore(ch).fractions.REM = frac_rem;    
end

%% Plot the fractions found above

plotvector_fraction = zeros(3, nCh);
for ch = 1:nCh
    plotvector_fraction(1, ch) =  SleepScore(ch).fractions.WAKE*100;
    plotvector_fraction(2, ch) =  SleepScore(ch).fractions.NREM*100;
    plotvector_fraction(3, ch) =  SleepScore(ch).fractions.REM*100;
end  
   
    
figure(1)
xaxis = categorical({'WAKE','NREM','REM'});
b = bar(xaxis, plotvector_fraction);
b(1).FaceColor = 'c';
b(2).FaceColor = 'k';
ylabel('Episodes (%)');
ylim([0 100])
legend('Location', 'northwest');
savefig('fig_fractions.fig');
movefile('fig_fractions.fig', savefigname);
disp('Figure saved: do you want to continue?')
pause;
clc; close;

%% Plot powers in different frequency ranges 

[nseg, BPnumber] = size(SleepScore(1).MSCO.POW);
averagepower_wake = zeros(2, BPnumber);
averagepower_nrem = zeros(2, BPnumber);
averagepower_rem = zeros(2, BPnumber);

for ch = 1:nCh
    idx_wake = find(SleepScore(ch).MSCO.S(:,4) == 1);
    idx_nrem = find(SleepScore(ch).MSCO.S(:,4) == 2);
    idx_rem = find(SleepScore(ch).MSCO.S(:,4) == 3);
    totalpower = sum(nanmean(SleepScore(ch).MSCO.POW(:,:)));
    for bp = 1:BPnumber
        averagepower_wake(ch, bp) = nanmean(SleepScore(ch).MSCO.POW(idx_wake,bp))/totalpower;
        averagepower_nrem(ch, bp) = nanmean(SleepScore(ch).MSCO.POW(idx_nrem,bp))/totalpower;
        averagepower_rem(ch, bp) = nanmean(SleepScore(ch).MSCO.POW(idx_rem,bp))/totalpower;
    end
end

% Save the average powers
for ch = 1:nCh
    SleepScore(ch).power.WAKE = averagepower_wake(ch, :);
    SleepScore(ch).power.NREM = averagepower_nrem(ch, :);
    SleepScore(ch).power.REM = averagepower_rem(ch, :);
end

figure(1)
if BPnumber == 6
    xaxiswake = categorical({'Delta','Theta','Alpha','Sigma','Beta','Gamma'});
    xaxiswake = reordercats(xaxiswake, {'Delta', 'Theta', 'Alpha', 'Sigma', 'Beta', 'Gamma'});
    b = bar(xaxiswake, averagepower_wake');
else 
    b = bar(averagepower_wake');
end
legend;
savefig('fig_wakepower.fig');

figure(2)
if BPnumber == 6
    xaxisnrem = categorical({'Delta','Theta','Alpha', 'Sigma', 'Beta', 'Gamma'});
    xaxisnrem = reordercats(xaxisnrem, {'Delta', 'Theta', 'Alpha', 'Sigma', 'Beta', 'Gamma'});
    b = bar(xaxisnrem, averagepower_nrem');
else 
    b = bar(averagepower_nrem');
end
legend;
savefig('fig_nrempower.fig');

figure(3)
if BPnumber == 6
    xaxisrem = categorical({'Delta','Theta','Alpha', 'Sigma', 'Beta', 'Gamma'});
    xaxisrem = reordercats(xaxisrem, {'Delta', 'Theta', 'Alpha', 'Sigma', 'Beta', 'Gamma'});
    b = bar(xaxisrem, averagepower_rem');
else 
    b = bar(averagepower_rem');
end
legend;
savefig('fig_rempower.fig');

movefile('fig_wakepower.fig', savefigname);
movefile('fig_nrempower.fig', savefigname);
movefile('fig_rempower.fig', savefigname);
disp('Figures saved: do you want to continue?')
pause;
clc; close all;


%% Fractionation of sleep
% First plot a histogram of duration epochs


for state = 1:3
    for ch = 1:nCh
        countwake = 0;
        wakelength = 0;
        macount = 0;
        for epoch = 1:length(SleepScore(ch).MSCO.S(:,4))
            if SleepScore(ch).MSCO.S(epoch,4)== state
                if epoch == 1
                    countwake = countwake + 1;
                    wakelength = wakelength + 1;
                elseif SleepScore(ch).MSCO.S(epoch-1,4)~= state
                    countwake = countwake+1;
                    wakelength = wakelength+1;
                else
                    wakelength = wakelength+1;
                end
                if epoch == length(SleepScore(ch).MSCO.S(:,4))
                    if state == 1
                        SleepScore(ch).epochlength_wake(countwake) = wakelength;
                    elseif state == 2 
                        SleepScore(ch).epochlength_nrem(countwake) = wakelength;
                    elseif state == 3
                        SleepScore(ch).epochlength_rem(countwake) = wakelength;
                    end
                    wakelength = 0;
                elseif SleepScore(ch).MSCO.S(epoch+1,4)~= state
                    if state == 1
                        if countwake > 1
                            if wakelength < 10 && SleepScore(ch).MSCO.S(epoch+1,4)==2 ...
                                    && SleepScore(ch).MSCO.S(epoch-wakelength,4)==2
                                macount = macount+1;
                                maduration(macount) = wakelength;
                                SleepScore(ch).epochlength_wake(countwake) = NaN;
                            else
                                SleepScore(ch).epochlength_wake(countwake) = wakelength;
                            end
                        else
                            SleepScore(ch).epochlength_wake(countwake) = wakelength;
                        end
                    elseif state == 2 
                        SleepScore(ch).epochlength_nrem(countwake) = wakelength;
                    elseif state == 3
                        SleepScore(ch).epochlength_rem(countwake) = wakelength;
                    end
                    wakelength = 0;
                end
            end
        end
        if state == 1
            ma(ch) = macount; % frequency of microarousals per hour
            madurations(ch) = mean(maduration)*4; % average duration ma in s
        else % do nothing
        end
    end
end

% Change the MA count to a frequency: per time spend in sleep

rectime = nseg*4;
for ch = 1:nCh
    fractionsleep = SleepScore(ch).fractions.NREM + SleepScore(ch).fractions.REM;
    sleeptime = rectime*fractionsleep;    
    SleepScore(ch).ma.freqtottime = (ma(ch)/rectime)*60; % frequency of ma in total recording per minute
    SleepScore(ch).ma.freqsleep =  (ma(ch)/sleeptime)*60; % frequency of ma in sleep time per minute
    SleepScore(ch).ma.durations = madurations(ch); 
end


% Save the output metrics: precentage epochs shorter than 5 epochs, and
% those lasting longer than 30. 

for ch = 1:nCh
    
    % For wake
    wakeshort = length(find(SleepScore(ch).epochlength_wake<=5));
    SleepScore(ch).fragm.perc_wake_short = (wakeshort/length(SleepScore(ch).epochlength_wake));
    wakelong = length(find(SleepScore(ch).epochlength_wake>=30));
    SleepScore(ch).fragm.perc_wake_long = (wakelong/length(SleepScore(ch).epochlength_wake));
    SleepScore(ch).fragm.av_bout_duration_wake = nanmean(SleepScore(ch).epochlength_wake);
    
    % NREM
    nremshort = length(find(SleepScore(ch).epochlength_nrem<=5));
    SleepScore(ch).fragm.perc_nrem_short = (nremshort/length(SleepScore(ch).epochlength_nrem));
    nremlong = length(find(SleepScore(ch).epochlength_nrem>=30));
    SleepScore(ch).fragm.perc_nrem_long = (nremlong/length(SleepScore(ch).epochlength_nrem));
    SleepScore(ch).fragm.av_bout_duration_nrem = nanmean(SleepScore(ch).epochlength_nrem);
    
    % REM
    remshort = length(find(SleepScore(ch).epochlength_rem<=5));
    SleepScore(ch).fragm.perc_rem_short = (remshort/length(SleepScore(ch).epochlength_rem));
    remlong = length(find(SleepScore(ch).epochlength_rem>=30));
    SleepScore(ch).fragm.perc_rem_long = (remlong/length(SleepScore(ch).epochlength_rem));
    SleepScore(ch).fragm.av_bout_duration_rem = nanmean(SleepScore(ch).epochlength_rem);
    
end

% plot the histogram

% save the values to be plotted
edges = [0:30 Inf];

figure(1)
for ch = 1:nCh
    hold on;
    ht = histogram(SleepScore(ch).epochlength_wake, edges);   
    legend
end

savefig('fig_wakelength.fig');

figure(2)
for ch = 1:nCh
    hold on;
    ht = histogram(SleepScore(ch).epochlength_nrem, edges);   
    legend
end

savefig('fig_nremlength.fig');

figure(3)
for ch = 1:nCh
    hold on;
    ht = histogram(SleepScore(ch).epochlength_rem, edges);   
    legend
end

savefig('fig_remlength.fig');

movefile('fig_wakelength.fig', savefigname);
movefile('fig_nremlength.fig', savefigname);
movefile('fig_remlength.fig', savefigname);
disp('Figures saved: do you want to continue?')
pause;
clc; close all;

%% Final part: sleep in day versus night

% find starttime:
day = str2double(header.startdate(1:2));
month = str2double(header.startdate(4:5));
year = str2double(header.startdate(7:8))+2000;

hour = str2double(header.starttime(1:2));
min = str2double(header.starttime(4:5));
sec = str2double(header.starttime(7:8));

dt = datetime(year, month, day, hour, min, sec); % dt now contains the start date and time of recording
disp('Is this the correct start date?:') % quick check if it makes sense
disp(dt)
pause;

lightout = datetime(year, month, day, 19, 0, 0); % all recordings start in light on cycle, we find the
                                    % time at which the light goes out
endrec = dt + seconds(nseg*4); % and the end of the recording
cycles2 = lightout:hours(12):endrec; % We define all the cycles: cycle 1 is start till light out, 
                               % additional cycles are 12 hour cycles,
                               % cycle 3 is from the last 12 hour cycle
                               % till the end of recording

durationcycle1 = lightout-dt;
durationcycles2 = hours(12);
ncycles2 = length(cycles2)-1;
if ncycles2 < 1
    error('recordings appears to be too short: no 12 hour cycle');
end
durationcycle3 = endrec-cycles2(end);

% Now to find the start and end epoch index of each cycle

% We will create two vectors, one containing al the start epoch numbers,
% and one all the end epoch numbers. Epochs are 4 seconds, so we divide the
% duration in seconds by 4

starts = zeros(ncycles2+2,1);
ends = zeros(ncycles2+2,1);
starts(1) = 1;
starts(2) = starts(1)+ceil((seconds(durationcycle1)/4));
if ncycles2>1
    for cy = 3:ncycles2+1
        starts(cy) = starts(cy-1)+ceil((seconds(durationcycles2)/4));
    end
end
starts(ncycles2+2) = starts(ncycles2+1)+ceil((seconds(durationcycles2)/4));

% Now for the end epoch numbers: that is simply the next start epoch -1
for cy=1:ncycles2+2
    if cy == ncycles2+2
        ends(cy) = nseg;
    else 
        ends(cy) = starts(cy+1)-1;
    end
end
    
% Now couple to scoring day and night
labels = zeros(ncycles2+2,1);
labels(1:2:end) = 1; % now all wake cycles are labeled one: we start recording during day. 
frac_wake_day = zeros(length(find(labels==1)), nCh);
frac_nrem_day = zeros(length(find(labels==1)), nCh);
frac_rem_day = zeros(length(find(labels==1)), nCh);
frac_wake_night = zeros(length(find(labels==0)), nCh);
frac_nrem_night = zeros(length(find(labels==0)), nCh);
frac_rem_night = zeros(length(find(labels==0)), nCh);

for ch = 1:nCh
    countday = 0;
    countnight = 0;
    for cy = 1:ncycles2+2
        st_epoch = starts(cy);
        en_epoch = ends(cy);
        if labels(cy) == 1 % it is day
            countday = countday+1;
            wake = length(find(SleepScore(ch).MSCO.S(st_epoch:en_epoch,4) == 1));             
            frac_wake_day(countday,ch) = wake/length(SleepScore(ch).MSCO.S(st_epoch:en_epoch,4));
            nrem = length(find(SleepScore(ch).MSCO.S(st_epoch:en_epoch,4) == 2));             
            frac_nrem_day(countday,ch) = nrem/length(SleepScore(ch).MSCO.S(st_epoch:en_epoch,4));
            rem = length(find(SleepScore(ch).MSCO.S(st_epoch:en_epoch,4) == 3));             
            frac_rem_day(countday,ch) = rem/length(SleepScore(ch).MSCO.S(st_epoch:en_epoch,4));
        else % it is night
            countnight = countnight+1;
            wake = length(find(SleepScore(ch).MSCO.S(st_epoch:en_epoch,4) == 1));             
            frac_wake_night(countnight,ch) = wake/length(SleepScore(ch).MSCO.S(st_epoch:en_epoch,4));
            nrem = length(find(SleepScore(ch).MSCO.S(st_epoch:en_epoch,4) == 2));             
            frac_nrem_night(countnight,ch) = nrem/length(SleepScore(ch).MSCO.S(st_epoch:en_epoch,4));
            rem = length(find(SleepScore(ch).MSCO.S(st_epoch:en_epoch,4) == 3));             
            frac_rem_night(countnight,ch) = rem/length(SleepScore(ch).MSCO.S(st_epoch:en_epoch,4));
        end
    end
end

% then we can save the values
for ch = 1:nCh
    SleepScore(ch).day.wake = nanmean(frac_wake_day(:,ch));
    SleepScore(ch).day.nrem = nanmean(frac_nrem_day(:,ch));
    SleepScore(ch).day.rem = nanmean(frac_rem_day(:,ch)); 
    SleepScore(ch).night.wake = nanmean(frac_wake_night(:,ch));
    SleepScore(ch).night.nrem = nanmean(frac_nrem_night(:,ch));
    SleepScore(ch).night.rem = nanmean(frac_rem_night(:,ch)); 
end


%% Save the output struct

% copy some stuff to the clipboard for easy analysis
table = sprintf('%f\t', SleepScore(1).ma.freqtottime, SleepScore(2).ma.freqtottime, SleepScore(1).ma.freqsleep, ...
    SleepScore(2).ma.freqsleep, SleepScore(1).ma.durations, SleepScore(2).ma.durations);
clipboard('copy', table);
disp('Ma table copied: ready to continue?')
pause;
clc;

tablewake = sprintf('%f\t', SleepScore(1).fragm.perc_wake_short, SleepScore(2).fragm.perc_wake_short, ...
    SleepScore(1).fragm.perc_wake_long, SleepScore(2).fragm.perc_wake_long,...
    SleepScore(1).fragm.av_bout_duration_wake, SleepScore(2).fragm.av_bout_duration_wake);

tablenrem = sprintf('%f\t', SleepScore(1).fragm.perc_nrem_short, SleepScore(2).fragm.perc_nrem_short, ...
    SleepScore(1).fragm.perc_nrem_long, SleepScore(2).fragm.perc_nrem_long,...
    SleepScore(1).fragm.av_bout_duration_nrem, SleepScore(2).fragm.av_bout_duration_nrem);

tablerem = sprintf('%f\t', SleepScore(1).fragm.perc_rem_short, SleepScore(2).fragm.perc_rem_short, ...
    SleepScore(1).fragm.perc_rem_long, SleepScore(2).fragm.perc_rem_long,...
    SleepScore(1).fragm.av_bout_duration_rem, SleepScore(2).fragm.av_bout_duration_rem);

clipboard('copy', tablewake);
disp('Fragmentation wake copied: ready to continue?')
pause;
clc;
clipboard('copy', tablenrem);
disp('Fragmentation nrem copied: ready to continue?')
pause;
clc;
clipboard('copy', tablerem);
disp('Fragmentation rem copied: ready to continue?')
pause;
clc;

tabledaynight = sprintf('%f\t', SleepScore(1).day.wake, SleepScore(2).day.wake, SleepScore(1).night.wake,...
    SleepScore(2).night.wake, SleepScore(1).day.nrem, SleepScore(2).day.nrem, SleepScore(1).night.nrem,...
    SleepScore(2).night.nrem, SleepScore(1).day.rem, SleepScore(2).day.rem, SleepScore(1).night.rem,...
    SleepScore(2).night.rem);
clipboard('copy', tabledaynight);
disp('Daynight table copied: ready to continue?')
pause;
clc;

save([savefigname '\output.mat'],'SleepScore');