%% To find sleepstage of beta oscillations

% It takes as an input the results from the SleepScore_mainscript.m
% program, and an excel file containing on each sheet one animal with beta
% oscillations. Each sheet consists of a row vector with the start times of
% all detected beta oscillations (as found by ZebraExplore). 

% The output is the percentage of beta oscillations in NREM sleep. 

% Note, I do this only for the injected hemisphere,and for treated mice.
% In all cases with beta oscillations, Ch1 was the injected hemisphere. So
% this script only analyses Ch1. 

%% Start by creating input files

% The sleep scoring data is on my hard disk
output_scoprism = ["F:\Data to be analyzed\EEG recordings\2019 05 06-12\2019-05-07 mouse 111 treated\figures\output.mat";
    "F:\Data to be analyzed\EEG recordings\2019 05 06-12\2019-05-10 mouse 112 treated\figures\output.mat";
    "F:\Data to be analyzed\EEG recordings\2019 06 02 07\2019-06-04 mouse 111 treated\figures\output.mat";
    "F:\Data to be analyzed\EEG recordings\2019 09 20 25\2019-09-20 topo 41 treated\figures\output.mat";
    "F:\Data to be analyzed\EEG recordings\2019 09 20 25\2019-09-24 topo 38 treated\figures\output.mat";
    "F:\Data to be analyzed\EEG recordings\2020 03 09 19\figures topo60\output.mat"];

animals = length(output_scoprism);

osc_file = "C:\Users\Didi\Desktop\beta.xlsx"; % refers to the excel file containing all beta osc start times
% Each sheet contains a new animal

%% Now loop over all animals to find their beta oscillations in NREM

perc = zeros(animals, 1); % initiate the output file with only zeroes

for animal_number = 1:animals
    load(output_scoprism(animal_number));
    [events] = xlsread(osc_file, animal_number); % events is now a row vector with the start time of each event
    number_oscillations = length(events);

    nrem = 0;
    for ev = 1:number_oscillations % loop through all events
        st_time = events(ev); % the start time of the event
        st_epoch = floor(st_time/4); % find the corresponding epoch
        sleepstage_vector = SleepScore(1).MSCO.S(st_epoch:st_epoch+3,4); % we take the start and next 3 epochs
        if length(find(sleepstage_vector == 2)) > 1 % if at least one of the epochs has been scored as NREM, i.e. 2            
            nrem = nrem + 1;
        end 
    end
    perc(animal_number) = (nrem/number_oscillations)*100;
end

mean_perc = mean(perc)






