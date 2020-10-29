# PCDH19_sleepscoring

## General Information 
This folder contains all the code I used during my PhD thesis on PCDH19 to score the EEG data acquired from mice with the Pinnacle setup.
It also contains my attempts at analysing patient EEG data from the Meyer hospital for sleep disturbances.

## Use in PCDH19 project
The main program, SleepScore_mainscript.m, was used to perform sleep scoring of mouse EEG data, measured with the Pinnacle setup. 
Additionally, I used the script event_in_sleepphase.m to find whether detected beta oscillations in the EEG recordings (as detected with ZebraExplore) occurred
in NREM sleep episodes. 

## Usage
All scripts developed by me contain detailed information in their top comments.

The sleep scoring algorithm SCOPRISM developed by Bastianini et al. 2014 was adapted to perform sleep scoring.

The script SleepScore_mainscript.m calls SCOPRISM and quantifies its output. Therefore, there is no need to run PRISM_scoring.m directly,
instead there is only need to run SleepScore_mainscript.m. It takes edf data as an input and outputs the fraction of time spent in each
sleep phase, power during each sleep phase in various frequency bands, and some sleep fragmentation information and circadian rhythm of sleep.

The program event_in_sleepphase.m takes as an input an excel file containing the start time (in seconds after the start of recording) of detected 
beta oscillations and the output struct of SleepScore_mainscript.m and outputs the percentage of beta oscillations that occured in NREM sleep episodes. 

## Contributors
The sleep scoring code is based on Bastianini et al. 2014. Their original code is contained in the file Appendix_SCOPRISMMatlab.m, while PRISM_scoring.m
contains my adaptation of it. 
edfread.m is an already existing Matlab function (not written by me), and so is findpeaks.m.
All other code in the repository has been written by me. 

## Dependencies
None