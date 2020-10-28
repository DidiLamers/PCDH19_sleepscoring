# PCDH19_sleepscoring
 
This folder contains all the code I used during my PhD thesis on PCDH19 to score the EEG data acquired from mice with the Pinnacle setup.
It also contains my attempts at analysing patient EEG data from the Meyer hospital for sleep disturbances.

All scripts developed by me contains detailed information in their top comments.edfread.m is an already existing Matlab function (not written by me), and
so is findpeaks.m 

The sleep scoring algorithm SCOPRISM developed by Bastianini et al. 2014 was adapted to perform sleep scoring.

The script Appendix_SCOPRISMMatlab.m contains their original code, while PRISM_scoring.m contains my adaptation of it.
The script SleepScore_mainscript.m calls SCOPRISM and quantifies its output. Therefore, there is no need to run PRISM_scoring.m directly,
instead there is only need to run SleepScore_mainscript.m. It takes edf data as an input and outputs the fraction of time spent in each
sleep phase, power during each sleep phase in various frequency bands, and some sleep fragmentation information and circadian rhythm of sleep.

The folder contains some other scripts as well that have been useful sometimes. Detail descriptions can be found in the comments on top of each program. 