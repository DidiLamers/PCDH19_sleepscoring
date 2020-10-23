%% My adaptation of Bastianini et al. 2014 code. It now suits our purposes

function MSCO = PRISM_scoring(eeg,emg,samprate, plotshow, emgpeak2)

if not(isequal(size(eeg),size(emg))), error('the eeg and emg inputs have different size'); end
qseg = floor(length(eeg)/(samprate*4));  % number of 4 sec segments 
% qseg is the number of 4-s data epochs in the recordings
length_epoch = samprate*4; % length of each epoch
eeg = eeg(1:samprate*4*qseg); emg = emg(1:samprate*4*qseg);   % cut down the tail of the data


%% Initialize BP filters to be included

% These are all to be filled out by the user according to personal preference
BPnumber = 6;
BP1 = [0.5; 4]; 
BP2 = [6; 9]; 
BP3 = [8; 12];
BP4 = [12; 15];
BP5 = [15; 30];
BP6 = [30; 100];
bp_delta = 1; %here indicate which band pass to use as delta in the calculation for sleep scoring
bp_theta = 2; %here indicate which band pass to use as theta in the calculation for sleep scoring
BP_complete = [BP1 BP2 BP3 BP4 BP5 BP6];

% Then we will find all indexes belong to these frequency bands, from now
% on all calculations are automated
freq = samprate*(0:(length_epoch/2))/length_epoch; % These are the frequency values the fft transform will calculate
idxBP = zeros(2,BPnumber); % initiate a vector with two rows: first row will contain start values of each BP, 
% second row the end value of each band pass in terms of index of freq
for bp = 1:BPnumber
    [value idx_st] = min(abs(freq-BP_complete(1,bp)));
    [value idx_en] = min(abs(freq-BP_complete(2,bp)));
    idxBP(1,bp) = idx_st;
    idxBP(2,bp) = idx_en;
end

[value notch1_idx] = min(abs(freq-49)); % find the start to filter out 50 Hz (49 Hz)
[value notch2_idx] = min(abs(freq-51)); % and the end

%% Automatic data rejection based on amplitude signal

% initialize the output vector and limit for amplitude
datarejection = zeros(qseg, 1); % this will for each segment contain a 1 if included, a 0 in excluded
limit = 300;
numberrejected = 0; % this variable is going to count the number of rejected epochs

% now start looping through the segments and label each segment as accepted
% or rejected
for ss = 1:qseg
    sigstart = (ss-1)*samprate*4+1; 
    sigstop = ss*samprate*4;  
    % sigstart and sigstop are the cell positions in the eeg and emg vectors that correspond 
    %to the beginning and end of the epoch nseg, respectively.
    seg_eeg = detrend(eeg(sigstart:sigstop));
    seg_emg = detrend(emg(sigstart:sigstop));    
    thr_eeg = find(seg_eeg > limit); thr_eeg_2 = find(seg_eeg < -limit);
    thr_emg = find(seg_emg > limit); thr_emg_2 = find(seg_emg < -limit);
    if isempty(thr_eeg) && isempty(thr_eeg_2) && isempty(thr_emg) && isempty(thr_emg_2)
        datarejection(ss) = 1;
    else 
        datarejection(ss) = 0;
        numberrejected = numberrejected + 1;
    end
end
perc_rejected = (numberrejected/qseg)*100;
text = [num2str(perc_rejected),  '% of the data is rejected'];
display(text);

pause; clc
MSCO.rejecteddata = perc_rejected;

%% Initialize output matrix

MSCO.S = zeros(qseg,2).*NaN;                                 
% initializes the matrix S, where each row corresponds to one epoch, in the output structure MSCO
MSCO.POW = zeros(qseg,BPnumber).*NaN;
powertoplot = zeros(qseg, length(freq));
% initializes the matrix POW, where each row corresponds to one epoch, in the output structure 
% MSCO.
%% Now start looping through all segments, and find power in frequency bands

for nseg = 1:qseg
    if datarejection(nseg) == 0
        td = NaN;
        em = NaN;
    else     
        % nseg is a 4-s data epoch
	    sigstart = (nseg-1)*samprate*4+1; sigstop = nseg*samprate*4;  
        % sigstart and sigstop are the cell positions in the eeg and emg vectors that correspond 
        %to the beginning and end of the epoch nseg, respectively.
        seg = eeg(sigstart:sigstop);
        dtdata = detrend(seg);
        dtfft = fft(dtdata);
        pow = abs(dtfft.^2)/(length_epoch^2)/(1/4);     
	    pow = pow(1:floor(length_epoch/2)+1)*2; % I changed starting from 2 to starting from 1 % Why times two??????           
        % pow is the power spectral density of the electroencephalographic signal in the epoch nseg
        pow(notch1_idx:notch2_idx) = (pow(notch1_idx-1)+pow(notch2_idx+1))/2; % Set power 49-51 Hz to average surrounding
        powertoplot(nseg, : ) = pow';
        td = sum(pow(idxBP(1,bp_theta):idxBP(2,bp_theta)))/sum(pow(idxBP(1,bp_delta):idxBP(2, bp_delta)));             
        % td is ratio between the electroencephalographic spectral power in the theta (BP2) and delta 
        % (BP1) frequency ranges (theta/delta) computed in the epoch nseg 
	    em = sqrt(mean(detrend(emg(sigstart:sigstop),'constant').^2));  % why detrending constant here, but not in eeg power??
        % em is the root mean square (emg rms) of the zero-mean electromyographic signal in the 
        % epoch	nseg
    end
    MSCO.S(nseg,1:2) = [td em];                              
    % The values of td and em corresponding to the epoch nseg are copied to the corresponding row of
    % the matrix S    
        
    for bp = 1:BPnumber
        if datarejection(nseg) == 0
            power = NaN;
        else
            power = sum(pow(idxBP(1,bp):idxBP(2,bp)));
        end
        MSCO.POW(nseg,bp) = power;
    end
    % Now the POW struct contains on each row the 4 second segment, and on each
    % column the calculated band passes
end

emgsort = sort(MSCO.S(:,2)); % sorts the emg in ascending order
MSCO.emg5c = nanmedian(emgsort(1:round(length(emgsort)/20)));
% The variable emg5c in the output structure MSCO reflects the lowest emg rms values recorded. 
%This index is computed as the median of the 5th centile of	emg rms values.
maxnaxis = prctile(MSCO.S(:,2),99); % divides emg into percentiles, results the value of the 99th percentile 
naxisteps = maxnaxis/60; % then displays the data in 60 steps
naxis = 0:naxisteps:maxnaxis;  % creates the x axis to display                           
% The vector naxis spans the range of emg rms values from 0 to the 99th percentile of emg values 
%(variable maxnaxis) with a resolution coded in the naxisteps variable. A value of 3 for
% maxnaxis may be adequate when an electromyographic signal expressed in units of volts.
naxis2 = naxis/MSCO.emg5c;                               
% The vector naxis2 recodes naxis in units relative to the lowest emg rms values recorded. The
% variable naxis2 is useful for graphical representations because it is directly comparable among 
% experiments with different amplification of the electromyographic signal.
maxisteps =  5/100; maxis = 0:maxisteps:5;                                  
% The vector maxis spans the theta/delta ratio values from 0 to an empirically determined value of 5, which is independent from electroencephalogram amplification. 
MSCO.Z = ones(length(maxis),length(naxis)).*NaN;                                  
% Initializes the matrix Z in the output structure MSCO. The matrix Z will include the fraction of
% epochs characterized by a given interval of theta/delta ratios (coded in the maxis vector and mapped to the rows of Z) and of emg rms values (coded in the naxis vector,  and mapped to
% the columns of Z).
for n = 1:length(naxis)
        %clc;
    	%disp(['computing 3D matrix, step ',mat2str(n),' of ',mat2str(length(naxis))])
    	for m = 1:length(maxis)
        		MSCO.Z(m,n) = (sum( MSCO.S(:,2) >= naxis(n) & MSCO.S(:,2) < (naxis(n)+naxisteps) & MSCO.S(:,1) >= maxis(m) & MSCO.S(:,1) < (maxis(m)+maxisteps)))/qseg;
    	end
end
MSCO.densi = sum(MSCO.Z);   % sums all rows (theta/delta part) for each column (emg part)                                      
% The vector densi in the output structure MSCO includes the fraction of epochs characterized by
%	a given interval of emg rms values, which is coded in the naxis vector.

%% The part below is aimed at automatically identifying the cutoff values for the emg, but for my data it does not 
% work very well. So I commented it out for now and wrote my own code

% cuthyp = 18; 

% The variable cuthyp codes the first guess of the position of densi corresponding to the minimum
%	between modes of the densi variable. The value of 18 is adequate for cuthyp in cases with 2
%	modes of the densi variable. In cases with 3 modes of the densi variable, it may be
%	necessary to set the value of cuthyp so that it corresponds to a position between the highest
%	and the intermediate modes of the densi vector.

% maxlow = find(MSCO.densi(1:cuthyp)==max(MSCO.densi(1:cuthyp))); % returns index of max value up to cuthyp
% maxhigh = find(MSCO.densi(cuthyp+1:end)==max(MSCO.densi(cuthyp+1:end)))+cuthyp; % returns index of max value starting from cuthyp
% bottom = find (MSCO.densi (maxlow:maxhigh) == min(MSCO.densi(maxlow:maxhigh)),1,'first') + maxlow-1; 
% find index of minimum value between maxlow and maxhigh
% [cutlow, cuthigh] = deal(bottom); % creates cutlow and cuthigh each with value bottom
% while abs((MSCO.densi(cutlow-1)-MSCO.densi(bottom))/MSCO.densi(bottom)) < 0.20
% 	cutlow = cutlow-1; 
% end
% while abs((MSCO.densi(cuthigh+1) - MSCO.densi(bottom))/MSCO.densi(bottom)) < 0.20
% 	cuthigh = cuthigh+1; 
% end

%% This is the code I wrote to solve the problem
[pks, idx_pks] = findpeaks(MSCO.densi); % find all peaks in the data

% We are interested in the highest and seond to highest peak assuming there
% are only two modes of EMG: wake versus sleep

% So first test whether there were at least two peaks in the data

if length(pks) < 2
    error('Only 1 peak in EMG, signal not good enough to score!');
end

% Find the highest and second higest peaks and their index in the densi
% variable
sorted_pks = sort(pks);
pks_high = sorted_pks(end);
pks_high_idx_pks = find(pks == pks_high);
pks_high_idx = idx_pks(pks_high_idx_pks);
pks_second_high = sorted_pks(end-1);
pks_second_high_idx_pks = find(pks == pks_second_high);
pks_second_high_idx = idx_pks(pks_second_high_idx_pks);
if emgpeak2 == 1
    if length(pks) < 3
         error('Only 2 peaks in EMG, cannot use third highest peak for analysis');
    end
    pks_third_high = sorted_pks(end-2);
    pks_third_high_idx_pks = find(pks == pks_third_high);
    pks_third_high_idx = idx_pks(pks_third_high_idx_pks);
end

% Then find the lowest value (bottom) and its idx between the peaks
if emgpeak2 == 0
    peakone_idx = min([pks_high_idx pks_second_high_idx]);
    peaktwo_idx = max([pks_high_idx pks_second_high_idx]);
    [bottom idx_bottom] = min(MSCO.densi(peakone_idx:peaktwo_idx));
    idx_bottom = peakone_idx+idx_bottom-1;
else 
    peakone_idx = min([pks_second_high_idx pks_third_high_idx]);
    peaktwo_idx = max([pks_second_high_idx pks_third_high_idx]);
    [bottom idx_bottom] = min(MSCO.densi(peakone_idx:peaktwo_idx));
    idx_bottom = peakone_idx+idx_bottom-1;
end

% This is the final solution as discussed with gimmi: simply take the
% bottom between the two largest peaks as cutoff
emgcut = naxis(idx_bottom);
MSCO.emgcut = emgcut;

%% This part below was my solutio before discussion with gimmi, but complicated and depends on arbritrary statistical values

% Now that we have found the 'correct' bottom, we can use the original code
% agaqin to find a minimum of 20% difference between bottom and cutlow and
% cuthigh values

%[cutlow, cuthigh] = deal(idx_bottom); 
%while abs((MSCO.densi(cutlow-1)-MSCO.densi(idx_bottom))/MSCO.densi(idx_bottom)) < 0.20
% 	  cutlow = cutlow-1; 
%end
%while abs((MSCO.densi(cuthigh+1) - MSCO.densi(idx_bottom))/MSCO.densi(idx_bottom)) < 0.20
% 	cuthigh = cuthigh+1; 
%end

% Sometimes if the peaks are not separated well, or if the slope from the
% bottom to the peaks is not so steep, the cutlow or cuthigh values pass by
% the peaks which is of course wrong. In this case, I try to find the first
% data point where the increase of the signal between the bottom towards
% the peak is 20% of the distance between the bottom and the lowest peak.
% If, because there are little data points between the bottom and the peak,
% the calculated cutlow and cuthigh values turn out to be the same as the
% bottom value, I simply take the index next to the bottom to avoid too
% much overlap betweent the wake and the sleep scoring. 

%if cutlow <= peakone_idx || cuthigh >= peaktwo_idx 
%    display('could not compute emg with a 20% increase slope, instead a take a 20% increase from bottom of bottom lowest peak distance')
%    distance = pks_second_high-bottom;
%    perc_distance = (distance/100)*20;
%    value_distance = bottom+perc_distance;
%    [low cutlow] = min(abs(MSCO.densi(peakone_idx:idx_bottom)-value_distance));
%    cutlow = (cutlow+peakone_idx)-1;
%    [high cuthigh] = min(abs(MSCO.densi(idx_bottom:peaktwo_idx)-value_distance));
%    cuthigh = (cuthigh+idx_bottom)-1;    
%end

%if cutlow == idx_bottom
%    display('warning: cutlow 20% increase still at bottom, probably peaks are not separated well. Cutlow is now one index from the bottom')
%    cutlow = cutlow-1;
%end
%if cuthigh == idx_bottom
%    display('warning: cuthigh 20% increase still at bottom, probably peaks are not separated well. Cutlhigh is now one index from the bottom')
%    cuthigh = cuthigh+1;
%end

%MSCO.emgcut = [naxis(cutlow) naxis(cuthigh)];   

% The vector emgcut in the output structure MSCO includes the lower and higher values of 
%	emg rms that define a boundary zone between the 2 modes (or between the highest and the 
%	intermediate modes) of the densi variable. 

%% Now start the scoring 

auto = ones(nseg,1).*NaN; 

% Initializes the column vector auto, which will include the first step of automatic sleep scoring 
%	based only on local properties (theta/delta ratio, emg rms) of each epoch. The wake-sleep
%	state is scored as 1 (wakefulness), 2 (non-rapid-eye-movement sleep, NREMS), or 3 (rapid-
%	eye-movement sleep, REMS). Values of 10 and 20 indicate epochs with values of
%	theta/delta ratio and emg rms intermediate between wakefulness and NREMS or between
%	NREMS and REMS, respectively.
auto(MSCO.S(:,2) > MSCO.emgcut) = 1;                      
% wakefulness is scored when emg rms is above the boundary region defined by emgcut
auto(MSCO.S(:,2) < MSCO.emgcut & MSCO.S(:,1) < 0.75) = 2;      
% NREMS is scored when emg rms is below the boundary region defined by emgcut and the
%	theta/delta ratio is lower than 0.75. The latter condition indicates that
%	electroencephalographic spectral power in the delta frequency range is prevalent compared
%	with that in the theta frequency range.
auto(MSCO.S(:,2) < MSCO.emgcut & MSCO.S(:,1) > 1.25) = 3;      
% REMS is scored when emg rms is below the boundary region defined by emgcut and the
%	theta/delta ratio is higher than 1.25. The latter condition indicates that
% 	electroencephalographic spectral power in the theta frequency range is prevalent compared
%	with that in the delta frequency range.
% auto(MSCO.S(:,2) <= MSCO.emgcut(2) & MSCO.S(:,2) >= MSCO.emgcut(1)) = 10
% % The above command was removed since we used only the bottom value (one cutoff value) and
% so this piece of code cannot be used anymore
% An epoch is considered as intermediate between wakefulness and NREMS (code 10) when its
%	emg rms value is included in the boundary region defined by emgcut.
auto(MSCO.S(:,2) < MSCO.emgcut & MSCO.S(:,1) >= 0.75 & MSCO.S(:,1) <= 1.25) = 20; 
% An epoch is considered as intermediate between NREMS and REMS (code 20) when the emg rms
%	is below the boundary region defined by emgcut and theta/delta ratio is between 0.75 and
%	1.25. The latter condition indicates that electroencephalographic spectral power in the delta
%	and theta frequency ranges are nearly equivalent.
auto(datarejection==0) = 20;
% finally label all rejected data as unknown

%% Now we rescore based on surrounding epochs

auto2 = auto;                                           
% Initializes the column vector auto2, which includes the second and final step of automatic sleep
%	scoring. This step is based on local properties of each epoch, as elaborated in the auto
%	vector, but also takes into account information on the automatic scoring of adjacent epochs. 
%	Codes for wakefulness (1), NREMS (2) and REMS (3) are the same as in the vector auto. The
%	codes 10 and 20 used in the vector auto are replaced by either codes 1-3 or code 4, which
%	indicates an indeterminate state.
for n = 1:length(auto2)-3
	if auto2(n) ~= 3
        continue
    else
        %clc
	    %disp( ['phase one: progress ', mat2str(round (n/length(auto2) * 100)) , ' %'])
        if ismember(auto2(n+1),[2 20]) && auto2(n+2)==3
			auto2(n+1) = 3;
            % Single epochs scored as either NREMS or indeterminate state (codes 10 or 20) are re-scored as
            % REMS in case they are preceded and followed by at least one REMS epoch 
        elseif sum (ismember(auto2(n+1:n+2), [2 20])) == 2 && sum (auto2(n+1:n+2) == 2) <2 && auto2(n+3) == 3
			auto2(n+1:n+2) = [3;3];
            % sum (ismember(auto2(n+1:n+2), [2 20])) == 2 --> means if n+1
            % and n+2 are either 2 or 20 (not 1)
            % sum (auto2(n+1:n+2) == 2) <2 --> AND if only one of those
            % values is 2 (NREM sleep), not both
            % auto2(n+3) == 3 --> AND if n+3 is REM sleep again!
            
            % Then we will rescore n+1 and n+2 to REM sleep (3)
            
            % Couple of epochs, which are both scored as indeterminate states (codes 10 or 20), or which are 
            % scored as one epoch in indeterminate state and one epoch in NREMS, are re-scored as 2
            % epochs of REMS in case they are preceded and followed by at least one REMS epoch.
        end
    end
end
for n = 1:length(auto2)-3
	if auto2(n) ~= 2
        continue
    else
        %clc
		%disp(['phase two: progress ',mat2str(round(n/length(auto2)*100)),' %'])
        if ismember(auto2(n+1),[3 20]) && auto2(n+2)==2
			auto2(n+1) = 2;
            % Single epochs scored as either REMS or indeterminate state (codes 10 or 20) are re-scored as
            % NREMS in case they are preceded and followed by at least one NREMS epoch 
        elseif auto2(n+1)==20 && auto2(n+2)==3 && auto2(n+3)==3
			auto2(n+1) = 4;
            % Epochs scored as indeterminate state (code 20) are confirmed as such (code 4) in case they are
            % preceded by at least one NREMS epoch and followed by at least 2 REMS epochs
        elseif sum (ismember (auto2(n+1:n+2), [3 20] )) == 2 && sum(auto2(n+1:n+2)==3)<2 && auto2(n+3)==2
			auto2(n+1:n+2) = [2;2];
            % Couple of epochs, which are both scored as indeterminate states (codes 10 or 20), or which are
            % scored as one epoch in indeterminate state and one epoch in REMS, are re-scored as 
            % epochs of NREMS in case they are preceded and followed by at least one epoch of NREMS.
        end
    end
end

auto2((auto2==20)) = 2; % auto2((auto2==10)) = 4; % we no longer need this, 10 no longer exists
% Epochs that after the previous substitutions remain scored as indeterminate state are re-scored as
% NREMS, if their previous code was 20, or are confirmed as indeterminate state (code 4), if
% their previous code was 10.

for n = 1:length(auto2)-3
	if auto2(n) ~= 1
        continue
    else
        %clc
		%disp(['phase three: progress ',mat2str(round(n/length(auto2)*100)),' %'])
        if (ismember(auto2(n+1),[2 3 4]) && auto2(n+2)==1)
			auto2(n+1) = 1;
            % Epochs scored as NREMS, REMS, or indeterminate state (code 4) are re-scored as wakefulness in
            % case they are preceded by at least one epoch of wakefulness.
        elseif (auto2(n+1)==3 && auto2(n+2)==2)
			auto2(n+1) = 4;
            % Epochs scored as REMS are re-scored as indetermined state (code 4) in case they are preceded by
            % at least one epoch of wakefulness and followed by at least one epoch of NREMS
        elseif sum(ismember(auto2(n+1:n+2),[2 3 4]))==2 && auto2(n+3)==1
			auto2(n+1:n+2) = [4;4];
            % Couples of epochs, which are scored as either NREMS, REMS, or indeterminate states (code 4) in
            % any combination, are re-scored as indeterminate state (code 4) in case they are preceded and
            % followed by at least one epoch of wakefulness.
        end
    end
end

posrem = find(auto2==3);
% The vector posrem indicates the position of the cells in the vector auto2, which are scored as
% REMS (code 3)
if not(isempty(posrem))
	for n = 1:length(posrem)
        if (auto2(max(1,posrem(n)-1)) ~= 3) && (auto2(max(1,posrem(n)-1)) == auto2(min(nseg,posrem(n)+1)))
            % (auto2(max(1,posrem(n)-1)) ~= 3) --> if the value of auto2
            % before the current value was not also 3 (REM)
            % (auto2(max(1,posrem(n)-1)) == auto2(min(nseg,posrem(n)+1)))
            % --> AND if that value before rem has the same score as the
            % value after rem
            auto2(posrem(n)) = auto2(max(1,posrem(n)-1));
        end
    end
end
% If an isolated epoch is scored as REMS, and if its adjacent epochs are scored in the same state,
% and if such a state is different from REMS, that epoch is re-scored with the same state as its
% adjacent epochs.

auto((auto==20)) = 4;
% Indetermined states in the auto vector (codes 10 and 20) are attributed the code 4, for the sake of
% consistency and comparison with the auto2 vector

MSCO.S(:,3) = auto; MSCO.S(:,4) = auto2;
% Copies the vectors auto and auto2 in the third and fourth columns, respectively, of the matrix S

%% Plot figure 1

if plotshow == 1
    znew = MSCO.Z; znew(znew==0)=NaN;
    % The matrix znew is a copy of the matrix Z, but with null values substituted by NaNs for the
    % purpose of graphical representation.
    figure(1);
    set(gcf,'Renderer','zbuffer')
    surf(naxis2,maxis,log10(znew)), shading interp, colormap ('jet')
    xlabel('emgrms normalized'), ylabel('eeg theta vs. delta ratio'), zlabel('fraction of epochs')
    savefig('scoprism_fig1.fig');
    % Figure 1 plots the matrix Z as a 3D graph
end

%% plot figure 2

if plotshow == 1
    figure(2);
    powtot= zeros(3,length(freq)).*NaN;
    powtot(1,:)= mean(powertoplot(MSCO.S(:,4)==1,:)); % take rows of scoring wake of pow vector, calculate mean
    powtot(2,:)= mean(powertoplot(MSCO.S(:,4)==2,:));
    powtot(3,:)= mean(powertoplot(MSCO.S(:,4)==3,:));
    
    % powtot is a matrix built with the same logic of MSCO.POW. The power spectral density values
    %reported in the 81 columns of powtot correspond to the average values calculated during 
    % wakefulness (first row), %%NREMS (second row) or REMS (third row).

    [tmp maxasX] = min(abs(freq-100));
    asX = freq(1:maxasX);
    plot(asX,powtot(1,1:maxasX),'r',asX,powtot(2,1:maxasX),'g',asX,powtot(3,1:maxasX),'b')
    legend('wake', 'NREM', 'REM');
    xlabel('Frequency (Hz)'); ylabel('EEG Power Spectral Density (units^2/Hz)');
    % Figure 2 plots the mean power spectral density (units^2/Hz) during wakefulness (red line),
    % NREMS (green line) and REMS (blue line) in the frequency range 0-200 Hz. 
    savefig('scoprism_fig2.fig');
end

%% plot figure 3

if plotshow == 1
    figure(3);   
    plot(naxis2,MSCO.densi,'rs',naxis2(idx_bottom),MSCO.densi(idx_bottom),'ks')
    xlabel('emgrms normalized'), ylabel('fraction of epochs')  
    savefig('scoprism_fig3.fig');  
    display('All figures are plotted: do you want to continue?')
    pause;
    close all;   
    % Figure 3 plots the matrix densi as a 2D graph, and marks the limits of the boundary zone between
    % the higher and lower modes of the densi variable. If the variable densi has 3 modes, users
    % should check on figure 3 whether the boundary zone is correctly positioned between the
    % highest mode and the intermediate mode. If this is not the case, the whole scoring routine
    % should be repeated after setting the value of the variable cuthyp to any position of the densi 
    % vector that is intermediate between the positions of the highest and the intermediate modes.
end

%% Validation of the data

if plotshow == 1
    posrem = find(MSCO.S(:,4)==3);
    % The vector posrem indicates the position of the rows in the matrix S, which are scored as REMS
    % (coded 3 in column 4, corresponding to the final step of automatic scoring).
    qremtd = length(find(MSCO.S(posrem,1) > 2.5))/length(posrem);
    % The variable qremtd is the fraction of REMS epochs with a theta/delta ratio larger than 2.5
    if qremtd < 0.1
        disp('WARNING: less than 10% of REMS epochs have a theta/delta EEG power ratio 	higher than 2.5');
        disp('this indicates that the quality of the raw EEG signal is insufficient for automatic sleep scoring')
        disp('Manual sleep scoring is thus recommended for this record')
        disp('Do you want to continue anyways?')
        pause;
    end
end
display('scoring complete!')
