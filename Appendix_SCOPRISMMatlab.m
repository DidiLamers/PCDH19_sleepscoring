%% Original code by Bastianini et al. 2014 Not used in my code

% PRISM_scoring.m
% July 8th, 2014

% TECHNICAL NOTE:
% This routine is aimed at clarifying the computational logic of the algorithm for automatic 	scoring of mouse sleep with a simplified Matlab implementation. 
% There are 3 input variables (eeg, emg, and samprate) and 1 output variable (MSCO).
% The following assumptions are made on the input variables: 
%	- eeg and emg include synchronized electroencephalographic and electromyographic signals, %	respectively, in the whole recording period.
%	- eeg and emg are structured as column vectors, each cell of which corresponds to one
%	sampling period (e.g., to 10 ms if the sampling rate is 100 Hz). 	
%        the length of the eeg and emg vectors is not limited by this routine. These vectors may be 
%	built starting, e.g., from .txt files using the “importdata” function or from .xls files using 
%        the “xlsread” function (built-in Matlab functions). The only limitation in applying this routine
%	concerns the computational power of the machine by which it is run. 
% 	- samprate is the sampling rate (in Hz) of the electroencephalographic and 
%	electromyographic signals, which is assumed to be 	the same.
%	- the duration of electroencephalographic and electromyographic recordings is assumed
%	 to be the same integer multiple of 4 s, which is the time resolution of sleep scoring
% 	in the present algorithm. This 4-s time unit will be referred to as epoch in the following 
% 	notes.

% TERMS OF USE:
% Users are free to adapt and apply this routine for research purposes provided that the source of 
%	the original routine and underlying algorithm is properly acknowledged by citation of this 
%	paper. 

% DISCLAIMER:
% This software routine is experimental in nature and is provided without any warranty of
%	merchantability or fitness for any purpose, or any other warranty expressed or implied. The %	Authors of this routine shall not be liable for any loss, claim or demand by any other party, %	due to or arising from the use of this software routine.

function MSCO = PRISM_scoring(eeg,emg,samprate)

if not(isequal(size(eeg),size(emg))), error('the eeg and emg inputs have different size'); end
qseg = floor(length(eeg)/(samprate*4));                        
% qseg is the number of 4-s data epochs in the recordings
eeg = eeg(1:samprate*4*qseg); emg = emg(1:samprate*4*qseg);
MSCO.S = zeros(qseg,2).*NaN;                                 
% initializes the matrix S, where each row corresponds to one epoch, in the output structure MSCO
MSCO.POW = zeros(qseg,81).*NaN;
% initializes the matrix POW, where each row corresponds to one epoch, in the output structure 
% MSCO.   
for nseg = 1:qseg                 
% nseg is a 4-s data epoch
	sigstart = (nseg-1)*samprate*4+1; sigstop = nseg*samprate*4;  
% sigstart and sigstop are the cell positions in the eeg and emg vectors that correspond 
%to the beginning and end of the epoch nseg, respectively.
    	pow = (abs(fft(detrend(eeg(sigstart:sigstop)))).^2)/((samprate*4)^2)/(1/4); 
	pow = pow(2:floor(samprate*4/2)+1)*2; 
% pow is the power spectral density of the electroencephalographic signal in the epoch nseg
    	td = sum(pow(24:36))/sum(pow(2:16));                
% td is ratio between the electroencephalographic spectral power in the theta (6-9 Hz) and delta 
%	(0.5-4 Hz) frequency ranges (theta/delta) computed in the epoch nseg 
	em = sqrt(mean(detrend(emg(sigstart:sigstop),'constant').^2));  
% em is the root mean square (emg rms) of the zero-mean electromyographic signal in the 
%	epoch	nseg
    	MSCO.S(nseg,1:2) = [td em];                              
% The values of td and em corresponding to the epoch nseg are copied to the corresponding row of
%	the matrix S
% MSCO.POW(nseg,1:81) = [[pow]',sum(pow)]; %% This does not work: pow is 800, not 80: pow(1:10:end) works!!
  MSCO.POW(nseg, 1:81) = [[pow(1:10:end)]', sum(pow)];
% The first 80 columns of the matrix POW include the power spectral density of the 
% electroencephalogram in the epoch nseg at frequencies from 0.25 Hz to 20 Hz (extremes 
% included) at steps of 0.25 Hz. The 81st column is the total power spectral density of the 
% electroencephalogram in the epoch nseg in the frequency range 0.25-20 Hz. 
end

emgsort = sort(MSCO.S(:,2)); MSCO.emg5c = median(emgsort(1:round(length(emgsort)/20)));
% The variable emg5c in the output structure MSCO reflects the lowest emg rms values recorded. 
%This index is computed as the median of the 5th centile of	emg rms values.
maxnaxis = prctile(MSCO.S(:,2),99);
naxisteps = maxnaxis/60; 
naxis = 0:naxisteps:maxnaxis;                             
% The vector naxis spans the range of emg rms values from 0 to the 99th percentile of emg values 
%(variable maxnaxis) with a resolution coded in the naxisteps variable. A value of 3 for
% maxnaxis may be adequate when an electromyographic signal expressed in units of volts.
naxis2 = naxis/MSCO.emg5c;                               
% The vector naxis2 recodes naxis in units relative to the lowest emg rms values recorded. The
% variable naxis2 is useful for graphical representations because it is directly comparable among 
% experiments with different amplification of the electromyographic signal.
maxisteps =  5/100; maxis = 0:maxisteps:5;                                  
% The vector maxis spans the theta/delta ratio values from 0 to an empirically determined value of %	5, which is independent from electroencephalogram amplification. 
MSCO.Z = ones(101,61).*NaN;                                  
% Initializes the matrix Z in the output structure MSCO. The matrix Z will include the fraction of
%	epochs characterized by a given interval of theta/delta ratios (coded in the maxis vector and %	mapped to the rows of Z) and of emg rms values (coded in the naxis vector,  and mapped to
%	the columns of Z).
for n = 1:length(naxis)
    	clc, disp(['computing 3D matrix, step ',mat2str(n),' of ',mat2str(length(naxis))])
    	for m = 1:length(maxis)
        		MSCO.Z(m,n) = (sum(MSCO.S(:,2) >= naxis(n) & MSCO.S(:,2) < (naxis(n)+naxisteps) & MSCO.S(:,1) >= maxis(m) & MSCO.S(:,1) < (maxis(m)+maxisteps)))/nseg;
    	end
end
MSCO.densi = sum(MSCO.Z);                                         
% The vector densi in the output structure MSCO includes the fraction of epochs characterized by
%	a given interval of emg rms values, which is coded in the naxis vector.
cuthyp = 18;                                            
% The variable cuthyp codes the first guess of the position of densi corresponding to the minimum
%	between modes of the densi variable. The value of 18 is adequate for cuthyp in cases with 2
%	modes of the densi variable. In cases with 3 modes of the densi variable, it may be
%	necessary to set the value of cuthyp so that it corresponds to a position between the highest
%	and the intermediate modes of the densi vector.
maxlow = find(MSCO.densi(1:cuthyp)==max(MSCO.densi(1:cuthyp))); 
maxhigh = find(MSCO.densi(cuthyp+1:end)==max(MSCO.densi(cuthyp+1:end)))+cuthyp;
bottom = find (MSCO.densi (maxlow:maxhigh) == min(MSCO.densi(maxlow:maxhigh)),1,'first') + maxlow-1; 
[cutlow, cuthigh] = deal(bottom);
while abs((MSCO.densi(cutlow-1)-MSCO.densi(bottom))/MSCO.densi(bottom)) < 0.20
	cutlow = cutlow-1; 
end
while abs((MSCO.densi(cuthigh+1) - MSCO.densi(bottom))/MSCO.densi(bottom)) < 0.20
	cuthigh = cuthigh+1; 
end
MSCO.emgcut = [naxis(cutlow) naxis(cuthigh)];           
% The vector emgcut in the output structure MSCO includes the lower and higher values of 
%	emg rms that define a boundary zone between the 2 modes (or between the highest and the 
%	intermediate modes) of the densi variable. The values of emgcut are defined arbitrarily as
%	the positions in the densi vector, which are nearest to the position defined by bottom, and
%	that have a value at least 20% higher than that at bottom. In turn, the bottom is computed 
%	based on the cuthyp variable.
auto = ones(nseg,1).*NaN;                               
% Initializes the column vector auto, which will include the first step of automatic sleep scoring 
%	based only on local properties (theta/delta ratio, emg rms) of each epoch. The wake-sleep
%	state is scored as 1 (wakefulness), 2 (non-rapid-eye-movement sleep, NREMS), or 3 (rapid-
%	eye-movement sleep, REMS). Values of 10 and 20 indicate epochs with values of
%	theta/delta ratio and emg rms intermediate between wakefulness and NREMS or between
%	NREMS and REMS, respectively.
auto(MSCO.S(:,2) > MSCO.emgcut(2)) = 1;                      
% wakefulness is scored when emg rms is above the boundary region defined by emgcut
auto(MSCO.S(:,2) < MSCO.emgcut(1) & MSCO.S(:,1) < 0.75) = 2;      
% NREMS is scored when emg rms is below the boundary region defined by emgcut and the
%	theta/delta ratio is lower than 0.75. The latter condition indicates that
%	electroencephalographic spectral power in the delta frequency range is prevalent compared
%	with that in the theta frequency range.
auto(MSCO.S(:,2) < MSCO.emgcut(1) & MSCO.S(:,1) > 1.25) = 3;      
% REMS is scored when emg rms is below the boundary region defined by emgcut and the
%	theta/delta ratio is higher than 1.25. The latter condition indicates that
% 	electroencephalographic spectral power in the theta frequency range is prevalent compared
%	with that in the delta frequency range.
auto(MSCO.S(:,2) <= MSCO.emgcut(2) & MSCO.S(:,2) >= MSCO.emgcut(1)) = 10; 
% An epoch is considered as intermediate between wakefulness and NREMS (code 10) when its
%	emg rms value is included in the boundary region defined by emgcut.
auto(MSCO.S(:,2) < MSCO.emgcut(1) & MSCO.S(:,1) >= 0.75 & MSCO.S(:,1) <= 1.25) = 20; 
% An epoch is considered as intermediate between NREMS and REMS (code 20) when the emg rms
%	is below the boundary region defined by emgcut and theta/delta ratio is between 0.75 and
%	1.25. The latter condition indicates that electroencephalographic spectral power in the delta
%	and theta frequency ranges are nearly equivalent.
auto2 = auto;                                           
% Initializes the column vector auto2, which includes the second and final step of automatic sleep
%	scoring. This step is based on local properties of each epoch, as elaborated in the auto
%	vector, but also takes into account information on the automatic scoring of adjacent epochs. 
%	Codes for wakefulness (1), NREMS (2) and REMS (3) are the same as in the vector auto. The
%	codes 10 and 20 used in the vector auto are replaced by either codes 1-3 or code 4, which
%	indicates an indeterminate state.
for n = 1:length(auto2)-3
	if auto2(n) ~= 3, continue
	else clc
	disp( ['phase one: progress ', mat2str(round (n/length(auto2) * 100)) , ' %'])
        	if ismember(auto2(n+1),[2 10 20]) && auto2(n+2)==3
			auto2(n+1) = 3;
% Single epochs scored as either NREMS or indeterminate state (codes 10 or 20) are re-scored as
%	REMS in case they are preceded and followed by at least one REMS epoch 
        	elseif sum (ismember(auto2(n+1:n+2), [2 10 20])) == 2 && sum (auto2(n+1:n+2) == 2) <2 && auto2(n+3) == 3
			auto2(n+1:n+2) = [3;3];
% Couple of epochs, which are both scored as indeterminate states (codes 10 or 20), or which are 
%	scored as one epoch in indeterminate state and one epoch in NREMS, are re-scored as 2
%	epochs of REMS in case they are preceded and followed by at least one REMS epoch.
        		end
    	end
end
for n = 1:length(auto2)-3
	if auto2(n) ~= 2, continue
    	else clc
		disp(['phase two: progress ',mat2str(round(n/length(auto2)*100)),' %'])
        		if ismember(auto2(n+1),[3 10 20]) && auto2(n+2)==2
			auto2(n+1) = 2;
% Single epochs scored as either REMS or indeterminate state (codes 10 or 20) are re-scored as
% 	NREMS in case they are preceded and followed by at least one NREMS epoch 
        		elseif auto2(n+1)==20 && auto2(n+2)==3 && auto2(n+3)==3
			auto2(n+1) = 4;
% Epochs scored as indeterminate state (code 20) are confirmed as such (code 4) in case they are
%	preceded by at least one NREMS epoch and followed by at least 2 REMS epochs
        		elseif sum (ismember (auto2(n+1:n+2), [3 10 20] )) == 2 && sum(auto2(n+1:n+2)==3)<2 && auto2(n+3)==2
			auto2(n+1:n+2) = [2;2];
% Couple of epochs, which are both scored as indeterminate states (codes 10 or 20), or which are
%	scored as one epoch in indeterminate state and one epoch in REMS, are re-scored as 
% 	epochs of NREMS in case they are preceded and followed by at least one epoch of NREMS.
        		end
    	end
end
auto2((auto2==20)) = 2; auto2((auto2==10)) = 4;
% Epochs that after the previous substitutions remain scored as indeterminate state are re-scored as
%	NREMS, if their previous code was 20, or are confirmed as indeterminate state (code 4), if
%	their previous code was 10.
for n = 1:length(auto2)-3
	if auto2(n) ~= 1, continue
    	else clc
		disp(['phase three: progress ',mat2str(round(n/length(auto2)*100)),' %'])
        		if (ismember(auto2(n+1),[2 3 4]) && auto2(n+2)==1)
			auto2(n+1) = 1;
% Epochs scored as NREMS, REMS, or indeterminate state (code 4) are re-scored as wakefulness in
%	case they are preceded by at least one epoch of wakefulness.
        		elseif (auto2(n+1)==3 && auto2(n+2)==2)
			auto2(n+1) = 4;
% Epochs scored as REMS are re-scored as indetermined state (code 4) in case they are preceded by
%	at least one epoch of wakefulness and followed by at least one epoch of NREMS
        		elseif sum(ismember(auto2(n+1:n+2),[2 3 4]))==2 && auto2(n+3)==1
			auto2(n+1:n+2) = [4;4];
% Couples of epochs, which are scored as either NREMS, REMS, or indeterminate states (code 4) in
%	any combination, are re-scored as indeterminate state (code 4) in case they are preceded and
%	followed by at least one epoch of wakefulness.
        		end
    	end
end

posrem = find(auto2==3);
% The vector posrem indicates the position of the cells in the vector auto2, which are scored as
%	REMS (code 3)
if not(isempty(posrem))
	for n = 1:length(posrem)
        		if (auto2(max(1,posrem(n)-1)) ~= 3) && (auto2(max(1,posrem(n)-1)) == auto2(min(nseg,posrem(n)+1)))
            		auto2(posrem(n)) = auto2(max(1,posrem(n)-1));
        		end
    	end
end
% If an isolated epoch is scored as REMS, and if its adjacent epochs are scored in the same state,
%	and if such a state is different from REMS, that epoch is re-scored with the same state as its
%	adjacent epochs.
auto((auto==10 | auto==20)) = 4;
% Indetermined states in the auto vector (codes 10 and 20) are attributed the code 4, for the sake of
%	consistency and comparison with the auto2 vector
MSCO.S(:,3) = auto; MSCO.S(:,4) = auto2;
% Copies the vectors auto and auto2 in the third and fourth columns, respectively, of the matrix S

znew = MSCO.Z; znew(znew==0)=NaN;
% The matrix znew is a copy of the matrix Z, but with null values substituted by NaNs for the
%	purpose of graphical representation.
figure(1)
set(gcf,'Renderer','zbuffer')
surf(naxis2,maxis,log10(znew)), shading interp, colormap ('jet')
xlabel('emgrms normalized'), ylabel('eeg theta vs. delta ratio'), zlabel('fraction of epochs')
% Figure 1 plots the matrix Z as a 3D graph

figure(2)
powtot= zeros(3,80).*NaN;
powtot(1,1:80)= mean(MSCO.POW(MSCO.S(:,4)==1,1:80));
powtot(2,1:80)= mean(MSCO.POW(MSCO.S(:,4)==2,1:80));
powtot(3,1:80)= mean(MSCO.POW(MSCO.S(:,4)==3,1:80));
%% powtot is a matrix built with the same logic of MSCO.POW. The power spectral density values
%% reported in the 81 columns of powtot correspond to the average values calculated during 
%% wakefulness (first row), %%NREMS (second row) or REMS (third row).

asX = [0.25:0.25:20];
plot(asX,powtot(1,:),'r',asX,powtot(2,:),'g',asX,powtot(3,:),'b')
xlabel('Frequency (Hz)'); ylabel('EEG Power Spectral Density (units^2/Hz)');
%% Figure 2 plots the mean power spectral density (units^2/Hz) during wakefulness (red line),
%% NREMS (green line) and REMS (blue line) in the frequency range 0.25-20 Hz. 

figure(3)
plot(naxis2,MSCO.densi,'rs',naxis2(cutlow),MSCO.densi(cutlow),'ks',naxis2(cuthigh),MSCO.densi(cuthigh),'ks')
xlabel('emgrms normalized'), ylabel('fraction of epochs')
% Figure 3 plots the matrix densi as a 2D graph, and marks the limits of the boundary zone between
%	the higher and lower modes of the densi variable. If the variable densi has 3 modes, users
%	should check on figure 3 whether the boundary zone is correctly positioned between the
%	highest mode and the intermediate mode. If this is not the case, the whole scoring routine
%	should be repeated after setting the value of the variable cuthyp to any position of the densi 
%	vector that is intermediate between the positions of the highest and the intermediate modes.

posrem = find(MSCO.S(:,4)==3);
% The vector posrem indicates the position of the rows in the matrix S, which are scored as REMS
%	(coded 3 in column 4, corresponding to the final step of automatic scoring).
qremtd = length(find(MSCO.S(posrem,1) > 2.5))/length(posrem);
% The variable qremtd is the fraction of REMS epochs with a theta/delta ratio lower than 2.5
if qremtd < 0.1
	disp('WARNING: less than 10% of REMS epochs have a theta/delta EEG power ratio 	higher than 2.5');
    	disp('this indicates that the quality of the raw EEG signal is insufficient for automatic sleep scoring')
    	disp('Manual sleep scoring is thus recommended for this record')
end
% If qremtd is lower than 0.1, users are warned that the automatic sleep scoring may be inadequate
%	because of problems with quality of the raw electroencephalogram 

