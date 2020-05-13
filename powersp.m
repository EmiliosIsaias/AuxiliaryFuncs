function [P,features,w_axis] = powersp(signal,fs,freq_band)
%% HELP
% This function returns the power spectrum in a n*m*2 hypermatrix, where N
% is the number of channels, M the signal length and 2 is due to the
% magnitud and the phase of the spectrum. This function also
% returns the mean power of each channel and stores it into a vector of
% length N. Finally, the function returns the frequency axis of the
% spectrum into w_axis. It can also be specified the specific frequency
% band for which the power is desired.
%% SYNTHAXIS
%
%     [Power_hypermatrix,Mean_power_vector,Frequency_axis] = ...
%               powersp(Matrix_of_signals,sampling_frequency,freq_band)

%% FREQUENCY BAND REVIEW
if sum(size(freq_band)) ~= 3
    error('The frequency band is not well structured. Example: [7,13]')
elseif freq_band(1) > freq_band(2)
    warning(['Shifting the frequency band limits order: [',...
        num2str(freq_band(2)),', ',num2str(freq_band(1)),']'])
    aux = freq_band(1);
    freq_band(1) = freq_band(2);
    freq_band(2) = aux;
end
[channels, length_] = size(signal);
% w_axis = 0:fs/lenght_:(fs/2)-(fs/lenght_);
%% ZERO PADDING AND WINDOW CREATION
Z = (2^ceil(log2(length_)))-length_;
L = ceil((length_+Z)/2)+1;
P = zeros(channels,ceil((length_+Z)/2)+1,2);
features = zeros(channels,16);
%Wh = parzenwin(length_);
w_axis = 0:fs/(length_+Z):fs/2;
idx = find(w_axis>=freq_band(1) & w_axis<=freq_band(2));
idx = [idx(1) idx(end)];
%% FEATURES COMPUTATION
for i=1:channels
    display(['Processing channel ' num2str(i)])
    %     s_spect = fftshift(fft([Wh'.*signal(i,:) zeros(1,Z)]));
    %     L = length(s_spect);
    %     L=1;
    % Power of the signal in frequency domain.
    %     s_magn = ((20 * log10(abs(s_spect((L/2):L)))).^2)/L^2;
    %     P(i,:,1) = s_magn;
    %     % Phase [rad]
    %     meanpwr(i) = sum(P(i,:,1));%/L;
    %     P(i,:,2) = unwrap(angle(s_spect((L/2):L)));
    %     P(i,:,1) = s_magn;
    %     P(i,:,2) = s_phase;
    
    spectre=fftshift(fft([signal(i,:) zeros(1,Z)]));
    L2=length(spectre);
    P(i,:,1)=abs(spectre(L-1:L2));
    P(i,:,2)=unwrap(angle(spectre(L-1:L2)));
    features(i,1) = mean(P(i,idx,1).^2);
    features(i,2) = median(P(i,idx,1).^2);
    features(i,3) = mode(P(i,idx,1).^2);
    features(i,4) = min(P(i,idx,1).^2);
    features(i,5) = max(P(i,idx,1).^2);
    features(i,6) = std(P(i,idx,1).^2);
    features(i,7) = var(P(i,idx,1).^2);
    features(i,8) = rms(P(i,idx,1).^2);
    features(i,9) = mean(signal(i,:),2);
    features(i,10) = median(signal(i,:),2);
    features(i,11) = mode(signal(i,:),2);
    features(i,12) = min(signal(i,:));
    features(i,13) = max(signal(i,:));
    features(i,14) = std(signal(i,:));
    features(i,15) = var(signal(i,:));
    features(:,16) = rms(signal(i,:));
%     meanpwr(i,8)=corrcoef(P(i,idx,1).^2);
    
end
end