%% clear data
clc;
clear;

%% load data
load("ecg.mat");
load("ecg2.mat");

measure_freq_ecg = 1000; % Hz 
measure_freq_ecg2 = 204.73; % Hz 
%% get amount ecg data points
[ecg_rijen, ~] = size(ecg);
[ecg2_rijen, ~] = size(ecg2);

%% Task 1:
    % tijdswaarden toevoegen aan ecg data
    %voeg kolom met eenen toe
    ecg_tijdskolom = [ones(ecg_rijen, 1), ecg];
    ecg2_tijdskolom = [ones(ecg2_rijen, 1), ecg2];
    %verander enen naar tijdswaarden
    for i = 1: ecg_rijen
        ecg_tijdskolom(i,1) = i * 1/measure_freq_ecg;
    end
    for i = 1: ecg2_rijen
        ecg2_tijdskolom(i,1) = i * 1/measure_freq_ecg2;
    end
    
    % plot ifv tijd
    figure
    plot(ecg_tijdskolom(:,1), ecg_tijdskolom(:,2));
    title("ECG data infv de tijd in seconden");
    xlabel("Tijd(s)");
    ylabel("ECG amplitude data");
    figure
    plot(ecg2_tijdskolom(:,1), ecg2_tijdskolom(:,2));
    title("ECG2 data infv de tijd in seconden");
    xlabel("Tijd(s)");
    ylabel("ECG amplitude data");

%% Task 2:
% powerline spectrum:
    %ecg
    FFT_ecg = fft(ecg);
    P2_ecg = abs(FFT_ecg / ecg_rijen);
    P1_ecg = P2_ecg(1:ecg_rijen/2+1)/length(ecg_rijen);
    P1_ecg(2:end-1) = 2 * P1_ecg(2:end-1);
    f_ecg = measure_freq_ecg *(0:(ecg_rijen/2))/ecg_rijen;
    figure
    hold on
    plot(f_ecg, P1_ecg)
    [FFT_amp_ecg, freq_ecg] = findpeaks(P1_ecg, f_ecg, 'MinPeakDistance',9 , 'MinPeakHeight', 0.02);
    plot(freq_ecg, FFT_amp_ecg, 'x');
    xlabel('Frequency [Hz]');
    ylabel('Amplitude');
    legend('FFT');
    hold off 
    
    %ecg2
    FFT_ecg2 = fft(ecg2);
    P2_ecg2 = abs(FFT_ecg2 / ecg2_rijen);
    P1_ecg2 = P2_ecg2(1:ecg2_rijen/2+1)/length(ecg2_rijen);
    P1_ecg2(2:end-1) = 2 * P1_ecg2(2:end-1);
    f_ecg2 = measure_freq_ecg2 *(0:(ecg2_rijen/2))/ecg2_rijen;
    figure
    hold on
    plot( f_ecg2,P1_ecg2)
    [FFT_amp_ecg2, freq_ecg2] = findpeaks(P1_ecg2, f_ecg2, 'MinPeakDistance',5, 'MinPeakHeight', 0.02);
    plot(freq_ecg2, FFT_amp_ecg2, 'x');
    xlabel('Frequency [Hz]');
    ylabel('Amplitude');
    legend('FFT');
    hold off
    
    %find pwl noise frequency ecg and ecg2
    thresh = 44; %Hz
    PLN_freq_degrees_ecg = freq_ecg(find(freq_ecg > thresh,1));
    PLN_freq_degrees_ecg2 = freq_ecg2(find(freq_ecg2 > thresh,1)); 
    
    %calculate notch filters:
    a = 0.9;
    %ecg:
    PLN_freq_radials_ecg = (PLN_freq_degrees_ecg / 180) * pi;
    PLN_z1_ecg = cos(PLN_freq_radials_ecg) + 1j * sin(PLN_freq_radials_ecg);
    PLN_z2_ecg = conj(PLN_z1_ecg);
    
    Numerator_ecg = {[1 -PLN_z1_ecg] ; [1 -PLN_z2_ecg]};
    Denominator_ecg = {[1 -a * PLN_z1_ecg] ; [1 - a * PLN_z2_ecg]};
    
    H_notch_ecg = filt(Numerator_ecg, Denominator_ecg, measure_freq_ecg);
    %ecg2
    PLN_freq_radials_ecg2 = (PLN_freq_degrees_ecg2 / 180) * pi;
    PLN_z1_ecg2 = cos(PLN_freq_radials_ecg2) + 1j * sin(PLN_freq_radials_ecg2);
    PLN_z2_ecg2 = conj(PLN_z1_ecg2);
    
    Numerator_ecg2 = {[1 -PLN_z1_ecg2] ; [1 -PLN_z2_ecg2]};
    Denominator_ecg2 = {[1 -a * PLN_z1_ecg2] ; [1 - a * PLN_z2_ecg2]};
    
    H_notch_ecg2 = filt(Numerator_ecg2, Denominator_ecg2, measure_freq_ecg2);
    
    % filter data
    
    
    
    
    
    
    
    
    