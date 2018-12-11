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
%visualize sample data:
    % tijdswaarden toevoegen aan ecg data
    %voeg kolom met eenen toe
    ecg_tijdskolom = [ones(ecg_rijen, 1), ecg];
    ecg2_tijdskolom = [ones(ecg2_rijen, 1), ecg2];
    %verander enen naar tijdswaarden en plot
    ecg_tijdskolom = convert_to_time_and_plot(ecg_tijdskolom, measure_freq_ecg, ecg_rijen);
    ecg2_tijdskolom = convert_to_time_and_plot(ecg2_tijdskolom, measure_freq_ecg2, ecg2_rijen);
    
    
%% Task 2:
% powerline spectrum:
    %visualization
    [FFT_amp_ecg, freq_ecg] = calculate_FFT(ecg, ecg_rijen, measure_freq_ecg);
    [FFT_amp_ecg2, freq_ecg2] = calculate_FFT(ecg2, ecg2_rijen, measure_freq_ecg2);
    
    %find PL noise frequency ecg and ecg2
    thresh = 44; %Hz 50 are 60 Hz are sources for the PLN
    PLN_freq_degrees_ecg = freq_ecg(find(freq_ecg > thresh,1));
    PLN_freq_degrees_ecg2 = freq_ecg2(find(freq_ecg2 > thresh,1)); 
    
    %calculate notch filters Transfer Functions:
    a = 0.9;
%     H_notch_ecg = calculate_notch_with_conj(PLN_freq_degrees_ecg, a, measure_freq_ecg);    
%     H_notch_ecg2 = calculate_notch_with_conj(PLN_freq_degrees_ecg2, a, measure_freq_ecg2);   


    radials = (PLN_freq_degrees_ecg *2*pi / measure_freq_ecg);
    z1 = cos(radials) + 1j * sin(radials);
    z2 = conj(z1);
    
%     Y/X = (z - z1)*(z-z2) / (z - a z1)*(z - a z2);
%     Y (z - a z1)*(z- a z2) = X (z - z1)*(z-z2) / (z - z1)*(z-z2);
%     Y (z^2 - (a*z1+a*z2)*z + a^2 * z1*z2)= X (z^2 - (z1+z2)*z + z1*z2)
%     Y z^2 - Y *(a*z1+a*z2)*z +Y* a^2 * z1*z2= X z^2 - X*(z1+z2)*z + X* z1*z2
%     Y(z) - Y(z-1) *(a*z1+a*z2) +Y(z-2)* a^2 * z1*z2= X(z) - X(z-1)*(z1+z2) + X(z-2)* z1*z2
%implement:
%Y(z) =  X(z) - X(z-1)*(z1+z2) + Y(z-1) *(a*z1+a*z2)+ X(z-2)* z1*z2 - Y(z-2)* a^2 * z1*z2;
    
    % filter data
    Y = [ecg_tijdskolom(:,1), zeros(ecg_rijen,1)];
    Y(1,2) = ecg_tijdskolom(1,2);
    Y(2,2) = ecg_tijdskolom(2,2) - ecg_tijdskolom(1,2)*(z1+z2) + Y(1,2) *(a*z1+a*z2);
    Y(3:end, 2) =  ecg_tijdskolom(3:end,2) - ecg_tijdskolom(2:end-1,2)*(z1+z2) + Y(2:end-1,2) *(a*z1+a*z2)+ ( ecg_tijdskolom(1:end-2,2) - Y(1:end-2,2) * a^2 )* z1*z2;
    
    [FFT_amp_ecg, freq_ecg] = calculate_FFT(Y(:,2), ecg_rijen, measure_freq_ecg);
    
    
    
    
    
    
    
%% personal function:
function output = convert_to_time_and_plot(input, freq, rows)
    output = input;
    for i = 1: rows
        output(i,1) = i * 1/freq;
    end
    
    figure
    plot(output(:,1), output(:,2));
    title("ECG data infv de tijd in seconden");
    xlabel("Tijd(s)");
    ylabel("ECG amplitude data");
end

function [amp, freq] = calculate_FFT(input, rijen, mfreq)
    FFT = fft(input);
    P2 = abs(FFT / rijen);
    P1 = P2(1:rijen/2+1)/length(rijen);
    P1(2:end-1) = 2 * P1(2:end-1);
    f = mfreq *(0:(rijen/2))/rijen;
    figure
    hold on
    plot(f, P1)
    [amp, freq] = findpeaks(P1, f, 'MinPeakDistance',9 , 'MinPeakHeight', 0.02);
    
    plot(freq, amp, 'x');
    xlabel('Frequency [Hz]');
    ylabel('Amplitude');
    legend('FFT');
    hold off 
end

function Notch = calculate_notch_with_conj(degrees, a, freq)
    radials = (degrees / 180) * pi;
    z1 = cos(radials) + 1j * sin(radials);
    z2 = conj(z1);

    Numerator_ecg = {[1 -z1] ; [1 -z2]};
    Denominator_ecg = {[1 -a * z1] ; [1 - a * z2]};

    Notch = filt(Numerator_ecg, Denominator_ecg, freq);
end