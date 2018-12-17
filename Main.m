%% clear data
clc;
clear;
close all;
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
    [FFT_amp_ecg_0, freq_ecg_0] = calculate_FFT(ecg, ecg_rijen, measure_freq_ecg);
    [FFT_amp_ecg2_0, freq_ecg2_0] = calculate_FFT(ecg2, ecg2_rijen, measure_freq_ecg2);
    
    %find PL noise frequency ecg and ecg2
    thresh_0 = 44; %Hz 50 are 60 Hz are sources for the PLN
    PLN_freq_degrees_ecg_0 = freq_ecg_0(find(freq_ecg_0 > thresh_0,1));
    PLN_freq_degrees_ecg2_0 = freq_ecg2_0(find(freq_ecg2_0 > thresh_0,1)); 
    
    %calculate notch and apply filters Transfer Functions:
    a = 0.9;
    [Notch_ecg_1, FFT_amp_ecg_1, freq_ecg_1, z1_ecg, z2_ecg] = calculate_notch_with_conj(ecg_tijdskolom, PLN_freq_degrees_ecg_0, a, measure_freq_ecg, ecg_rijen);    
    [Notch_ecg_clean, FFT_amp_ecg_clean, freq_ecg_clean] = clean_notch(100, freq_ecg_1, Notch_ecg_1, a, measure_freq_ecg, ecg_rijen);
 
    [Notch2_ecg, FFT_amp_ecg2_1, freq_ecg2_1, z1_ecg2, z2_ecg2] = calculate_notch_with_conj(ecg2_tijdskolom, PLN_freq_degrees_ecg2_0, a, measure_freq_ecg2, ecg2_rijen);
    
%% Task 4
    %impulse respons zie notch function
    %0
    radials = 0;
    z1 = cos(radials) + 1j * sin(radials);
    z2 = conj(z1);
    impulse_response(z1,z2);
    
    %fs/4
    radials = pi/2;
    z1 = cos(radials) + 1j * sin(radials);
    z2 = conj(z1);
    impulse_response(z1,z2);
    
    %fs/2
    radials = pi;
    z1 = cos(radials) + 1j * sin(radials);
    z2 = conj(z1);
    impulse_response(z1,z2);
%% personal functions:
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

function [Notch, amp, freq_array_out, z1, z2] = calculate_notch_with_conj(input, degrees, a, freq, rows)
    radials = degrees *2*pi / freq;
    z1 = cos(radials) + 1j * sin(radials);
    z2 = conj(z1);
    impulse_response(z1,z2);
    Notch = [input(:,1), zeros(rows,1)];
    Notch(1,2) = input(1,2);
    Notch(2,2) = input(2,2) - input(1,2)*(z1+z2) + Notch(1,2) *(a*z1+a*z2);
    Notch(3:end, 2) =  input(3:end,2) - input(2:end-1,2)*(z1+z2) + Notch(2:end-1,2) *(a*z1+a*z2)+ ( input(1:end-2,2) - Notch(1:end-2,2) * a^2 )* z1*z2;
    
    [amp, freq_array_out] = calculate_FFT(Notch(:,2), rows, freq);
end

function [Notch_ecg, FFT_amp_ecg, freq_ecg] = clean_notch(thresh, ecg_freq, Notch_ecg, a, freq, rows)
    freq_degrees = ecg_freq(find(ecg_freq > thresh,1));
    while( freq/2 > thresh)
        [Notch_ecg, FFT_amp_ecg, freq_ecg] = calculate_notch_with_conj(Notch_ecg, freq_degrees, a, freq, rows);
        freq_degrees = freq_ecg(find(freq_ecg > thresh,1));    
        thresh = freq_degrees - 50;
        [~,n] = size( freq_degrees);
        if n == 0
            thresh = freq;
        end
    end
end

function impulse_response(z1 ,z2)
    sys = filt([1 -z1-z2 z1*z2],[1 -0.9*z1-0.9*z2 0.81*z1*z2]);
    figure
    impulse(sys);
    title('impulse response');
    figure
    [Amp,freq_norm] = freqz([1 -z1-z2 z1*z2], [1 -0.9*z1-0.9*z2 0.81*z1*z2],2000);
    plot(freq_norm/pi, 20*log10(abs(Amp)));
    title('frequency response');
    figure
    bode(sys)
end