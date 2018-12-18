%% clear workspace and commandline and close all figures
clc;
clear;
close all;

%% load provided data
load("ecg.mat");
load("ecg2.mat");
measure_freq_ecg = 1000; % Hz 
measure_freq_ecg2 = 204.73; % Hz 

%% get amount ecg data points
[ecg_rows, ~] = size(ecg);
[ecg2_rows, ~] = size(ecg2);

%% Task 1:    
%visualize sample data:
    % add a column with time values to ecg data
    ecg_timecolumn = [ones(ecg_rows, 1), ecg];
    ecg_timecolumn = convert_to_time_and_plot(ecg_timecolumn, measure_freq_ecg, ecg_rows);
		title("ECG.mat data infv de tijd in seconden");
    
	ecg2_timecolumn = [ones(ecg2_rows, 1), ecg2];
    ecg2_timecolumn = convert_to_time_and_plot(ecg2_timecolumn, measure_freq_ecg2, ecg2_rows);
		title("ECG2.mat data infv de tijd in seconden");
    
    pause;
%% Task 2:
% powerline spectrum:
    %visualization
    [FFT_amp_ecg_0, freq_ecg_0] = calculate_FFT(ecg, ecg_rows, measure_freq_ecg);
    title('FFT van ecg.mat')
    [FFT_amp_ecg2_0, freq_ecg2_0] = calculate_FFT(ecg2, ecg2_rows, measure_freq_ecg2);
    title('FFT van ecg2.mat')
    
    %find PL noise frequency ecg and ecg2
    thresh_0 = 44; %Hz 50 are 60 Hz are sources for the PLN
    PLN_freq_degrees_ecg_0 = freq_ecg_0(find(freq_ecg_0 > thresh_0,1));
    PLN_freq_degrees_ecg2_0 = freq_ecg2_0(find(freq_ecg2_0 > thresh_0,1)); 

    %calculate notch and apply filters Transfer Functions:
    a = 0.9;
    [Notch_ecg_1, z1_ecg_1, z2_ecg_1] = calculate_notch_with_conj(ecg_timecolumn, PLN_freq_degrees_ecg_0, a, measure_freq_ecg, ecg_rows);    
    [~, freq_ecg_1] = calculate_FFT(Notch_ecg_1(:,2), ecg_rows, measure_freq_ecg);
    title('ecg.mat dataset na eenmalig te filteren met de Notch filter.')
    thresh_1 = 100;
    [Notch_ecg_2, z1_ecg_2, z2_ecg_2] = clean_notch(thresh_1, freq_ecg_1, Notch_ecg_1, a, measure_freq_ecg, ecg_rows);

    
    [Notch_ecg2_1, z1_ecg2, z2_ecg2] = calculate_notch_with_conj(ecg2_timecolumn, PLN_freq_degrees_ecg2_0, a, measure_freq_ecg2, ecg2_rows);
    [~, ~] = calculate_FFT(Notch_ecg2_1(:,2), ecg2_rows, measure_freq_ecg2);
    title('ecg2.mat dataset na eenmalig te filteren met de Notch filter.')    
    pause;
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
    
%% Task 6
%removing baseline wander 
    HP_filter = designfilt('highpassfir',...       % Response type
       'StopbandFrequency',0.3, ...     % Frequency constraints
       'PassbandFrequency',0.5, ...
       'StopbandAttenuation',55, ...    % Magnitude constraints
       'PassbandRipple',2, ...
       'DesignMethod','kaiserwin', ...  % Design method
       'ScalePassband',false, ...       % Design method options
       'SampleRate',measure_freq_ecg2);
    ecg2_HP = filtfilt(HP_filter, ecg2);
    ecg2_HP_tijdskolom = [ones(ecg2_rows, 1), ecg2_HP];
    %verander enen naar tijdswaarden en plot
    ecg2_HP_tijdskolom = convert_to_time_and_plot(ecg2_HP_tijdskolom, measure_freq_ecg2, ecg2_rows);

% removing high frequency noise:
    LP_filter = designfilt('lowpassfir', ...        % Response type
       'PassbandFrequency',40, ...     % Frequency constraints
       'StopbandFrequency',40.2, ...
       'StopbandAttenuation',55, ...    % Magnitude constraints
       'PassbandRipple',2, ...
       'DesignMethod','kaiserwin', ...  % Design method
       'ScalePassband',false, ...
       'SampleRate',measure_freq_ecg2);
    ecg2_LP = filtfilt(LP_filter, ecg2);
    ecg2_LP_tijdskolom = [ones(ecg2_rows, 1), ecg2_LP];
    %verander enen naar tijdswaarden en plot
    ecg2_LP_tijdskolom = convert_to_time_and_plot(ecg2_LP_tijdskolom, measure_freq_ecg2, ecg2_rows);
    
    % removing both at the same time:
    BP_filter = designfilt('bandpassfir', ...       % Response type
       'StopbandFrequency1',0.3, ...    % Frequency constraints
       'PassbandFrequency1',0.5, ...
       'PassbandFrequency2',40, ...
       'StopbandFrequency2',40.2, ...
       'StopbandAttenuation1',55, ...    % Magnitude constraints
       'PassbandRipple',2, ...
       'StopbandAttenuation2',55, ...    % Magnitude constraints
       'DesignMethod','kaiserwin', ...  % Design method
       'ScalePassband',false, ...
       'SampleRate',measure_freq_ecg2); 
    ecg2_BP = filtfilt(BP_filter, ecg2);
    ecg2_with_bp = [ones(ecg2_rows, 1), ecg2_BP];
    %verander enen naar tijdswaarden en plot
    ecg2_with_bp = convert_to_time_and_plot(ecg2_with_bp, measure_freq_ecg2, ecg2_rows);
    
%% Task 7
    LP_filter =designfilt('lowpassiir', ...        % Response type
       'PassbandFrequency',55, ...     % Frequency constraints
       'StopbandFrequency',55.2, ...
       'PassbandRipple',2, ...          % Magnitude constraints
       'StopbandAttenuation',55, ...
       'DesignMethod','ellip', ...      % Design method
       'MatchExactly','both', ...   % Design method options
       'SampleRate',measure_freq_ecg);
    ecg_LP = filtfilt(LP_filter, ecg);
    ecg_LP_t = [ones(length(ecg_LP), 1), ecg_LP];
    %verander enen naar tijdswaarden en plot
    ecg_LP_t = convert_to_time_and_plot(ecg_LP_t, measure_freq_ecg, length(ecg_LP));
    
    
    ecg_downsampled = downsample(ecg_LP, 10);
    ecg_downsampled_tijdskolom = [ones(length(ecg_downsampled), 1), ecg_downsampled];
    %verander enen naar tijdswaarden en plot
    ecg_downsampled_tijdskolom = convert_to_time_and_plot(ecg_downsampled_tijdskolom, 6, length(ecg_downsampled));
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

function [Notch, z1, z2] = calculate_notch_with_conj(input, degrees, a, freq, rows)
    radials = degrees *2*pi / freq;
    z1 = cos(radials) + 1j * sin(radials);
    z2 = conj(z1);
    Notch = [input(:,1), zeros(rows,1)];
    Notch(1,2) = input(1,2);
    Notch(2,2) = input(2,2) - input(1,2)*(z1+z2) + Notch(1,2) *(a*z1+a*z2);
    Notch(3:end, 2) =  input(3:end,2) - input(2:end-1,2)*(z1+z2) + Notch(2:end-1,2) *(a*z1+a*z2)+ ( input(1:end-2,2) - Notch(1:end-2,2) * a^2 )* z1*z2;
end

function [Notch_ecg, z1, z2] = clean_notch(thresh, ecg_freq, Notch_ecg, a, freq, rows)
    freq_degrees = ecg_freq(find(ecg_freq > thresh,1));
    while( freq/2 > thresh)
        [Notch_ecg, z1,z2] = calculate_notch_with_conj_clean(Notch_ecg, freq_degrees, a, freq, rows);
        [~,freq_ecg] = calculate_FFT(Notch_ecg(:,2), rows, freq);
        title('filtered more than once with the Notch filter')
        freq_degrees = freq_ecg(find(freq_ecg > thresh,1));    
        thresh = freq_degrees - 50;
        [~,n] = size( freq_degrees);
        if n == 0
            thresh = freq;
        end
    end
end

function [Notch, z1, z2] = calculate_notch_with_conj_clean(input, degrees, a, freq, rows)
    radials = degrees *2*pi / freq;
    z1 = cos(radials) + 1j * sin(radials);
    z2 = conj(z1);
    Notch = [input(:,1), zeros(rows,1)];
    Notch(1,2) = input(1,2);
    Notch(2,2) = input(2,2) - input(1,2)*(z1+z2) + Notch(1,2) *(a*z1+a*z2);
    Notch(3:end, 2) =  input(3:end,2) - input(2:end-1,2)*(z1+z2) + Notch(2:end-1,2) *(a*z1+a*z2)+ ( input(1:end-2,2) - Notch(1:end-2,2) * a^2 )* z1*z2;
   
end

function impulse_response(z1 ,z2)
    Num = [1 -z1-z2 z1*z2];
    Denum = [1 -0.9*z1-0.9*z2 0.81*z1*z2];
    figure
    impz(Num,Denum,30);
    title('impulse response');
    figure
    freqz(Num,Denum);
end