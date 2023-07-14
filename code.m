clear; clc; close all;

%% ---- QUESTION A ----
disp('---- QUESTION A ----')
[data, sample_rate] = audioread('guit1.wav');

figure('Name', 'A & B knowing raw data');
subplot(3,3,1); plot(data); title('Input signal');

% Greating white gaussian noise with 0.01 deviation
deviation = 0.01; 
noise = deviation * randn(size(data));  
n = length(noise) ;
clearvars deviation

subplot(3,3,4); plot(noise); title('Noise');

% Adding noise to my data
data_noise = data + noise;
audiowrite('guit1_noise.wav', data_noise, sample_rate);

subplot(3,3,7); plot(data_noise); title('Signal with noise');

% Wiener Filter --> FIR
p = [10, 20, 30];
plot_order = [2 5 8];
for i = 1:length(p)
    
    z = wiener_filter(p(i), data, data_noise, n);
    
    subplot(3,3, plot_order(i)); plot(z); 
    title(['Wiener Filter order ' num2str(p(i))]);
        
    fileName  = strcat('guit1_filtered_'  , num2str(p(i)), '.wav');
    fileName1 = strcat('guit1_filtered_A_', num2str(p(i)), '.wav');
    audiowrite(fileName , z, sample_rate);
    audiowrite(fileName1, z, sample_rate);
end

% Calculating SNR of noised and denoised signal
snr_A = SNR_calculator(data, data_noise, noise, p);

disp(' ')

%% ---- QUESTION B ----
disp('---- QUESTION B ----')

frameSize = 256;
overlap = 0.5;
windows        = frame_wind(data, frameSize, overlap);
windows_noised = frame_wind(data_noise, frameSize, overlap);

% Appling Wiener filter to each window
filteredWindows = zeros(size(windows_noised));
for i = 1:length(p)
    fileName  = strcat('guit1_filtered_', num2str(p(i)), '.wav');
    fileName1 = strcat('guit1_filtered_B_', num2str(p(i)), '.wav');
    
    for j = 1:size(windows_noised, 2)
        window       = windows(:, j);
        window_noise = windows_noised(:, j);
                
        filteredWindow = wiener_filter(p(i), window, window_noise, frameSize);
        filteredWindows(:, j) = filteredWindow;
    end
    
    % Convert the filtered windows back to a single denoised signal
    z = frame_recon(filteredWindows, overlap);
    
    subplot(3,3, i*3); plot(z); title(['Windows + Wiener order ' num2str(p(i))]);
    
    audiowrite(fileName , z, sample_rate);
    audiowrite(fileName1, z, sample_rate);
end
% sound(z, sample_rate);

% Calculating SNR of noised and denoised signal
temp = frame_wind(data, frameSize, overlap);
temp = frame_recon(temp, overlap);
snr_B = SNR_calculator(temp', data_noise, noise, p);

clearvars i j z fileName1 temp
clearvars filteredWindow filteredWindows filterLength 
clearvars window windows window_noise windows_noise
disp(' ')

%% ---- QUESTION C ----
disp('---- QUESTION C ----')

correlation = xcorr(noise);
% Plot the autocorrelation function
lag = -(length(noise)-1):(length(noise)-1);
figure('Name', 'C : White Gaussian Noise')
subplot(2,1,1);plot(lag, correlation); title('Correlation of Noise');
xlabel('Lag');  ylabel('Correlation');

subplot(2,1,2); histogram(correlation); title('Histogram');

disp('For white noise, the autocorrelation should be close to zero at all lags, except for a peak at lag 0 due to the noise''s variance.');
disp('According to figure 2, I have white noise.');

clearvars correlation lag correlation
disp(' ')

%% ---- QUESTION D ----
disp('---- QUESTION D ----')

figure('Name', 'Data with noise --> Finding an array with jusy noise');
subplot(2,1,1);    plot(data_noise);   % noise in [95440, 102600]
disp('Noised is located in [95440 102600]');
noise_d = data_noise(95440 : 102600 , :);
subplot(2,1,2); plot(noise_d);

figure('Name', 'D & E without knowing raw data');
subplot(3,3,1); plot(data); title('Input signal');
subplot(3,3,4); plot(noise); title('Noise');
subplot(3,3,7); plot(data_noise); title('Signal with noise');

% Finding each correlation
corr_noised = xcorr(data_noise);
corr_noise  = xcorr(noise_d);

% Finding the mean of each vector 
mean_noised = round(length(corr_noised)/2, 0);
mean_noise  = round(length(corr_noise )/2, 0);

% Repeat for raw input in order to evaluate my results
corr_signal = xcorr(data);
mean_signal = round(length(corr_signal)/2, 0);
disp('The mean square error of the index I have found and the true indexes is:');
disp(['Filter''s order'  '   ' 'mse' '      ' 'max_value' '     ' 'min_value'] );

% Getting the most important indexes (order of the filter) and evaluating
indexes = zeros(3, max(p));
mse_d = zeros(1, length(p));
for i = 1:length(p)
    indexes_noised = corr_noised(mean_noised : mean_noised+p(i)-1 );    % rxx
    indexes_noise  = corr_noise (mean_noise  : mean_noise +p(i)-1 );    % rnn
    
    temp = (indexes_noised - indexes_noise)';                           % rxd
    
    % zero padding
    if length(temp) < max(p)
        temp = [temp, zeros(1,max(p)-p(i))];
    end
    
    indexes(i,:) = temp;
    
    % Evaluating my results with Mean Square Error
    indexes_signal = (corr_signal(mean_signal : mean_signal+p(i)-1 ))';
    
%     display(['Wiener Filter order: ', num2str(p(i))]);
%     display('First my indexes, then the true indexes of the raw signal');
%     display(indexes(i,1:p(i)));
%     display(indexes_signal);
    
    mse_d(i) = immse( indexes(i,1:p(i)), indexes_signal  );
    
    display( ['     ' num2str(p(i)) '        ' num2str(mse_d(i)) '     ' num2str(max(indexes_signal)) '     ' num2str(min(indexes_signal)) ] );
    
    % Using the above colloration in Wiener Filter
    Rxx = toeplitz(indexes_noised);
    w = pinv(Rxx, 0.0001)*temp(1:p(i))';
    z = filter(w,1,data_noise);
    
    subplot(3,3, plot_order(i)); plot(z); 
    title(['Wiener Filter order ' num2str(p(i))]);
    
    fileName  = strcat('guit1_filtered_'  , num2str(p(i)), '.wav');
    fileName1 = strcat('guit1_filtered_D_', num2str(p(i)), '.wav');

    audiowrite(fileName , z, sample_rate);
    audiowrite(fileName1, z, sample_rate);
end
disp(' ');

% Calculating SNR
snr_D = SNR_calculator(data, data_noise, noise, p);

clearvars indexes indexes_noised indexes_noise indexes_signal temp i j n 
clearvars mean_noise mean_noised mean_signal noise_d
clearvars corr_noise corr_noised corr_signal 
clearvars rxd rxx Rxx w z telos fileName1
disp(' ')

%% ---- QUESTION Ε ----
disp('---- QUESTION Ε ----')
threshold = 0.0001;

windows        = frame_wind(data, frameSize, overlap);
windows_noised = frame_wind(data_noise, frameSize, overlap);

% Appling Wiener filter to each window
filteredWindows = zeros(size(windows_noised));
for i = 1:length(p)
    fileName = strcat('guit1_filtered_', num2str(p(i)), '.wav');
    fileName1 = strcat('guit1_filtered_E_', num2str(p(i)), '.wav');
    
    window_noise  = windows_noised(:, 1);
    
    for j = 1:size(windows_noised, 2)
        
        if mean(abs(windows_noised(:,j))) < threshold
            window_noise = windows_noised(:,j);
        end
        
        window_noised = windows_noised(:,j);
        
        corralation_noised = xcorr(window_noised);      % rxx
        corralation_noise  = xcorr(window_noise );      % rnn
        
        temp = corralation_noised - corralation_noise;  % rxd
        
        % Appling Winer Filter
        telos = frameSize + p(i) -1 ;
        rxd = temp(frameSize : telos);
        rxx = corralation_noised(frameSize : telos);
        Rxx = toeplitz(rxx);
        w = pinv(Rxx, 0.0001)*rxd;
        filteredWindow = filter(w,1,window_noised);
                
        filteredWindows(:, j) = filteredWindow;
    end
    
    % Convert the filtered windows back to a single denoised signal
    z = frame_recon(filteredWindows, overlap);
    
    subplot(3,3, i*3); plot(z); 
    title(['Windows + Wiener order ' num2str(p(i))]);
    
    audiowrite(fileName , z, sample_rate);
    audiowrite(fileName1, z, sample_rate);
end
% sound(z, sample_rate);

% Calculating SNR of noised and denoised signal
temp = frame_wind(data, frameSize, overlap);
temp = frame_recon(temp, overlap);
snr_E = SNR_calculator(temp', data_noise, noise, p);

clearvars telos rxd rxx Rxx w i j z fileName1 threshold
clearvars corralation_noised corralation_noise
clearvars filteredWindow filteredWindows  
clearvars window windows window_noised windows_noised windows_noise
clearvars temp overlap frameSize
disp(' ')

%% Deleteing temporary files
for i = 1:length(p)
    fileName = strcat('guit1_filtered_'  , num2str(p(i)), '.wav');
    delete(fileName);
end

%% ---- QUESTION F ----
disp('---- QUESTION F ----')
values = [2 10 15];
snr_F = zeros(1, length(p));

figure('Name', 'F: Prediction')
    
frameSize = 256;
overlap = 0.5;
windows = frame_wind(data, frameSize, overlap);
data2   = frame_recon(windows, overlap);

% Appling Wiener filter to each window
filteredWindows = zeros(size(windows));
for i = 1:length(p)
    for v = 1:length(values)
    % Wiener Filter
    for j = 1:size(windows, 2)
        window = windows(:, j);
        
        telos = frameSize+p(i)-1;
        rss = xcorr(data_noise);      
        rxx = rss(frameSize   : telos  );
        rss = rss(frameSize+values(v) : telos+values(v));
        Rxx = toeplitz(rxx);
        w = pinv(Rxx,0.0001)*rss;
        filteredWindow = filter(w,1,window);
        
        filteredWindows(:, j) = filteredWindow;
    end
    
    % Convert the filtered windows back to a single denoised signal
    z = frame_recon(filteredWindows, overlap);
    
    subplot(3,3,(v*3-2)+i-1); plot(z); 
    hold on; plot(data); plot(z,'r'); hold off;
    title(['Windows + Wiener order ' num2str(p(i)) ' ' num2str(values(v)) ' elements']);
    
    % Calculating SNR of noised and denoised signal
    noise_temp = data2 - z;
    snr_F(1, i) = 10 * log10(sum(z.^2)/ sum(noise_temp.^2));
    disp(['SNR (Wiener ' num2str(p(i)) ' rank) with ' num2str(values(v)) ' elements: ' num2str(snr_F(1, i)) ' dB']);
    
    end
end

clearvars i j z fileName overlap frameSize 
clearvars telos rss rxx Rxx w data2 noise_temp
clearvars window windows window_noise filteredWindow filteredWindows

clearvars sample_rate plot_order p
disp(' ')

disp('----------------- Done ----------------');
disp('--------- Bourcha Ioanna 58019 --------');
disp('---------------- Summer <3 ------------');