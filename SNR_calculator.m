function snr = SNR_calculator(data, data_noise, noise, p)
    snr= zeros(1, length(p)+1);
    snr(1) = 10 * log10(sum(data_noise.^2)/ sum(noise.^2));
    disp(['SNR (noised signal): ' num2str(snr(1)) ' dB']);
    
    for i = 1:length(p)
        fileName = strcat('guit1_filtered_', num2str(p(i)), '.wav');

        [z, ~] = audioread(fileName);
        
        noise_temp = z - data;
        snr(1, i+1) = 10 * log10(sum(z.^2)/ sum(noise_temp.^2));
        disp(['SNR (' fileName '): ' num2str(snr(1, i+1)) ' dB']);
    end
end

%     fileName = strcat('guit1_filtered_', num2str(p(i)), '.wav');
%     fileName1 = strcat('guit1_filtered_F_', num2str(p(i)), '.wav');
% 
%     audiowrite(fileName , z, sample_rate);
%     audiowrite(fileName1, z, sample_rate);
   