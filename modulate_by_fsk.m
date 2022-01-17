% %{
%  -------------------------- 
%  FSK Modulation             
%  Name   : modulate_by_fsk.m 
%  Author : Benjamin Herber   
%  Date   : Fall 2021         
%  --------------------------
% %}
% 
% function [ reconstructed_data, modulated_signal ] = modulate_by_fsk(t, carrier, data, channel)
% % MODULATE_BY_FSK Modulate a given signal by frequency-shift keying.
% %   Given carrier and data signals, modulate the carrier wave according to
% %   frequency-shift keying as a backscatter tag would through a given channel.
% 
% %% Load Contants
% 
% simulation_constants
% 
% %% Frequency modulation
% 
% % 1. Add channel noise headed to tag
% noisy_carrier = channel(carrier) * 4;
% 
% % 2. Modulate
% f_modulation = zeros(size(data));
% sq_one = square((FB1) * 2 * pi * t);
% sq_zero = square((FB0) * 2 * pi * t);
% all_ones = sq_one .* noisy_carrier;
% all_zeros = sq_zero .* noisy_carrier;
% for idx = 1:length(f_modulation)
%     if data(idx) == 1
%         f_modulation(idx) = all_ones(idx);
%     else
%         f_modulation(idx) = all_zeros(idx);
%     end
% end
% 
% % 3. Add channel noise headed back to basestation
% noisy_f_modulation = f_modulation;
% 
% %% Frequency Demodulation
% 
% % 1. Freq mix
% mixed_one = noisy_f_modulation .* all_ones;
% mixed_zero = noisy_f_modulation .* all_zeros;
% 
% % 2. Integrate and sharpen
% filtered_one = zeros(size(mixed_one));
% for idx = (0:(NUM_SAMPLES - 1))
%     left_bound = int64(idx * LENGTH_OF_SYMBOL * Fs + 1);
%     right_bound = int64((idx + 1) * LENGTH_OF_SYMBOL * Fs);
%     filtered_one(left_bound:right_bound) = ...
%         trapz(t(left_bound:right_bound), ...
%         mixed_one(left_bound:right_bound));
% end
% norm_filtered_one = (filtered_one - min(filtered_one)) / (max(filtered_one) - min(filtered_one));
% for idx = (1:numel(filtered_one))
%     if norm_filtered_one(idx) > 0.5
%         norm_filtered_one(idx) = 1;
%     else
%         norm_filtered_one(idx) = 0;
%     end
% end
% 
% filtered_zero = zeros(size(mixed_zero));
% for idx = (0:(NUM_SAMPLES - 1))
%     left_bound = int64(idx * LENGTH_OF_SYMBOL * Fs + 1);
%     right_bound = int64((idx + 1) * LENGTH_OF_SYMBOL * Fs);
%     filtered_zero(left_bound:right_bound) = ...
%         trapz(t(left_bound:right_bound), ...
%         mixed_zero(left_bound:right_bound));
% end
% norm_filtered_zero = (filtered_zero - min(filtered_zero)) / (max(filtered_zero) - min(filtered_zero));
% for idx = (1:numel(filtered_zero))
%     if norm_filtered_zero(idx) > 0.5
%         norm_filtered_zero(idx) = 1;
%     else
%         norm_filtered_zero(idx) = 0;
%     end
% end
% 
% % 3. Combine Streams
% combined_streams = norm_filtered_one - norm_filtered_zero;
% 
% % 4. Decide
% lambda = 0;
% reconstructed_data = zeros(size(t));
% for idx = 1:length(reconstructed_data)
%     if combined_streams(idx) > lambda
%         reconstructed_data(idx) = 1;
%     else
%         reconstructed_data(idx) = 0;
%     end
% end
% 
% % 5. Reconstruct bitstream
% resbits = downsample(reconstructed_data, (LENGTH_OF_SYMBOL * Fs), 2).';
% 
% %% Return
% 
% reconstructed_data = 0;
% modulated_signal = noisy_f_modulation;
% 
% end
% 
