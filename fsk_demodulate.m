function res_bits = fsk_demodulate(signal, carrier, time, f1, f0, params)
%FSK_DEMODULATE a signal
%   Demodulate a signal modulated by frequency-shift keying

% 1. Freq mix
all_ones = square((f1) * 2 * pi * time) .* carrier;
all_zeros = square((f0) * 2 * pi * time) .* carrier;
mixed_one = signal .* all_ones;
mixed_zero = signal .* all_zeros;

% 2. Integrate and sharpen
correlated_one = trapz(time, mixed_one);
correlated_zero = trapz(time, mixed_zero);

% 3. Combine Streams
combined_streams = abs(correlated_one) - abs(correlated_zero);

% 4. Decide
lambda = 0;
if combined_streams > lambda
    res_bits = 1;
else
    res_bits = 0;
end

% % 1. Freq mix
% sq_one = square((f1) * 2 * pi * time);
% sq_zero = square((f0) * 2 * pi * time);
% all_ones = sq_one .* carrier;
% all_zeros = sq_zero .* carrier;
% mixed_one = signal.' .* all_ones;
% mixed_zero = signal.' .* all_zeros;
% 
% % 2. Integrate and sharpen
% filtered_one = zeros(size(mixed_one));
% for idx = (0:(params.num_symbs - 1))
%     left_bound = int64(idx * params.symb_sz + 1);
%     right_bound = int64((idx + 1) * params.symb_sz);
%     filtered_one(left_bound:right_bound) = ...
%         trapz(time(left_bound:right_bound), ...
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
% for idx = (0:(params.num_symbs - 1))
%     left_bound = int64(idx * params.symb_sz + 1);
%     right_bound = int64((idx + 1) * params.symb_sz);
%     filtered_zero(left_bound:right_bound) = ...
%         trapz(time(left_bound:right_bound), ...
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
% reconstructed_data = zeros(size(time));
% for idx = 1:length(reconstructed_data)
%     if combined_streams(idx) > lambda
%         reconstructed_data(idx) = 1;
%     else
%         reconstructed_data(idx) = 0;
%     end
% end
% 
% % 5. Reconstruct bitstream
% res_bits = downsample(reconstructed_data, params.symb_sz, 2).';

end