% function res_bits = fsk_demodulate(signal, carrier, time, f1, f0, params)
% %FSK_DEMODULATE a signal
% %   Demodulate a signal modulated by frequency-shift keying
% 
% % 1. Freq mix
% all_ones = square((f1) * 2 * pi * time) .* carrier;
% all_zeros = square((f0) * 2 * pi * time) .* carrier;
% mixed_one = signal .* all_ones;
% mixed_zero = signal .* all_zeros;
% 
% % 2. Integrate and sharpen
% correlated_one = trapz(time, mixed_one);
% correlated_zero = trapz(time, mixed_zero);
% 
% % 3. Combine Streams
% combined_streams = abs(correlated_one) - abs(correlated_zero);
% 
% % 4. Decide
% lambda = 0;
% if combined_streams > lambda
%     res_bits = 1;
% else
%     res_bits = 0;
% end
% 
% end