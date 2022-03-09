classdef FSKModulationUnit < bherber.thesis.modulation_units.ModulationUnit
    %FSK_MODULATION_UNIT extending the ModulationUnit class
    %   Simulate Frequency-Shift Keying modulation and demodulation for a 
    %   generic system.
    
    methods
        function obj = FSKModulationUnit(signal, data, carrier, params, ...
                distance, square_one, square_zero, freq_channel)
         obj = obj@bherber.thesis.modulation_units.ModulationUnit(...
             signal, data, carrier, params, distance, square_one, square_zero, freq_channel);
        end 

        function f_modulation = modulate(this)
            % MODULATE_BY_FSK a given signal.
            %   Given carrier and data signals, modulate the carrier wave
            %   according to frequency-shift keying as a backscatter tag
            %   would through a given channel.

            % Frequency modulation
            f_modulation = zeros(size(this.data));
            all_ones = this.square_one .* this.carrier * ...
                bherber.thesis.Tag.vanatta_gain(this.params.num_elements, ...
                this.params.Fc, this.params.Fc + this.freq_channel.f1);
            all_zeros = this.square_zero .* this.carrier * ...
                bherber.thesis.Tag.vanatta_gain(this.params.num_elements, ...
                this.params.Fc, this.params.Fc + this.freq_channel.f0);

            for idx = 1:length(f_modulation)
                if this.data(idx) == 1
                    f_modulation(idx) = all_ones(idx);
                else
                    f_modulation(idx) = all_zeros(idx);
                end
            end
        end

        function res_bits = demodulate(this)
        %FSK_DEMODULATE a signal
            %   Demodulate a signal modulated by frequency-shift keying with
            %   frequencies f1 and f0 corresponding to 'on' & 'off'
            %   respectively.

            % 1. Freq mix
            all_ones = this.square_one .* this.carrier;
            all_zeros = this.square_zero .* this.carrier;
            mixed_one = this.signal .* all_ones;
            mixed_zero = this.signal .* all_zeros;

            % 2. Filter
            filtered_one = lowpass(mixed_one, this.freq_channel.f1 * 2, this.params.Fs);
            filtered_zero = lowpass(mixed_zero, this.freq_channel.f1 * 2, this.params.Fs);

            % 3. Mean
            correlated_one = mean(filtered_one);
            correlated_zero = mean(filtered_zero);

            % 4. Combine Streams
            combined_streams = correlated_one - correlated_zero;
%             fprintf("%0.2f %0.2fj\n", real(combined_streams), imag(combined_streams));
%             hold on
%             scatter(gca, real(combined_streams), imag(combined_streams), "filled"); 

            % 4. Decide
            lambda = 0;
            if real(combined_streams) > lambda
                res_bits = 1;
            else
                res_bits = 0;
            end

        end
    end
end

