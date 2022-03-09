classdef OOKModulationUnit < bherber.thesis.modulation_units.ModulationUnit
    %OOK_MODULATION_UNIT extending the ModulationUnit class
    %   Simulate On-Off Keying modulation and demodulation for a generic
    %   system.
    
    methods
        function obj = OOKModulationUnit(signal, data, carrier, params, ...
                distance, square_one, square_zero, freq_channel)
         obj = obj@bherber.thesis.modulation_units.ModulationUnit(...
             signal, data, carrier, params, distance, square_one, square_zero, freq_channel);
        end 

        function ook_modulation = modulate(this)
            % MODULATE a given signal by on-off keying.
            %   Given carrier and data signals, modulate the carrier wave
            %   according to on-off keying as a backscatter tag would
            %   through a given channel.

            ook_modulation = this.carrier .* this.data .* ...
                bherber.thesis.Tag.vanatta_gain(params.num_elements, params.Fc, params.Fc);
        end

        function res_bits = demodulate(this)
            %OOK_DEMODULATE a signal
            %   Demodulate a signal modulated by on-off keying.

            % 1. Freq mix
            mixed_ook = (this.signal .* this.carrier);

            % 2. Filter
            filtered = lowpass(mixed_ook, this.params.symb_freq * 2, this.params.Fs);

            % 3. Mean Value
            correlated_ook = mean(filtered);
            
            % 4. Decide
            expected_amplitude = this.params.amplitude * ...
                ((this.params.wavelen / (4 * pi * this.distance)) ^ 2) * ...
                 bherber.thesis.Tag.vanatta_gain(this.params.num_elements, this.params.Fc, this.params.Fc);
            lambda = (expected_amplitude * this.params.amplitude) / 4;
            if real(correlated_ook) > lambda
                res_bits = 1;
            else
                res_bits = 0;
            end
        end
    end
end

