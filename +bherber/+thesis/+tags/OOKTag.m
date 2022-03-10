% ----------------------------- %
% On-Off Keying Backscatter Tag %
% Name   : OOKTag.m             %
% Author : Benjamin Herber      %
% Date   : Spring 2021          %
% ----------------------------- %

classdef OOKTag < bherber.thesis.tags.Tag
    %OOK_TAG model
    %   Backscatter tag model of a passive Van Atta configuration using
    %   On-Off Keying to communicate

    %% Public Instance Methods

    methods (Access = public)

        function this = OOKTag(x, y, z, mode, bits, sim_params)
            %OOK_TAG constructor
            %   Create a tag a specific distance from a basestation in 3D
            %   space.
            
            arguments
                x double {mustBeScalarOrEmpty} % x-position of tag relative to basestation
                y double {mustBeScalarOrEmpty} % y-position of tag relative to basestation
                z double {mustBeScalarOrEmpty} % z-position of tag relative to basestation
                mode bherber.thesis.TagType % Modulation mode
                bits (1, :) {mustBeNonmissing, mustBeInRange(bits, 0, 1)} % bitstream
                sim_params bherber.thesis.SimulationConstants % Simulation parameters
            end

            this@bherber.thesis.tags.Tag(x, y, z, mode, bits, sim_params);

        end

% ----------------------------------------------------------------------- %
    
        function res = step(this)
            %STEP one simulation step using OOK modulation.
    
            if ~this.is_active
                    res = zeros(1, this.params.simstep_sz);
                    return
            end
    
            details = this.step_data();
            res = details.carrier_slice .* details.data .* ...
                    this.vanatta_gain(this.params.num_elements, this.params.Fc, this.params.Fc);
        end
    
% ----------------------------------------------------------------------- %
    
    function res_bits = demodulate(this, symb_num, signal)
            arguments
                this bherber.thesis.tags.OOKTag
                symb_num double {mustBePositive}
                signal (1, :)
            end

            time_slice = (0:(1 / this.params.Fs):(this.params.symb_sz * (1 / this.params.Fs) - (1 / this.params.Fs))) ...
                + (symb_num * double(this.params.symb_sz) * (1 / this.params.Fs));

            carrier_slice = this.params.amplitude * cos(2 * pi * this.params.Fc * time_slice);

            % 1. Freq mix
            mixed_ook = (signal .* carrier_slice);

            % 2. Filter
            filtered = lowpass(mixed_ook, this.params.symb_freq * 2, this.params.Fs);

            % 3. Mean Value
            correlated_ook = mean(filtered);
            
            % 4. Decide
            expected_amplitude = this.params.amplitude * ...
                ((this.params.wavelen / (4 * pi * this.distance)) ^ 2) * ...
                 this.vanatta_gain(this.params.num_elements, this.params.Fc, this.params.Fc);
            lambda = (expected_amplitude * this.params.amplitude) / 4;
            if real(correlated_ook) > lambda
                res_bits = 1;
            else
                res_bits = 0;
            end
        end
    end
end
