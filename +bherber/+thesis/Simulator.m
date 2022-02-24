% ------------------------ %
% Passive Tag Simulator    %
% Name   : Simulator.m     %
% Author : Benjamin Herber %
% Date   : Fall 2021       %
% ------------------------ %

classdef Simulator < handle
    %SIMULATOR object of backscatter tags.
    %   Simulator of system of backscatter tags communcating with a
    %   basestation. Effective range is ~300m.

    properties (GetAccess = public, SetAccess = private)
        tags (:, 1)bherber.thesis.Tag
        tag_preambles (:, 8) {mustBeFinite}
        curr_step int64
        total_steps int64
        time (1, :) {mustBeFinite, mustBeReal}
        carrier (1, :) {mustBeFinite}
        channel function_handle
        tag_delays (:, 1)
        params bherber.thesis.SimulationConstants
        bits {mustBeReal}
    end

    methods
        function this = Simulator(tags, tag_modes, channel, sim_params)
            %SIMULATOR constructor
            %   Initialize and start simulator.

            this.params = sim_params;

            % Get time
            this.time = (0:(1 / this.params.Fs):(this.params.total_time - (1 / this.params.Fs)));

            % Init tags
            tag_sz = size(tags, 2);

            % Init Channel
            this.channel = channel;

            % Carrier signal
            this.carrier = complex(this.params.amplitude * sin(complex(2 * pi * this.params.Fc * this.time)));
            noisy_carrier = channel(this.carrier);

            % Init Tag ID's
            this.tag_preambles = zeros(tag_sz, 8);
            for idx = (1:tag_sz)
                tag_id = dec2bin(randi([0, 2^8 - 1], 1), 8);
                this.tag_preambles(idx, :) = double(tag_id - '0');
            end

            % Init Data
            this.bits = zeros(tag_sz, this.params.num_symbs);
            tags_objs = repmat(bherber.thesis.Tag(0, 0, 0, bherber.thesis.TagType.OOK, ...
                1, 1, 1, this.channel, this.params), 1, tag_sz);
            for idx = 1:tag_sz
                % Get random data signal
                bits = [this.tag_preambles(idx, :), randi([0, 1], 1, this.params.num_symbs - 8)];
                this.bits(idx, :) = bits;
                data = repelem(bits, this.params.symb_sz);
                data = [data, zeros(1, length(this.time) - length(data))].';

                tags_objs(idx) = bherber.thesis.Tag(tags(1, idx), tags(2, idx), tags(3, idx), ...
                    tag_modes(idx), this.time, noisy_carrier, data, this.channel, this.params);
            end
            this.tags = tags_objs;

            % Init Tag Delay Vector
            this.tag_delays = NaN(tag_sz, 1);

            % Init Simulation Time
            this.curr_step = int64(0);
            this.total_steps = int64((this.params.num_symbs + 1) * this.params.sim_sym_ratio);

            gcp;
        end

        function res = step(this)
            %STEP through one frame of the simulation
            %   Step one 1us frame in time through the simulation for the
            %   entire system of tags.
            %   NOTE: Requires noise to be applied to signal post-steps.

            if this.curr_step > this.total_steps
                error("ATTEMPTED TO SIMULATE TOO MANY STEPS");
            end

            res = zeros(length(this.tags), this.params.simstep_sz);

            for idx = 1:numel(this.tags)
                curr = this.tags(idx);
                res(idx, :) = curr.step();
            end


            this.curr_step = this.curr_step + 1;
        end

        function res = auto_align(this)
            %AUTO_ALIGN each received signal.
            %   Find delays of each received tag signal based off of given
            %   preambles.

            if (this.curr_step / this.params.sim_sym_ratio) < 8
                res = false;
                return;
            end

            calcdelay = @(distance) floor((distance / physconst("Lightspeed")) * this.params.Fs);

            for idx = 1:length(this.tags)
                fprintf("idx:%d dist:%f, %f\n", idx, this.tags(idx).distance, calcdelay(this.tags(idx).distance));
                this.tag_delays(idx) = calcdelay(this.tags(idx).distance);
            end

            res = isempty(this.tag_delays(isnan(this.tag_delays)));
        end
    end

    methods (Static, Access = public)
        function gray = dec2gray(num)
            %DEC2GRAY conversion.
            %   Convert a decimal number to a gray encoding. Return array.

            arguments
                num {mustBeFinite, mustBePositive}
            end

            % Convert to binary
            bits = dec2bin(num);

            % Convert to Gray
            gray = zeros(size(bits), "logical");
            gray(1) = bits(1);
            for idx = (2:numel(bits))
                gray(idx) = xor(logical(str2double(bits(idx - 1))), ...
                    logical(str2double(bits(idx))));
            end
        end

        function dec = gray2dec(gray)
            %GRAY2DEC conversion.
            %   Convert a decimal number to a gray encoding. Return array.

            arguments
                gray {mustBeFinite}
            end

            % Convert to Gray
            bits = zeros(size(gray), "like", gray);
            bits(1) = gray(1);
            for idx = (2:length(gray))
                bits(idx) = xor(bits(idx - 1), gray(idx));
            end

            % Convert to Decimal
            dec = bin2dec(num2str(bits));
        end
    end
end