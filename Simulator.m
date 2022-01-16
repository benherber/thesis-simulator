% ------------------------ %
% Passive Tag Simulator    %
% Name   : Simulator.m     %
% Author : Benjamin Herber %
% Date   : Fall 2021       %
% ------------------------ %

classdef Simulator < handle
    %SIMULATOR object of backscatter tags.
    %   Simulator of system of backscatter tags communcating with a
    %   basestation.

    properties (GetAccess = public, SetAccess = private)
        tags (:, 1) Tag
        tag_preambles (:, 8) {mustBeFinite}
        curr_step int64
        total_steps int64
        time (1, :) {mustBeFinite, mustBeReal}
        carrier (1, :) {mustBeFinite}
        channel function_handle
        tag_delays (:, 1)
        params SimulationConstants
    end

    methods
        function this = Simulator(tags, tag_modes, channel, sim_params)
            %SIMULATOR constructor
            %   Initialize and start simulator.

            this.params = sim_params;

            % Get time
            this.time = (0:(1 / this.params.Fs):(this.params.total_time - (1 / this.params.Fs)));

            % Carrier signal
            this.carrier = complex(this.params.amplitude * sin(complex(2 * pi * this.params.Fc * this.time)));

            % Init Channel
            this.channel = channel;

            % Init tags
            tag_sz = size(tags);
            tag_sz = tag_sz(2);
            for idx = 1:tag_sz
                % Get random data signal
                bits = randi([0, 1], this.params.num_symbs, 1);
                data = repelem(bits, this.params.symb_sz).';

                this.tags = [this.tags, Tag(tags(1, idx), tags(2, idx), tags(3, idx), ...
                    tag_modes(idx), this.time, this.carrier, data, this.channel, this.params)];
            end

            % Init Tag ID's
            this.tag_preambles = zeros(tag_sz, 8);
            for idx = (1:tag_sz)
                tag_id = dec2bin(idx - 1, 8);
                this.tag_preambles(idx, :) = double(tag_id - '0');
            end

            % Init Tag Delay Vector
            this.tag_delays = NaN(tag_sz, 1);

            % Init Simulation Time
            this.curr_step = int64(0);
            this.total_steps = int64(this.params.num_symbs * this.params.sim_sym_ratio);

            gcp;
        end

        function res = step(this)
            %STEP through one frame of the simulation
            %   Step one 1us frame in time through the simulation for the
            %   entire system of tags.

            if this.curr_step > this.total_steps
                error("ATTEMPTED TO SIMULATE TOO MANY STEPS");
            end

            res = zeros(1, this.params.simstep_sz);

            for idx = 1:numel(this.tags)
                curr = this.tags(idx);
                res = res + curr.step();
            end

            this.curr_step = this.curr_step + 1;
        end

        function res = auto_align(this, signal)
            if (this.curr_step * this.params.sim_sym_ratio) > 8
                res = false;
                return;
            end

            curr_samps = floor(this.curr_step / this.params.sim_sym_ratio) * ...
                this.params.symb_sz;
            t = this.time(1:curr_samps);
            carr = this.carrier(1:curr_samps);
            opts = this.params;
            symb_size = this.params.symb_sz;
            tag_preams = this.tag_preambles;
            thetags = this.tags;
            delays = this.tag_delays;
            parfor idx = 1:length(this.tags)
                preamble = repelem(tag_preams(idx, :), symb_size).';
                if thetags{idx}.mode == TagType.OOK
                    modded_preamble = Tag.modulate_by_ook(carr, preamble, opts);
                else
                    if thetags{idx}.mode == TagType.FSK_LO
                        f1 = opts.fsk_channel0.f1;
                        f0 = opts.params.fsk_channel0.f0;

                    elseif thetags{idx}.mode == TagType.FSK_HI
                        f1 = opts.fsk_channel1.f1;
                        f0 = opts.fsk_channel1.f0;
                    else
                        error("BAD MODE");
                    end

                    modded_preamble = Tag.modulate_by_fsk(t, carr, opts, f1, f0);

                end

                delay = finddelay(modded_preamble, signal);
                if delay ~= 0
                    delays(idx) = delay;
                end

            end
            this.tag_delays = delays;

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