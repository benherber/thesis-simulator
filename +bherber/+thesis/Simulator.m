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
        tags (:, 1)
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
        function this = Simulator(tags, tag_modes, channel, sim_params, options)
            %SIMULATOR constructor
            %   Initialize and start simulator.

            arguments
                tags
                tag_modes
                channel
                sim_params
                options.bitstream = [];
                options.snr_db = NaN;
                options.complex_noise = false;
            end

            this.params = sim_params;

            % Get time
            this.time = []; %(0:(1 / this.params.Fs):(this.params.total_time - (1 / this.params.Fs)));

            % Init tags
            tag_sz = size(tags, 2);

            % Init Channel
            this.channel = channel;

            % Carrier signal
            this.carrier = []; %this.params.amplitude * cos(2 * pi * this.params.Fc * this.time);
%             noisy_carrier = channel(this.carrier);

            % Init Tag ID's
            this.tag_preambles = zeros(tag_sz, 8);
            for idx = (1:tag_sz)
                tag_id = dec2bin(randi([0, 2^8 - 1], 1), 8);
                this.tag_preambles(idx, :) = double(tag_id - '0');
            end

            % Set up Hop-Patterns
            hopsets = bherber.thesis.HopsetsGenerator.generate(tag_sz, this.params.num_channels, 1);

            % Init Data
            this.bits = zeros(tag_sz, (this.params.num_symbs) * log2(this.params.m_ary_modulation));
            tags_objs = [];
            for idx = 1:tag_sz
                % Get random data signal
                if isempty(options.bitstream)
                    this.bits(idx, :) = randi([0, 1], 1, (this.params.num_symbs) * log2(this.params.m_ary_modulation));
                else
                    this.bits(idx, :) = options.bitstream;
                end

                switch tag_modes(idx)
                    case bherber.thesis.TagType.OOK
                        curr_tag = bherber.thesis.tags.OOKTag(...
                            tags(1, idx), tags(2, idx), tags(3, idx), tag_modes(idx), ...
                            this.bits(idx, :), this.params, snr_db=options.snr_db, complex_noise=options.complex_noise);
                    case bherber.thesis.TagType.FSK_LO
                        curr_tag = bherber.thesis.tags.FSKTag(...
                            tags(1, idx), tags(2, idx), tags(3, idx), tag_modes(idx), ...
                            this.bits(idx, :), this.params, this.params.freq_channels(1), ...
                            snr_db=options.snr_db, complex_noise=options.complex_noise);
                    case bherber.thesis.TagType.FSK_HI
                        curr_tag = bherber.thesis.tags.FSKTag(...
                            tags(1, idx), tags(2, idx), tags(3, idx), tag_modes(idx), ...
                            this.bits(idx, :), this.params, this.params.freq_channels(2), ...
                            snr_db=options.snr_db, complex_noise=options.complex_noise);
                    case bherber.thesis.TagType.FREQ_HOP
                        curr_tag = bherber.thesis.tags.FreqHopTag(...
                            tags(1, idx), tags(2, idx), tags(3, idx), tag_modes(idx), ...
                            this.bits(idx, :), this.params, this.params.freq_channels, hopsets(idx, :));
                end

                if isempty(tags_objs)
                    tags_objs = repmat(curr_tag, 1, tag_sz);
                else
                    tags_objs(idx) = curr_tag; 
                end
                %bherber.thesis.Tag(tags(1, idx), tags(2, idx), tags(3, idx), ...
%                     tag_modes(idx), this.time, this.carrier, bits, this.channel, this.params);
            end
            this.tags = tags_objs;

            % Init Tag Delay Vector
            this.tag_delays = NaN(tag_sz, 1);

            % Init Simulation Time
            this.curr_step = int64(0);
            this.total_steps = int64((this.params.num_symbs + 1) * this.params.sim_sym_ratio);
            
        end

        function res = step(this)
            %STEP through one frame of the simulation
            %   Step one 1us frame in time through the simulation for the
            %   entire system of tags.
            %   NOTE: Requires noise to be applied to signal post-steps.

            if this.curr_step > this.total_steps
                error("ATTEMPTED TO SIMULATE TOO MANY STEPS");
            end

            res = zeros(1, int64(this.params.simstep_sz));

            for idx = 1:length(this.tags)
                res = res + this.tags(idx).step();
            end


            this.curr_step = this.curr_step + 1;
        end

        function res = auto_align(this)
            %AUTO_ALIGN each received signal.
            %   Find delays of each received tag signal based off of given
            %   preambles.

            calcdelay = @(distance) ceil((2 * distance / physconst("Lightspeed")) * this.params.Fs);

            for idx = 1:length(this.tags)
%                 fprintf("idx:%d dist:%f, %f\n", idx, this.tags(idx).distance, calcdelay(this.tags(idx).distance));
                this.tag_delays(idx) = calcdelay(this.tags(idx).distance);
            end

            res = isempty(this.tag_delays(isnan(this.tag_delays)));
        end
    end
end