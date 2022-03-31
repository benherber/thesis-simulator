classdef PacketSimulator < handle
    %PACKETSIMULATOR object of backscatter tags using system-level TDM.
    %   Simulator of system of backscatter tags communcating with a
    %   basestation employing packet-level 

    properties (GetAccess = public, SetAccess = private)
        tags (1, :) bherber.thesis.PacketTimeDivisionOrchestrator
        tag_locs (3, :) double
        custom_bitstream (1, :) double = []
        tag_mode (1, 1) bherber.thesis.TagType 
        slot_size (1, 1) double
        num_slots (1, 1) double
        tag_preambles (:, 16) {mustBeFinite}
        curr_frame int64
        total_frames int64
        tag_delays (:, 1)
        params bherber.thesis.SimulationConstants
    end

% ----------------------------------------------------------------------- %

    methods
        function this = PacketSimulator(num_packets, tag_locs, tag_mode, ...
                slot_size, num_slots, sim_params, options)
            %SIMULATOR constructor
            %   Initialize and start simulator.

            arguments
                num_packets (1, 1) double
                tag_locs (3, :) double
                tag_mode (1, 1) bherber.thesis.TagType 
                slot_size (1, 1) double
                num_slots (1, 1) double
                sim_params (1, 1) bherber.thesis.SimulationConstants
                options.bitstream (1, :) double = [];
                options.snr_db (1, 1) double = NaN;
                options.complex_noise (1, 1) logical = false;
            end

            this.params = sim_params;
            this.slot_size = slot_size;
            this.num_slots = num_slots;
            this.tag_mode = tag_mode;
            this.tag_locs = tag_locs;
            this.custom_bitstream = options.bitstream;

            % Init Tag ID's
            preamble_len = 16;
            this.tag_preambles = zeros(tag_sz, preamble_len);
            for idx = (1:tag_sz)
                tag_id = dec2bin(randi([0, 2 ^ preamble_len - 1], 1), preamble_len);
                this.tag_preambles(idx, :) = double(tag_id - '0');
            end

            % Init Tag Delay Vector
            this.tag_delays = NaN(tag_sz, 1);

            % Init Simulation Time
            this.curr_frame = int64(0);
            this.total_frames = num_packets;
            
        end

% ----------------------------------------------------------------------- %

        function generate_orchestrated_tags(this)

            tag_sz = size(this.tag_locs, 2);
            tag_pos = this.tag_locs;

            % Init Data
            bits = zeros(tag_sz, (this.params.num_symbs) * log2(this.params.m_ary_modulation));
            tags_objs = [];
            for idx = 1:tag_sz
                % Get random data signal
                if isempty(options.bitstream)
                    bits(idx, :) = randi([0, 1], 1, ...
                        (this.params.num_symbs) * log2(this.params.m_ary_modulation));
                else
                    bits(idx, :) = this.custom_bitstream;
                end

                bits(idx, 1:16) = this.tag_preambles(idx, :);

                switch this.tag_mode
                    case bherber.thesis.TagType.OOK
                        curr_tag = bherber.thesis.tags.OOKTag(...
                            tag_pos(1, idx), tag_pos(2, idx), tag_pos(3, idx), ...
                            tag_modes(idx), ...
                            bits(idx, :), ...
                            this.params, ...
                            snr_db=options.snr_db, ...
                            complex_noise=options.complex_noise);
                    case bherber.thesis.TagType.FSK_LO
                        curr_tag = bherber.thesis.tags.FSKTag(...
                            tag_pos(1, idx), tag_pos(2, idx), tag_pos(3, idx), ...
                            tag_modes(idx), ...
                            bits(idx, :), ...
                            this.params, ...
                            this.params.freq_channels(1, :), ...
                            snr_db=options.snr_db, ...
                            complex_noise=options.complex_noise);
                    case bherber.thesis.TagType.FSK_HI
                        curr_tag = bherber.thesis.tags.FSKTag(...
                            tag_pos(1, idx), tag_pos(2, idx), tag_pos(3, idx), ...
                            tag_modes(idx), ...
                            bits(idx, :), ...
                            this.params, ...
                            this.params.freq_channels(2, :), ...
                            snr_db=options.snr_db, ...
                            complex_noise=options.complex_noise);
%                     case bherber.thesis.TagType.FREQ_HOP
%                         curr_tag = bherber.thesis.tags.FreqHopTag(...
%                             tag_pos(1, idx), tag_pos(2, idx), tag_pos(3, idx), ...
%                             tag_modes(idx), ...
%                             bits(idx, :), ...
%                             this.params, ...
%                             this.params.freq_channels, ...
%                             hopsets(idx, :), ...
%                             snr_db=options.snr_db, ...
%                             complex_noise=options.complex_noise);
                    otherwise
                        error("Unsupported type: %s", this.tag_mode);
                end

%                 rand_slot = randi([1, this.num_slots], 1, 1);
                curr_orchestrator = bherber.thesis.PacketTimeDivisionOrchestrator(...
                    this.slot_size, this.num_slots, 1, curr_tag);

                if isempty(tags_objs)
                    tags_objs = repmat(curr_orchestrator, 1, tag_sz);
                else
                    tags_objs(idx) = curr_orchestrator; 
                end
            end
            this.tags = tags_objs;
        end

% ----------------------------------------------------------------------- %

        function res = step(this)
            %STEP through one frame of the simulation
            %   Step one packet-sized frame in time through the simulation
            %   for the entire system of tags.

            if this.curr_frame > this.total_frames
                error("Tried to simulate too many frames (%d of %d)", ...
                    this.curr_frame, this.total_frames);
            end

            this.generate_orchestrated_tags();

            res = zeros(int64(this.params.simstep_sz), this.params.sim_sym_ratio, ...
                this.params.num_symbs, this.num_slots);

            % FXIME: Make not so nested
            for idx = 1:this.num_slots
                for jdx = 1:this.params.num_symbs
                    for kdx = 1:this.params.sim_sym_ratio
                        for tag_idx = 1:length(this.tags)
                            res(:, kdx, jdx, idx) = res(:, kdx, jdx, idx) + ...
                                this.tags(tag_idx).step();
                        end
                    end
                end
            end

            res = res(:);
            this.curr_frame = this.curr_frame + 1;
        end

% ----------------------------------------------------------------------- %

        function res = auto_align(this)
            %AUTO_ALIGN each received signal.
            %   Find delays of each received tag signal based off of given
            %   preambles.

            calcdelay = @(distance) ...
                ceil((2 * distance / physconst("Lightspeed")) * this.params.Fs);

            for idx = 1:length(this.tags)
                this.tag_delays(idx) = calcdelay(this.tags(idx).tag.distance);
            end

            res = isempty(this.tag_delays(isnan(this.tag_delays)));
        end

% ----------------------------------------------------------------------- %

        function res = deframeify(this, idx, frame)
            arguments
                this (1, 1) bherber.thesis.PacketSimulator
                idx (1, 1) double
                frame (1, :)
            end

            res = this.tags(idx).tag.deframeify(frame);
        end
    end
end