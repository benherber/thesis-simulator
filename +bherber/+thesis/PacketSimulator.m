classdef PacketSimulator < handle
    %PACKETSIMULATOR object of backscatter tags using system-level TDM.
    %   Simulator of system of backscatter tags communcating with a
    %   basestation employing packet-level 

    properties (GetAccess = public, SetAccess = private)
        tags (1, :) bherber.thesis.tags.Tag
%         hopsets (:, :) double = []
        tag_locs (3, :) double
        tag_mode (1, 1) bherber.thesis.TagType 
        num_slots (1, 1) double
        slot_selects (1, :) double
        tag_preambles (:, :) {mustBeFinite}
        curr_frame int64
        total_frames int64
        tag_delays (:, 1)
        options (1, 1) struct
        params bherber.thesis.SimulationConstants
    end

% ----------------------------------------------------------------------- %

    methods
        function this = PacketSimulator(num_packets, tag_locs, tag_mode, ...
                num_slots, sim_params, options)
            %SIMULATOR constructor
            %   Initialize and start simulator.

            arguments
                num_packets (1, 1) double
                tag_locs (3, :) double
                tag_mode (1, 1) bherber.thesis.TagType 
                num_slots (1, 1) double
                sim_params (1, 1) bherber.thesis.SimulationConstants
                options.bitstream (1, :) double = [];
                options.snr_db (1, 1) double = NaN;
                options.complex_noise (1, 1) logical = false;
            end

            this.params = sim_params;
            this.num_slots = num_slots;
            this.tag_mode = tag_mode;
            this.tag_locs = tag_locs;
            this.options = options;

            % Init Tag ID's
            preamble_len = this.params.preamble_len;
            preambles = randperm(2 ^ preamble_len, size(tag_locs, 2)).';
            this.tag_preambles = double(dec2bin(preambles, preamble_len) - '0');
%             for idx = (1:size(tag_locs, 2))
%                 tag_id = dec2bin(randi([0, 2 ^ preamble_len - 1], 1), preamble_len);
%                 this.tag_preambles(idx, :) = double(tag_id - '0');
%             end

            % Init Tag Delay Vector
            this.tag_delays = NaN(size(tag_locs, 2), 1);

            % Init Simulation Time
            this.curr_frame = int64(0);
            this.total_frames = num_packets;
            
%             if tag_mode == bherber.thesis.TagType.FREQ_HOP
%                 this.hopsets = zeros(size(tag_locs, 2), this.params.pattern_len);
%                 chosen_idxs = randperm(this.params.num_channels ^ this.params.pattern_len, size(tag_locs, 2));
%                 for idx = 1:length(chosen_idxs)
%                     vals = dec2base(chosen_idxs(idx) - 1, this.params.num_channels, this.params.pattern_len);
%                     this.hopsets(idx, :) = arrayfun(@(val) str2double(val), vals) + 1;
%                 end
%             end
        end

% ----------------------------------------------------------------------- %

        function generate_tags(this)

            tag_sz = size(this.tag_locs, 2);
            tag_pos = this.tag_locs;

            % Init Data
            bits = zeros(tag_sz, (this.params.num_symbs) * log2(this.params.m_ary_modulation));
            tags_objs = [];
            for idx = 1:tag_sz
                % Get random data signal
                if isempty(this.options.bitstream)
                    bits(idx, :) = randi([0, 1], 1, ...
                        (this.params.num_symbs) * log2(this.params.m_ary_modulation));
                else
                    bits(idx, :) = this.options.bitstream;
                end

                bits(idx, 1:this.params.preamble_len) = this.tag_preambles(idx, :);

                switch this.tag_mode
                    case bherber.thesis.TagType.OOK
                        curr_tag = bherber.thesis.tags.OOKTag(...
                            tag_pos(1, idx), tag_pos(2, idx), tag_pos(3, idx), ...
                            this.tag_mode, ...
                            bits(idx, :), ...
                            this.params, ...
                            snr_db=this.options.snr_db, ...
                            complex_noise=this.options.complex_noise);
                    case bherber.thesis.TagType.FSK_LO
                        curr_tag = bherber.thesis.tags.FSKTag(...
                            tag_pos(1, idx), tag_pos(2, idx), tag_pos(3, idx), ...
                            this.tag_mode, ...
                            bits(idx, :), ...
                            this.params, ...
                            this.params.freq_channels(1, :), ...
                            snr_db=this.options.snr_db, ...
                            complex_noise=this.options.complex_noise);
                    case bherber.thesis.TagType.FSK_HI
                        curr_tag = bherber.thesis.tags.FSKTag(...
                            tag_pos(1, idx), tag_pos(2, idx), tag_pos(3, idx), ...
                            this.tag_mode, ...
                            bits(idx, :), ...
                            this.params, ...
                            this.params.freq_channels(2, :), ...
                            snr_db=this.options.snr_db, ...
                            complex_noise=this.options.complex_noise);
                    case bherber.thesis.TagType.FREQ_HOP
%                         curr_pattern_idx = mod(this.curr_frame + 1, this.params.pattern_len);
%                         if curr_pattern_idx == 0; curr_pattern_idx = this.params.pattern_len; end
                        freq_idx = randi(this.params.num_channels, 1, 1);
                        curr_tag = bherber.thesis.tags.FSKTag(...
                            tag_pos(1, idx), tag_pos(2, idx), tag_pos(3, idx), ...
                            this.tag_mode, ...
                            bits(idx, :), ...
                            this.params, ...
                            this.params.freq_channels(freq_idx, :), ...
                            snr_db=this.options.snr_db, ...
                            complex_noise=this.options.complex_noise);
                    otherwise
                        error("Unsupported type: %s", this.tag_mode);
                end

%                 rand_slot = randi([1, this.num_slots], 1, 1);
%                 curr_orchestrator = bherber.thesis.PacketTimeDivisionOrchestrator(...
%                     this.slot_size, this.num_slots, 1, curr_tag);

                if isempty(tags_objs)
                    tags_objs = repmat(curr_tag, 1, tag_sz);
                else
                    tags_objs(idx) = curr_tag; 
                end

                rand_slot = randi([1, this.num_slots], 1, 1);
                this.slot_selects(idx) = rand_slot;
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

            this.generate_tags();

            res = zeros(1, int32(this.num_slots * (this.params.num_symbs + 1) * ...
                        this.params.sim_sym_ratio * this.params.simstep_sz));
            for idx = 1:length(this.tags)
                res = res + bherber.thesis.orchestrate_tdm( ...
                    this.num_slots, this.slot_selects(idx), this.tags(idx));
            end

            this.curr_frame = this.curr_frame + 1;
        end

% ----------------------------------------------------------------------- %

        function auto_align(this)
            %AUTO_ALIGN each received signal.
            %   Find delays of each received tag signal based off of given
            %   preambles.

            calcdelay = @(tag) ...
                ceil((2 * tag.distance / physconst("Lightspeed")) * tag.params.Fs);

            this.tag_delays = arrayfun(calcdelay, this.tags);
        end

% ----------------------------------------------------------------------- %

        function res = deframeify(this, idx, frame)
            arguments
                this (1, 1) bherber.thesis.PacketSimulator
                idx (1, 1) double
                frame (1, :)
            end
            this.auto_align()

            slot_size = int32((this.params.num_symbs + 1) * ...
                this.params.sim_sym_ratio * this.params.simstep_sz);
            slots = reshape(frame, slot_size, this.num_slots);
            res = delayseq(slots(:, this.slot_selects(idx)), -1 * this.tag_delays(idx)).';
        end

% ----------------------------------------------------------------------- %

        function clash = clash_in_next_frame(this)
            arguments
                this bherber.thesis.PacketSimulator
            end
    
            frame = this.step();
            clash = 0;
            for idx = 1:size(this.tag_locs, 2)
                symbols = reshape(this.deframeify(idx, frame), [], this.params.num_symbs + 1);
                resbits = zeros(this.params.num_symbs, 1);
                for jdx = 1:this.params.num_symbs
                    resbits(jdx) = this.tags(idx).demodulate(jdx, symbols(:, jdx).');
                end

                clash = clash + logical(biterr(this.tags(idx).bits(1:this.params.num_symbs), resbits));
            end
        end

    end
end