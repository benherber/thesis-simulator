classdef PacketTimeDivisionOrchestrator < handle
    %PACKETTIMEORCHESTRATOR implementation.
    %   An Orchestrator that implements Time Division Multiplexing
    %   on-top of an exisitng modulation scheme.

    properties (SetAccess = immutable, GetAccess = public)
        tag (1, 1)                % Tag associated with modulation scheme
        slot_select (1, 1) double % Slot in frame to transmit packet
        num_slots (1, 1) double   % Number of total slots in a frame
        slot_size (1, 1) double   % Number of `simsteps` in a slot
        frame_size                % Number of `simsteps` in frame
    end

    properties (SetAccess = private, GetAccess = public)
        curr_step (1, 1) double
        prev_data (1, :)
    end

% ----------------------------------------------------------------------- %

    methods
        function this = PacketTimeDivisionOrchestrator(slot_size, num_slots, slot_select, tag)
            %PACKETTIMEORCHESTRATOR constructor.

            arguments
                slot_size (1, 1) double {mustBePositive}
                num_slots (1, 1) double {mustBePositive}
                slot_select (1, 1) double {mustBePositive}
                tag (1, 1) {mustBeA(tag, [ ...
                                          "bherber.thesis.tags.OOKTag", ...
                                          "bherber.thesis.tags.FSKTag", ...
                                          "bherber.thesis.tags.FreqHopTag" ...
                                         ])}
            end

            this.num_slots = num_slots;
            this.slot_size = slot_size;
            this.frame_size = num_slots * slot_size;
            this.slot_select = slot_select;

            this.curr_step = 0;
            this.tag = tag;

            this.prev_data = zeros(1, this.tag.params.simstep_sz);

        end

% ----------------------------------------------------------------------- %

        function res = step(this)
            % STEP through one simulation step with the orchestrator.

            this.curr_step = this.curr_step + 1;
            curr_step_in_frame = mod(this.curr_step, this.frame_size);
            if curr_step_in_frame == 0; curr_step_in_frame = this.frame_size; end
            curr_slot = ceil(curr_step_in_frame / this.slot_size);

            delay = 2 * this.tag.distance / physconst("Lightspeed");
            delay_samples = ceil(delay * this.tag.params.Fs);

            if (curr_slot ~= this.slot_select)
                this.tag.deactivate();
                res = this.tag.step();

                if (mod(curr_step_in_frame, this.slot_size) == 1) ...
                    && ((curr_slot - 1) == this.slot_select)
                    res(1:delay_samples) = this.prev_data(1:delay_samples);
                end

                return;
            end

            this.tag.activate();

            if (mod(curr_step_in_frame, this.slot_size) == 1)
                res = this.prev_data;
                res(1:delay_samples) = 0;
            elseif (mod(curr_step_in_frame, this.slot_size) == 0)
                res = this.tag.step();
                this.prev_data = this.tag.step();
            end
        end

% ----------------------------------------------------------------------- %

    function res = deframeify(this, frame)
        %DEFRAMEIFY a given frame into the packet that the assocated tag
        %   had sent. Assumes external adjustment for Time-of-Flight.

        arguments
            this (1, 1) bherber.thesis.tags.PacketTimeDivisionTag
            frame (:, 1)
        end

        slots = reshape(frame, [], this.num_slots);
        res = slots(:, this.slot_select).';

    end

% ----------------------------------------------------------------------- %

    end
end