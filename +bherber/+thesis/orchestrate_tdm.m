function [frame, frame_sz] = orchestrate_tdm(num_slots, slot_select, tag)
    %ORCHESTRATETDM utilizing a given tag implementation.
    %   Orchestrate a tag to transmit during a given time slot in order to
    %   implement a variable pack-length Time Division Multiplexing scheme.
    %   Packet size depends of parameters fed to tag object.
    %ARGUMENTS:
    %   - NUMSLOTS in a given multiplexed frame.
    %   - SLOTSELECT for the provided tag to be active in.
    %   - TAG implementation with a given modulation scheme and packet
    %     length.

    arguments
        num_slots (1, 1) double {mustBePositive}
        slot_select (1, 1) double {mustBePositive}
        tag (1, 1) bherber.thesis.tags.Tag
    end

    if (tag.curr_step ~= 0)
        error("Tag at non-zero step (%d): Tag must be freshly constucted", tag.curr_step);
    end

    if (slot_select > num_slots)
        error("Selected slot (%d) higher than that available (%d)", slot_select, num_slots);
    end

    num_steps = int32((tag.params.num_symbs + 1) * tag.params.sim_sym_ratio);
    slot_sz = int32(num_steps * tag.params.simstep_sz);
    all_steps = zeros(int32(tag.params.simstep_sz), num_steps);

    for idx = 1:num_steps
        all_steps(:, idx) = tag.step();
    end

    frame = zeros(slot_sz, num_slots);
    frame(:, slot_select) = all_steps(:);
    frame = frame(:).';
    frame_sz = slot_sz * num_slots;

end