clear; clc;
addpath(genpath(".."));
import bherber.thesis.PacketSimulator bherber.thesis.TagType
filename = sprintf("../data/tdmcollisions::%s.mat", datestr(now, "mm-dd-yy-HH:MM:SS"));
save(filename);

CollisionFactories;

if ~isempty(getenv("SLURM_ARRAY_TASK_ID"))
    num_slots = str2double(getenv("SLURM_ARRAY_TASK_ID"));
else
    num_slots = 100;
end

num_packets = 1000;
tag_type = TagType.FREQ_HOP;
num_channels = 1;
params = make_params(num_channels);

max_activity = 1;

for num_tags = ceil(linspace(1, max_activity * num_slots, 10))

    clashes = zeros(1, num_threads);
    packets_per_thread = ceil(num_packets / num_threads);
    
    parfor idx = 1:num_threads
        sim = make_sim(packets_per_thread, num_tags, tag_type, num_slots, params);
        for ignore__ = 1:packets_per_thread
            clashes(idx) = clashes(idx) + double(sim.clash_in_next_frame());
        end
    end

    clash_percentage = sum(clashes, "all") / (packets_per_thread * num_threads * num_tags);
    save_out(filename, struct(sprintf("tags%numslots%d", num_tags, num_channels), clash_percentage));
    
end

%%
function save_out(filename, structure)
    save(filename, "-append", "-struct", "structure");
end
