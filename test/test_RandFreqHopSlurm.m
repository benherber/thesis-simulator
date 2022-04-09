clear; clc;
addpath(genpath(".."));
import bherber.thesis.PPersistPacketSimulator bherber.thesis.TagType
filename = sprintf("../data/randfreqhopcollisions::%s.mat", datestr(now, "mm-dd-yy-HH:MM:SS"));
save(filename);

CollisionFactories;

% if ~isempty(getenv("SLURM_ARRAY_TASK_ID"))
%     num_channels = str2double(getenv("SLURM_ARRAY_TASK_ID"));
% else
%     num_channels = 350;
% end
num_channels = 350;

num_packets = 1000;
tag_type = TagType.FREQ_HOP;
ppersist = 1;
params = make_params(num_channels);

max_activity = 350;

for num_tags = ceil(linspace(1, max_activity, 50))

    clashes = zeros(1, num_threads);
    packets_per_thread = ceil(num_packets / num_threads);
    
    for idx = 1:num_threads
        sim = make_sim(packets_per_thread, num_tags, tag_type, ppersist, params);
        for ignore__ = 1:packets_per_thread
            clashes(idx) = clashes(idx) + double(sim.clash_in_next_frame());
        end
    end

    clash_prob = sum(clashes, "all") / (packets_per_thread * num_threads * num_tags);
    save_out(filename, struct(sprintf("tags%dnumchannels%d", num_tags, num_channels), clash_prob));
    
end

%%
function save_out(filename, structure)
    save(filename, "-append", "-struct", "structure");
end
