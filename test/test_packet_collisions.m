import bherber.thesis.PacketSimulator ...
    bherber.thesis.SimulationConstants ...
    bherber.thesis.TagType

clear; clc;

id_len = 16;
data_len = 0;

make_params = @(num_channels) bherber.thesis.SimulationConstants( ...
                2.4e9, ...   "Fc"
                100e3, ...   "Fb_base"
                100e3, ...   "Fb_step"
                0.5 / (1 / (100e3)), ...   "Fb_channel_spacing"
                num_channels, ...       "num_channels"
                4, ...       "fs_granularity"
                (id_len + data_len), ... "num_symbs"
                1e5, ...     "symb_freq"
                1, ...       "sim_sym_ratio"
                sqrt(8), ... "amplitude"
                4,  ...      "num_elements"
                3, ...        "pattern_len"
                2, ... "m_ary_modulation"
                id_len ... "preamble_len"
            );

make_sim = @(num_packets, num_tags, type, num_slots, params) ...
    PacketSimulator(num_packets, repelem([1;0;0], 1, num_tags), type, num_slots, params);

if isempty(gcp("nocreate"))
  pc = parcluster("local");
  pc.JobStorageLocation = "/tmp/";
  parpool(pc, 6)
end
pool = gcp;
num_threads = pool.NumWorkers;

%% OOK-TDM

clearvars -except make_sim make_params num_threads

import bherber.thesis.PacketSimulator ...
    bherber.thesis.SimulationConstants ...
    bherber.thesis.TagType

num_tags = 2;
num_channels = 1;
num_slots = 4;
num_packets = 10;
tag_type = TagType.OOK;

params = make_params(num_channels);
sim = make_sim(num_packets, num_tags, tag_type, num_slots, params);

clashes = 0;

for idx = 1:num_packets
    clashes = clashes + double(sim.clash_in_next_frame());
end

ook_clash_percentage = clashes / num_packets;


%% FH-TDM

clearvars -except make_sim make_params num_threads

import bherber.thesis.PacketSimulator ...
    bherber.thesis.SimulationConstants ...
    bherber.thesis.TagType

num_tags = 2;
num_channels = 2;
num_slots = 2;
num_packets = 60;
tag_type = TagType.FREQ_HOP;

params = make_params(num_channels);

clashes = zeros(1, num_threads);
packets_per_thread = ceil(num_packets / num_threads);

parfor idx = 1:num_threads
    sim = make_sim(packets_per_thread, num_tags, tag_type, num_slots, params);
    for ignore__ = 1:packets_per_thread
        clashes(idx) = clashes(idx) + double(sim.clash_in_next_frame());
    end
end

fsk_clash_percentage = sum(clashes, "all") / (packets_per_thread * num_threads);


