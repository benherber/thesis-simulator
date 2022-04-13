import bherber.thesis.PacketSimulator ...
    bherber.thesis.SimulationConstants ...
    bherber.thesis.TagType

clear; clc;

id_len = 16;
data_len = 0;

params = bherber.thesis.SimulationConstants( ...
                2.4e9, ...   "Fc"
                100e3, ...   "Fb_base"
                100e3, ...   "Fb_step"
                0.5 / (1 / (100e3)), ...   "Fb_channel_spacing"
                5, ...       "num_channels"
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

%% OOK-TDM

clearvars -except params

import bherber.thesis.PacketSimulator ...
    bherber.thesis.SimulationConstants ...
    bherber.thesis.TagType

tag_locs = [[1;0;0],[1;0;0]];
tag_types = TagType.OOK;
num_slots = 4;
num_packets = 8;

sim = PacketSimulator(num_packets, tag_locs, tag_types, num_slots, params);

frames = zeros(int32(sim.num_slots * (sim.params.num_symbs + 1) * ...
        sim.params.sim_sym_ratio * sim.params.simstep_sz), num_packets);

for idx = 1:num_packets
    frames(:, idx) = sim.step();
    symbols = reshape(sim.deframeify(1, frames(:, idx)), [], params.num_symbs + 1);
    resbits = zeros(params.num_symbs, 1);
    for jdx = 1:params.num_symbs
        resbits(jdx) = sim.tags(1).demodulate(jdx, symbols(:, jdx).');
    end
    fprintf("OOKTag BER: %d\n", biterr(sim.tags(1).bits(1:params.num_symbs), resbits))
end

frame_total = frames(:).';
time = (0:(length(frame_total) - 1)) / params.Fs;

figure();
plot(time, frame_total);
title("OOK-TDM Modulated Singal", FontSize=20);
xlabel("Transmission Frames", FontSize=18);
ylabel("Amplitude (V)", FontSize=18);
xlines = zeros(1, num_packets);
xlines_labels = strings(1, num_packets);
frame_sz = (sim.num_slots * (sim.params.num_symbs + 1) * ...
    sim.params.sim_sym_ratio * sim.params.simstep_sz) / params.Fs;
for idx = 1:num_packets
    xlines(idx) = idx * frame_sz;
    xlines_labels(idx) = sprintf("Frame%d", idx);
end
xline(xlines, "--k");
xlim([0, xlines(end)]);
xticks(xlines);
xticklabels(xlines_labels);

%% FH-TDM

clearvars -except params

import bherber.thesis.PacketSimulator ...
    bherber.thesis.SimulationConstants ...
    bherber.thesis.TagType

tag_locs = [[1;0;0],[1;0;0]];
tag_types = TagType.FREQ_HOP;
num_slots = 2;
num_packets = 10;

sim = PacketSimulator(num_packets, tag_locs, tag_types, num_slots, params);

frames = zeros(int32(sim.num_slots * (sim.params.num_symbs + 1) * ...
        sim.params.sim_sym_ratio * sim.params.simstep_sz), num_packets);

for idx = 1:num_packets
    frames(:, idx) = sim.step();
    symbols = reshape(sim.deframeify(1, frames(:, idx)), [], params.num_symbs + 1);
    resbits = zeros(params.num_symbs, 1);
    for jdx = 1:params.num_symbs
        resbits(jdx) = sim.tags(1).demodulate(jdx, symbols(:, jdx).');
    end
    fprintf("FSKTag BER: %d\n", biterr(sim.tags(1).bits(1:params.num_symbs), resbits))
end

frame_total = frames(:).';
time = (0:(length(frame_total) - 1)) / params.Fs;

% PLOT
figure();
plot(time, frame_total);
title("Slow FH-TDM Modulated Singal", FontSize=20);
xlabel("Transmission Frames", FontSize=18);
ylabel("Amplitude (V)", FontSize=18);
xlines = zeros(1, num_packets);
xlines_labels = strings(1, num_packets);
frame_sz = (sim.num_slots * (sim.params.num_symbs + 1) * ...
    sim.params.sim_sym_ratio * sim.params.simstep_sz) / params.Fs;
for idx = 1:num_packets
    xlines(idx) = idx * frame_sz;
    xlines_labels(idx) = sprintf("Frame%d", idx);
end
xline(xlines, "--k");
xlim([0, xlines(end)]);
xticks(xlines);
xticklabels(xlines_labels);

%% SPECTROGRAM
nfft = 10000000;
f = figure();
spectrogram(frame_total, hamming(frame_sz / ...
    (((params.num_symbs + 1) * params.symb_sz) / params.symb_freq)), ...
    0, nfft, params.Fs, "yaxis", "power");
axs = f.Children;
data_objs = axs(2).Children;
tmpYData = data_objs.YData;
data_objs(1).YData = (tmpYData - 2.4) * 1e6;
ylim([-120, -70]);
