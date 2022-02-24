clc; clear;

addpath(genpath(".."));
import bherber.thesis.*

consts_struct = struct( ...
    "Fc", 2.4e9, ...
    "Fb_base", 75e3, ...
    "Fb_step", 25e3, ...
    "Fb_channel_spacing", 25e3, ...
    "fs_granularity", 3, ...
    "num_symbs", 100, ...
    "symb_freq", 100e3, ...
    "sim_sym_ratio", 10, ...
    "amplitude", 0.1, ...
    "num_elements", 4 ...
    );
consts_cell = struct2cell(consts_struct);
simconsts = SimulationConstants(consts_cell{:});

channel = @(signal) awgn_channel(signal, 10);
% channel = @(signal, tag_modes, distances, params) signal;

tags = [1000;1000;1000];
tag_modes = [TagType.FSK_LO];

sim = Simulator(tags, tag_modes, channel, simconsts);

modded_steps = zeros(simconsts.simstep_sz, (simconsts.num_symbs * simconsts.sim_sym_ratio));
for idx = 1:((simconsts.num_symbs + 1) * simconsts.sim_sym_ratio)
    modded_steps(:, idx) = sim.step();
end

noisy_modded_steps = channel(modded_steps(:));
modded_symbs = reshape(noisy_modded_steps, [simconsts.symb_sz, simconsts.num_symbs + 1]);

sim.auto_align();
fprintf("delay: %f\n", sim.tag_delays(1))

correct = @(sig, samples) [sig(samples:end).'; zeros(samples, 1)];

encoded_bits = sim.bits(1, :);
res_bits = zeros(size(encoded_bits));
mod_delay = floor(1 * sim.tag_delays(1));
stop = simconsts.num_symbs * simconsts.symb_sz + mod_delay - 1;
t = sim.time(mod_delay:stop);
symb_times = reshape(t, [simconsts.symb_sz, simconsts.num_symbs]);
carr = sim.carrier(mod_delay:stop);
carrier_split = reshape(carr, [simconsts.symb_sz, simconsts.num_symbs]);
f1 = simconsts.fsk_channel0.f1;
f0 = simconsts.fsk_channel0.f0;
sigcomplete = modded_symbs(:);
sig = sigcomplete(mod_delay:stop);
advanced_sig = reshape(sig, [simconsts.symb_sz, simconsts.num_symbs]);
parfor idx = 1:simconsts.num_symbs
    res_bits(idx) = Tag.fsk_demodulate(advanced_sig(:, idx), carrier_split(:, idx), ... 
        symb_times(:, idx), f1, f0);
end

fprintf("FSK (lo) Bit Error \t: %%%0.2f\n", (biterr(encoded_bits, res_bits) / simconsts.num_symbs) * 100);