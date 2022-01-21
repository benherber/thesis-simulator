clc; clear;

consts_struct = struct( ...
    "Fc", 2.4e9, ...
    "Fb_base", 75e3, ...
    "Fb_step", 25e3, ...
    "fs_granularity", 3, ...
    "num_symbs", 20, ...
    "symb_freq", 100e3, ...
    "sim_sym_ratio", 20, ...
    "amplitude", 0.1, ...
    "num_elements", 4 ...
    );
consts_cell = struct2cell(consts_struct);
simconsts = SimulationConstants(consts_cell{:});

channel = @(given) awgn(given, 10, "measured", "linear");

tags = [1000;1000;1000];
tag_modes = [TagType.FSK_LO];

sim = Simulator(tags, tag_modes, channel, simconsts);

modded_steps = zeros(simconsts.simstep_sz, (simconsts.num_symbs * simconsts.sim_sym_ratio));
for idx = 1:(simconsts.num_symbs * simconsts.sim_sym_ratio)
    modded_steps(:, idx) = sim.step();
end

modded_symbs = reshape(modded_steps, [simconsts.symb_sz, simconsts.num_symbs]);

sim.auto_align(modded_symbs(:));

encoded_bits = sim.bits(1, :).';
res_bits = zeros(numel(encoded_bits), 1);
symb_times = reshape(sim.time, [simconsts.symb_sz, simconsts.num_symbs]);
carrier_split = reshape(sim.carrier, [simconsts.symb_sz, simconsts.num_symbs]);
f1 = simconsts.fsk_channel0.f1;
f0 = simconsts.fsk_channel0.f0;
advanced_sig = delayseq(modded_symbs(:), ceil(sim.tag_delays(1) / 2));
advanced_sig = reshape(advanced_sig, [simconsts.symb_sz, simconsts.num_symbs]);
parfor idx = 1:simconsts.num_symbs
    res_bits(idx) = Tag.fsk_demodulate(advanced_sig(:, idx), carrier_split(:, idx), ... 
        symb_times(:, idx), f1, f0);
end

fprintf("delay: %f\n", sim.tag_delays(1))
fprintf("FSK (lo) Bit Error \t: %%%0.2f\n", (biterr(encoded_bits, res_bits) / simconsts.num_symbs) * 100);