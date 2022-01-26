%% Set-Up

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

% Channel
channel = @(given) awgn(given, 10, "measured", "linear");

% Tag Locations
num_tags = 10;
max_distance = 100; % Metes
tags = zeros(3, num_tags);
tag_modes = repmat(TagType.FSK_LO, 1, num_tags);
for idx = 1:num_tags
    x = max_distance;
    y = max_distance;
    z = max_distance;
    while (x^2 + y^2 + z^2) > max_distance
        x = rand * max_distance;
        y = rand * max_distance;
        z = rand * max_distance;
    end
    tags(:, idx) = [x, y, z];
    
    if logical(randi([0,1])) % Randomize frequency channels
        tag_modes(idx) = TagType.FSK_HI;
    end
end

% Create Simulator
sim = Simulator(tags, tag_modes, channel, simconsts);

%% Run Simulation

% Evaluate
modded_steps = zeros(simconsts.simstep_sz, (simconsts.num_symbs * simconsts.sim_sym_ratio));
for idx = 1:(simconsts.num_symbs * simconsts.sim_sym_ratio)
    modded_steps(:, idx) = sim.step();
end

modded_symbs = reshape(modded_steps, [simconsts.symb_sz, simconsts.num_symbs]);

% Find delays
sim.auto_align(modded_symbs(:));

% Demodulate each stream
running_ber = 0;
for jdx = 1:num_tags
    fprintf("Tag delay %d: %f\n", jdx, sim.tag_delays(jdx))
    encoded_bits = sim.bits(jdx, :).';
    res_bits = zeros(numel(encoded_bits), 1);
    symb_times = reshape(sim.time, [simconsts.symb_sz, simconsts.num_symbs]);
    carrier_split = reshape(sim.carrier, [simconsts.symb_sz, simconsts.num_symbs]);
    if sim.tags(jdx).mode == TagType.FSK_LO
        f1 = simconsts.fsk_channel0.f1;
        f0 = simconsts.fsk_channel0.f0;
    else
        f1 = simconsts.fsk_channel1.f1;
        f0 = simconsts.fsk_channel1.f0;
    end
    advanced_sig = delayseq(modded_symbs(:), -1 * sim.tag_delays(jdx));
    advanced_sig = reshape(advanced_sig, [simconsts.symb_sz, simconsts.num_symbs]);
    parfor idx = 1:simconsts.num_symbs
        res_bits(idx) = Tag.fsk_demodulate(advanced_sig(:, idx), carrier_split(:, idx), ... 
            symb_times(:, idx), f1, f0);
    end
    running_ber = running_ber + biterr(encoded_bits, res_bits);
end

fprintf("FSK (lo) Bit Error \t: %%%0.2f\n", (biterr(encoded_bits, res_bits) / simconsts.num_symbs) * 100);