clc; clear;

consts_struct = struct( ...
    "Fc", 2.4e9, ...
    "Fb_base", 75e3, ...
    "Fb_step", 25e3, ...
    "fs_granularity", 3, ...
    "num_symbs", 1000, ...
    "symb_freq", 100e3, ...
    "sim_sym_ratio", 10, ...
    "amplitude", 0.1, ...
    "num_elements", 4 ...
    );
consts_cell = struct2cell(consts_struct);
simconsts = SimulationConstants(consts_cell{:});

channel = @(given) awgn(given, 10, "measured", "linear");

tags = [[1;1;1],[1;-1;1]];
tag_modes = [TagType.FSK_LO, TagType.FSK_HI];

sim = Simulator(tags, tag_modes, channel, simconsts);