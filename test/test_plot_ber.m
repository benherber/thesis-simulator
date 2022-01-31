%% Setup

clear; clc;

addpath(genpath(".."));
import bherber.thesis.*

consts_struct = struct( ...
    "Fc", 2.4e9, ...
    "Fb_base", 75e3, ...
    "Fb_step", 25e3, ...
    "Fb_channel_spacing", 50e3, ...
    "fs_granularity", 3, ...
    "num_symbs", 1000, ...
    "symb_freq", 100e3, ...
    "sim_sym_ratio", 10, ...
    "amplitude", 0.1, ...
    "num_elements", 4 ...
    );
consts_cell = struct2cell(consts_struct);
simconsts = SimulationConstants(consts_cell{:});

NUM_RUNS = int8(100);

% Get time
time = (0:(1 / simconsts.Fs):(simconsts.total_time - (1 / simconsts.Fs)));

% Carrier signal
carrier = complex(simconsts.amplitude * sin(complex(2 * pi * simconsts.Fc * time)));

%% Calculate BERs

snrs = (0:0.001:0.05);

ook_bers = NaN(NUM_RUNS, length(snrs));
fsk_bers = NaN(NUM_RUNS, length(snrs));
dual_bers = NaN(NUM_RUNS, length(snrs));

for run = 1:NUM_RUNS
    % Get random data signal
    bits = randi([0, 1], simconsts.num_symbs, 1);
    data = repelem(bits, simconsts.symb_sz).';

    % Calculate partial BER
    parfor idx = 1:length(snrs)
        ook_bers(run, idx) = ook_ber(snrs(idx), time, carrier, data, bits, simconsts);
        fsk_bers(run, idx) = fsk_ber(snrs(idx), TagType.FSK_LO, time, carrier, data, bits, simconsts);
    end
end
save(sprintf("data/bers::%s.mat", datestr(now, "mm-dd-yy-HH:MM:SS")), ...
    "ook_bers", "fsk_bers");

%% Plot

fig = figure();
hold on
semilogy(snrs, mean(ook_bers, 1, "omitnan"));
semilogy(snrs, mean(fsk_bers, 1, "omitnan"));
legend(["OOK", "FSK", "Dual FSK"]);
savefig(sprintf("plots/ber::%s.fig", datestr(now, "mm-dd-yy-HH:MM:SS")));
hold off
close(fig);

%% OOK BER Function
function res = ook_ber(snr, time, carrier, data, bits, params)
    
    arguments
        snr {mustBeFinite}
        time {mustBeNonmissing, mustBeFinite}
        carrier {mustBeNonmissing, mustBeFinite}
        data {mustBeNonmissing, mustBeFinite}
        bits
        params bherber.thesis.SimulationConstants
    end

    import bherber.thesis.*

    channel = @(signal) awgn(signal, snr, "measured", "linear");

    tag = Tag(0, 0, 0, TagType.OOK, time, carrier, data, channel, params);

    % Modulate
    modded_steps = zeros(params.simstep_sz, (params.num_symbs * params.sim_sym_ratio));
    for idx = 1:(params.num_symbs * params.sim_sym_ratio)
        modded_steps(:, idx) = tag.step();
    end
    
    % Demodulate
    modded_symbs = reshape(modded_steps, [params.symb_sz, params.num_symbs]);
    
    res_bits = zeros(numel(bits), 1);
    symb_times = reshape(time, [params.symb_sz, params.num_symbs]);
    carrier_split = reshape(carrier, [params.symb_sz, params.num_symbs]);
    parfor idx = 1:params.num_symbs
        res_bits(idx) = Tag.ook_demodulate(modded_symbs(:, idx), carrier_split(:, idx), ... 
             symb_times(:, idx));
    end

    % Calculate BER
    res = (biterr(bits, res_bits) / params.num_symbs);

end

%% FSK BER Function
function res = fsk_ber(snr, tag_type, time, carrier, data, bits, params)
    
    arguments
        snr {mustBeFinite}
        tag_type bherber.thesis.TagType
        time {mustBeNonmissing, mustBeFinite}
        carrier {mustBeNonmissing, mustBeFinite}
        data {mustBeNonmissing, mustBeFinite}
        bits
        params bherber.thesis.SimulationConstants
    end

    import bherber.thesis.*

    channel = @(signal) awgn(signal, snr, "measured", "linear");

    tag = Tag(0, 0, 0, tag_type, time, carrier, data, channel, params);

    % Modulate
    modded_steps = zeros(params.simstep_sz, (params.num_symbs * params.sim_sym_ratio));
    for idx = 1:(params.num_symbs * params.sim_sym_ratio)
        modded_steps(:, idx) = tag.step();
    end
    
    % Demodulate
    modded_symbs = reshape(modded_steps, [params.symb_sz, params.num_symbs]);
    
    res_bits = zeros(numel(bits), 1);
    symb_times = reshape(time, [params.symb_sz, params.num_symbs]);
    carrier_split = reshape(carrier, [params.symb_sz, params.num_symbs]);

    if tag_type == TagType.FSK_HI
        f1 = params.fsk_channel1.f1;
        f0 = params.fsk_channel1.f0;
    else
        f1 = params.fsk_channel0.f1;
        f0 = params.fsk_channel0.f0;
    end

    parfor idx = 1:params.num_symbs
        res_bits(idx) = Tag.fsk_demodulate(modded_symbs(:, idx), carrier_split(:, idx), ... 
            symb_times(:, idx), f1, f0);
    end

    % Calculate BER
    res = (biterr(bits, res_bits) / params.num_symbs);

end

%% Dual FSK BER Function
function res = dual_fsk_ber(snr, time, carrier, data, bits, params)
    
    arguments
        snr {mustBeFinite}
        time {mustBeNonmissing, mustBeFinite}
        carrier {mustBeNonmissing, mustBeFinite}
        data {mustBeNonmissing, mustBeFinite}
        bits
        params bherber.thesis.SimulationConstants
    end

    import bherber.thesis.*

    channel = @(signal) awgn(signal, snr, "measured", "linear");

    % Get random data signals
    bits1 = bits;
    data1 = data;
    bits2 = double(~bits);
    data2 = repelem(bits2, params.symb_sz).'; 
    
    tag1 = Tag(0, 0, 0, TagType.FSK_LO, time, carrier, data1, channel, params);
    tag2 = Tag(0, 0, 0, TagType.FSK_HI, time, carrier, data2, channel, params);
    
    % Modulate
    modded_steps = zeros(params.simstep_sz, (params.num_symbs * params.sim_sym_ratio));
    for idx = 1:(params.num_symbs * params.sim_sym_ratio)
        modded_steps(:, idx) = tag2.step() + tag1.step();
    end
    
    % Demodulate
    modded_symbs = reshape(modded_steps, [params.symb_sz, params.num_symbs]);
    
    res_bits1 = zeros(numel(bits1), 1);
    res_bits2 = zeros(numel(bits2), 1);
    symb_times = reshape(time, [params.symb_sz, params.num_symbs]);
    carrier_split = reshape(carrier, [params.symb_sz, params.num_symbs]);
    channel0f1 = params.fsk_channel0.f1;
    channel0f0 = params.fsk_channel0.f0;
    channel1f1 = params.fsk_channel1.f1;
    channel1f0 = params.fsk_channel1.f0;
    parfor idx = 1:params.num_symbs
        res_bits1(idx) = Tag.fsk_demodulate(modded_symbs(:, idx), carrier_split(:, idx), ... 
            symb_times(:, idx), channel0f1, channel0f0);
        res_bits2(idx) = Tag.fsk_demodulate(modded_symbs(:, idx), carrier_split(:, idx), ... 
            symb_times(:, idx), channel1f1, channel1f0);
    end
    
    % Calculate BER
    dualchannel_fsk_ber1 = (biterr(bits1, res_bits1) / params.num_symbs);
    dualchannel_fsk_ber2 = (biterr(bits2, res_bits2) / params.num_symbs);
    
    res = mean([dualchannel_fsk_ber1, dualchannel_fsk_ber2]);

end
