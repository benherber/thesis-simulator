clear; clc;

addpath(genpath(".."));
import bherber.thesis.*

consts_struct = struct( ...
    "Fc", 2.4e9, ...
    "Fb_base", 75e3, ...
    "Fb_step", 25e3, ...
    "Fb_channel_spacing", 50e3, ...
    "fs_granularity", 3, ...
    "num_symbs", 100, ...
    "symb_freq", 100e3, ...
    "sim_sym_ratio", 10, ...
    "amplitude", 0.1, ...
    "num_elements", 4 ...
    );
consts_cell = struct2cell(consts_struct);
simconsts = SimulationConstants(consts_cell{:});

% Get time
time = (0:(1 / simconsts.Fs):(simconsts.total_time - (1 / simconsts.Fs)));

% Carrier signal
carrier = complex(simconsts.amplitude * sin(complex(2 * pi * simconsts.Fc * time)));

% Channel Definition
channel = @(given) awgn(given, 0.002, "measured", "linear");

% Get random data signal
bits = randi([0, 1], 1, simconsts.num_symbs);
data = repelem(bits, simconsts.symb_sz); 

%% OOK Tag

tag = Tag(0, 0, 0, TagType.OOK, time, carrier, data, channel, simconsts);

% Modulate
modded_steps = zeros(simconsts.simstep_sz, (simconsts.num_symbs * simconsts.sim_sym_ratio));
for idx = 1:(simconsts.num_symbs * simconsts.sim_sym_ratio)
    modded_steps(:, idx) = tag.step();
end

% Demodulate
modded_symbs = reshape(modded_steps, [simconsts.symb_sz, simconsts.num_symbs]);

res_bits = zeros(1, numel(bits));
symb_times = reshape(time, [simconsts.symb_sz, simconsts.num_symbs]);
carrier_split = reshape(carrier, [simconsts.symb_sz, simconsts.num_symbs]);
parfor idx = 1:simconsts.num_symbs
    res_bits(idx) = Tag.ook_demodulate(modded_symbs(:, idx), carrier_split(:, idx), ... 
         symb_times(:, idx));
end

fprintf("OOK Bit Error \t\t: %%%0.2f\n", (biterr(bits, res_bits) / simconsts.num_symbs) * 100);

%% FSK Tag (LO)

tag = Tag(0, 0, 0, TagType.FSK_LO, time, carrier, data, channel, simconsts);

% Modulate
modded_steps = zeros(simconsts.simstep_sz, (simconsts.num_symbs * simconsts.sim_sym_ratio));
for idx = 1:(simconsts.num_symbs * simconsts.sim_sym_ratio)
    modded_steps(:, idx) = tag.step();
end

% Demodulate
modded_symbs = reshape(modded_steps, [simconsts.symb_sz, simconsts.num_symbs]);

res_bits = zeros(1, numel(bits));
symb_times = reshape(time, [simconsts.symb_sz, simconsts.num_symbs]);
carrier_split = reshape(carrier, [simconsts.symb_sz, simconsts.num_symbs]);
f1 = simconsts.fsk_channel0.f1;
f0 = simconsts.fsk_channel0.f0;
parfor idx = 1:simconsts.num_symbs
    res_bits(idx) = Tag.fsk_demodulate(modded_symbs(:, idx), carrier_split(:, idx), ... 
        symb_times(:, idx), f1, f0);
end

fprintf("FSK (lo) Bit Error \t: %%%0.2f\n", (biterr(bits, res_bits) / simconsts.num_symbs) * 100);

%% FSK Tag (HI)

tag = Tag(0, 0, 0, TagType.FSK_HI, time, carrier, data, channel, simconsts);

% Modulate
modded_steps = zeros(simconsts.simstep_sz, (simconsts.num_symbs * simconsts.sim_sym_ratio));
for idx = 1:(simconsts.num_symbs * simconsts.sim_sym_ratio)
    modded_steps(:, idx) = tag.step();
end

% Demodulate
modded_symbs = reshape(modded_steps, [simconsts.symb_sz, simconsts.num_symbs]);

res_bits = zeros(1, numel(bits));
symb_times = reshape(time, [simconsts.symb_sz, simconsts.num_symbs]);
carrier_split = reshape(carrier, [simconsts.symb_sz, simconsts.num_symbs]);
f1 = simconsts.fsk_channel1.f1;
f0 = simconsts.fsk_channel1.f0;
parfor idx = 1:simconsts.num_symbs
    res_bits(idx) = Tag.fsk_demodulate(modded_symbs(:, idx), carrier_split(:, idx), ... 
        symb_times(:, idx), f1, f0);
end

fprintf("FSK (hi) Bit Error \t: %%%0.2f\n", (biterr(bits, res_bits) / simconsts.num_symbs) * 100);

%% Dual-Channel FSK Tags

% Get random data signals
bits1 = randi([0, 1], 1, simconsts.num_symbs);
data1 = repelem(bits1, simconsts.symb_sz); 
bits2 = randi([0, 1], 1, simconsts.num_symbs);
data2 = repelem(bits2, simconsts.symb_sz); 

tag1 = Tag(0, 0, 0, TagType.FSK_LO, time, carrier, data1, channel, simconsts);
tag2 = Tag(0, 0, 0, TagType.FSK_HI, time, carrier, data2, channel, simconsts);

% Modulate
modded_steps = zeros(simconsts.simstep_sz, (simconsts.num_symbs * simconsts.sim_sym_ratio));
for idx = 1:(simconsts.num_symbs * simconsts.sim_sym_ratio)
    modded_steps(:, idx) = tag2.step() + tag1.step();
end
f_modulation = modded_steps(:);

% Demodulate
modded_symbs = reshape(modded_steps, [simconsts.symb_sz, simconsts.num_symbs]);

res_bits1 = zeros(1, numel(bits));
res_bits2 = zeros(1, numel(bits));
symb_times = reshape(time, [simconsts.symb_sz, simconsts.num_symbs]);
carrier_split = reshape(carrier, [simconsts.symb_sz, simconsts.num_symbs]);
channel0f1 = simconsts.fsk_channel0.f1;
channel0f0 = simconsts.fsk_channel0.f0;
channel1f1 = simconsts.fsk_channel1.f1;
channel1f0 = simconsts.fsk_channel1.f0;
parfor idx = 1:simconsts.num_symbs
    res_bits1(idx) = Tag.fsk_demodulate(modded_symbs(:, idx), carrier_split(:, idx), ... 
        symb_times(:, idx), channel0f1, channel0f0);
    res_bits2(idx) = Tag.fsk_demodulate(modded_symbs(:, idx), carrier_split(:, idx), ... 
        symb_times(:, idx), channel1f1, channel1f0);
end

dualchannel_fsk_ber1 = 100 * (biterr(bits1, res_bits1) / simconsts.num_symbs);
dualchannel_fsk_ber2 = 100 * (biterr(bits2, res_bits2) / simconsts.num_symbs);

fprintf("FSK (dual) Bit Error \t: %%%0.2f ('lo': %%%0.2f, 'hi': %%%0.2f)\n", ...
    mean([dualchannel_fsk_ber1, dualchannel_fsk_ber2]), dualchannel_fsk_ber1, dualchannel_fsk_ber2);

%% Calculate FFT

nfft = 100000000;
f_fft = fftshift(fft(f_modulation, nfft));

% Carrier fft
carrier_fft = fftshift(fft(carrier, nfft));

% Frequency
f = simconsts.Fs * (-nfft/2 : nfft/2 - 1) / nfft;
scaled_f = (f - 2.4e9) / 1e3;

%% FFT Plot

f1 = figure();
hold on
scaled_fft = 2 * abs(f_fft) / length(f_modulation);
scaled_carrier_fft = 2 * abs(carrier_fft) / length(carrier);
plot(scaled_f, scaled_fft);
plot(scaled_f, scaled_carrier_fft);
label_pts = zeros(1, 4);
[ ~, label_pts(1) ] = min(abs(scaled_f - simconsts.fsk_channel0.f0 * 1e-3));
[ ~, label_pts(2) ] = min(abs(scaled_f - simconsts.fsk_channel0.f1 * 1e-3));
[ ~, label_pts(3) ] = min(abs(scaled_f - simconsts.fsk_channel1.f0 * 1e-3));
[ ~, label_pts(4) ] = min(abs(scaled_f - simconsts.fsk_channel1.f1 * 1e-3));
plot(scaled_f(label_pts), scaled_fft(label_pts), ".r");
text(scaled_f(label_pts) + 10, scaled_fft(label_pts), ["ch0(f0)", "ch0(f1)", "ch1(f0)", "ch1(f1)"], ...
    FontWeight="Bold", Rotation=60)
[ ~, label_pts(1) ] = min(abs(scaled_f + simconsts.fsk_channel0.f0 * 1e-3));
[ ~, label_pts(2) ] = min(abs(scaled_f + simconsts.fsk_channel0.f1 * 1e-3));
[ ~, label_pts(3) ] = min(abs(scaled_f + simconsts.fsk_channel1.f0 * 1e-3));
[ ~, label_pts(4) ] = min(abs(scaled_f + simconsts.fsk_channel1.f1 * 1e-3));
plot(scaled_f(label_pts), scaled_fft(label_pts), ".r");
text(scaled_f(label_pts) - 10, scaled_fft(label_pts), ["ch0(-f0)", "ch0(-f1)", "ch1(-f0)", "ch1(-f1)"], ...
    FontWeight="Bold", Rotation=60)
xlim([-500, 500])
title("FFT of a Frequency-Modulated Carrier Wave (2.4GHz)");
legend("Modulated Signal", "Carrier","","");
ylabel("Amplitude");
xlabel("Offset from Carrier Frequency (kHz)");
xticks(-500:50:500);