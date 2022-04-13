
params = BerPlotterConstants(6000, 8);

f1 = 100e3;

time = double(1:(100 * (params.symb_sz - 1))) / params.Fs;

s1 = square(2 * pi * f1 * time) .* params.amplitude .* cos(2 .* pi .* params.Fc * time);
f_modulation = s1;

nfft = 100000000;
f_fft = fftshift(fft(f_modulation, nfft));

% Frequency
f = params.Fs * (-nfft/2 : nfft/2 - 1) / nfft;
scaled_f = (f - 2.4e9) / 1e3;

figure();
hold on
scaled_fft = 2 * abs(f_fft) / length(f_modulation);
plot(scaled_f, scaled_fft);
carrier_line = xline(0, "r", "Carrier", DisplayName="Carrier Frequency");
carrier_line.LabelVerticalAlignment = "middle";
carrier_line.LabelHorizontalAlignment = "center";
xlim([-6000, 6000])
title("FFT of an FSK signal", FontSize=20);
ylabel("Amplitude", FontSize=18);
xlabel("Offset from Carrier Frequency (kHz)", FontSize=18);
xticks(-6000:500:6000);

%% Packed Spectrum

clear; clc;


params = BerPlotterConstants(6000, 8);

f_spacing = 0.5 / (1 / params.symb_freq);
f1 = 50e3;
num_peaks = 400;

time = double(1:(100 * (params.symb_sz - 1))) / params.Fs;

f_modulation = zeros(size(time));
for idx = 1:num_peaks
    disp(idx);
    f_modulation = f_modulation + ...
        square(2 * pi * (f1 + (f_spacing * (idx - 1))) * time) .* ...
        params.amplitude .* cos(2 .* pi .* params.Fc * time);
end

%%

nfft = 100000000;
f_fft = fftshift(fft(f_modulation, nfft));

% Frequency
f = params.Fs * (-nfft/2 : nfft/2 - 1) / nfft;
scaled_f = (f - 2.4e9) / 1e3;

%%
figure();
hold on
scaled_fft = 2 * abs(f_fft) / length(f_modulation);
plot(scaled_f, scaled_fft);
carrier_line = xline(0, "r", "Carrier", DisplayName="Carrier Frequency");
carrier_line.LabelVerticalAlignment = "middle";
carrier_line.LabelHorizontalAlignment = "center";
xlim([-20000, 20000])
title("FFT of a FDM-FSK signal", FontSize=20);
ylabel("Amplitude", FontSize=18);
xlabel("Offset from Carrier Frequency (kHz)", FontSize=18);
% xticks(-6000:500:6000);