clear; clc;

simulation_constants

% Get time
time = (0:1/Fs:(TOTAL_TIME - 1/Fs));

% Carrier signal
carrier = complex(A * sin(complex(2 * pi * FC * time)));

channel = @(given) awgn(given, 10, "measured", "linear");

% Get random data signal
bits = randi([0, 1], NUM_SYMBS, 1);
data = repelem(bits, SYMB_SIZE).'; 

tag = Tag(1, 1, 1, "lo", time, carrier, data, channel);

f_modulation = [];
for idx = 1:(NUM_SYMBS * 10)
    f_modulation = [f_modulation, tag.step()];
end

%% Calculate FFT

nfft = 100000000;
f_fft = fftshift(fft(f_modulation, nfft));

% Carrier fft
carrier_fft = fftshift(fft(carrier, nfft));

% Frequency
f = Fs * (-nfft/2 : nfft/2 - 1) / nfft;
scaled_f = (f - 2.4e9) / 1e3;

%% FFT Plot

f1 = figure();
hold on
scaled_fft = 2 * abs(f_fft) / length(f_modulation);
scaled_carrier_fft = 2 * abs(carrier_fft) / length(carrier);
plot(scaled_f, scaled_fft);
% plot(scaled_f, scaled_carrier_fft);
label_pts = zeros(1, 2);
[ ~, label_pts(1) ] = min(abs(scaled_f - 75));
[ ~, label_pts(2) ] = min(abs(scaled_f - 100));
plot(scaled_f(label_pts), scaled_fft(label_pts), ".r");
text(scaled_f(label_pts) + 10, scaled_fft(label_pts), ["75kHz", "100kHz"], "FontWeight", "Bold", "HorizontalAlignment", "left")
[ ~, label_pts(1) ] = min(abs(scaled_f + 75));
[ ~, label_pts(2) ] = min(abs(scaled_f + 100));
plot(scaled_f(label_pts), scaled_fft(label_pts), ".r");
text(scaled_f(label_pts) - 10, scaled_fft(label_pts), ["-75kHz", "-100kHz"], "FontWeight", "Bold", "HorizontalAlignment", "right")
xlim([-500, 500])
title("FFT of a Frequency-Modulated Carrier Wave (2.4GHz)");
legend("Modulated Signal", "Carrier","","");
ylabel("Amplitude");
xlabel("Offset from Carrier Frequency (kHz)");
xticks(-500:50:500);