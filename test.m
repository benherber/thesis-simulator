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