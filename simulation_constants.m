FC = 2.4e9; % Carrier Frequency

FB0 = 75e3; % First bit=0 Frequency Sideband
FB1 = 100e3; % First bit=1 Frequency Sideband
FB01 = 125e3; % Second bit=0 Frequency Sideband
FB11 = 150e3; % Second bit=1 Frequency Sideband

WAVELEN = physconst("Lightspeed") / FC; % Wavelength of Carrier
GRANULARITY = 3; % Sampling Granularity
Fs = FC * GRANULARITY; % Sampling frequency

NUM_SYMBS = 10; % Total number of symbols in simulation
SYMB_FREQ = 100e3; % Frequency of symbols
SYMB_SIZE = (1 / SYMB_FREQ) * Fs; % Samples per symbol

SIM_SYMB_RATIO = 10; % Simulation steps per symbol window

SIMSTEP_FREQ = SYMB_FREQ * SIM_SYMB_RATIO; % Frequency of simulation-steps
SIMSTEP_SIZE = (1 / SIMSTEP_FREQ) * Fs; % Number of samples per simulation-step

TOTAL_TIME = NUM_SYMBS * (1 / SYMB_FREQ); % Total simulation time
TOTAL_SAMPLES = SYMB_SIZE * NUM_SYMBS; % Total simulation samples

A = 0.1; % Amplitude of Carrier
NUM_ELEMENTS = 4; % Number of elements in Van Atta Array
% TRANSMIT_POWER = 20e-3; % Typical average power of a 2.4GHz WiFi signal (dB)