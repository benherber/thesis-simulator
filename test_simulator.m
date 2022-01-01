clc; clear;

channel = @(given) awgn(given, 10, "measured", "linear");

tags = [[1;1;1],[1;-1;1]];
tag_modes = ["lo", "hi"];

sim = Simulator(tags, tag_modes, channel, NUM_SYMBS);