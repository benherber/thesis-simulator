% -------------------------- %
% Friis Path Loss            %
% Name   : friis_path_loss.m %
% Author : Benjamin Herber   %
% Date   : Fall 2021         %
% -------------------------- %

function [res] = friis_path_loss(wave, d)
%FRIIS_PATH_LOSS of signal through free space
%   Using Friis Formula for path loss through free space, calculate the
%   remaining wave.

% Load Constants
simulation_constants

% Calculate power loss
power_loss = WAVELEN / (16 * pi * pi * d * d);

% Calculate the resultant wave
res = wave * sqrt(2) * sqrt(power_loss);
end
    
    