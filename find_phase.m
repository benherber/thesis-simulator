% -------------------------- %
% Find Phase in Phased Array %
% Name   : find_phase.m      %
% Author : Benjamin Herber   %
% Date   : Fall 2021         %
% -------------------------- %

function phase = find_phase(theta, wavelen, spacing, n)
    % FIND_PHASE Find phase shift of a given element in a phased array.
    %   Given an impinging angle, wavelength, inter-element spacing, and
    %   element number in a phased array, calculate the phase shift.
    
        phase = ((2.0 * pi) / wavelen) * spacing * n * cos(theta);
    
    end
    