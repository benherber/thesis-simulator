classdef TagType
    %TAG_TYPE of a backscatter tag.
    %   Type of backscatter tag in a passive wireless sensor network.

    properties (Constant)
        FSK_LO = 0;
        FSK_HI = 1;
        FSK_RANDOM = 2;
        OOK = 3;
    end
end