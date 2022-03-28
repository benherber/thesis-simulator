function [params, total_symbs] = BerPlotterConstants(num_symbs, modulation_order)
    arguments
        num_symbs double {mustBePositive} = 1;
        modulation_order double {mustBePositive} = 2;
    end

    total_symbs = 600;

    params = bherber.thesis.SimulationConstants( ...
                2.4e9, ...   "Fc"
                100e3, ...   "Fb_base"
                100e3, ...   "Fb_step"
                1 / 2 * (1 / (4 * 2.4e9)), ...   "Fb_channel_spacing"
                5, ...       "num_channels"
                4, ...       "fs_granularity"
                num_symbs, ... "num_symbs"
                1e5, ...     "symb_freq"
                1, ...       "sim_sym_ratio"
                sqrt(8), ... "amplitude"
                4,  ...      "num_elements"
                3, ...        "pattern_len"
                modulation_order ... "m_ary_modulation"
            );
end