clear; clc;
import bherber.thesis.Simulator bherber.thesis.TagType
addpath(genpath(".."));
filename = sprintf("data/fskbers::%s.mat", datestr(now, "mm-dd-yy-HH:MM:SS"));
save(filename);

errors = [];
try
    pool = gcp;
    num_threads = pool.NumWorkers;
    modulation_order = 2;
    [params, num_symbs] = BerPlotterConstants(1, modulation_order);
    
    symbs_per_worker = ceil(num_symbs / num_threads);
    num_symbs = symbs_per_worker * num_threads;
    bits_per_symb = log2(params.m_ary_modulation);
    
    tags = [1;0;0];
    modes = [TagType.FSK_LO];

    bers = NaN(1, num_threads);

    for snr_db = linspace(0, 18, 7)
        parfor thread = 1:num_threads
            scaled_params = BerPlotterConstants(symbs_per_worker, modulation_order);
            sim = Simulator(tags, modes, @(a) a, scaled_params, snr_db=snr_db, complex_noise=true);
            sim.auto_align()
            delay = sim.tag_delays(1);
            symbols = zeros(int32(scaled_params.simstep_sz), int32(scaled_params.sim_sym_ratio), 2);
            num_errors = 0;
    
            for jdx = 1:scaled_params.sim_sym_ratio
                    symbols(:, jdx, 2) = sim.step();
            end
    
            for idx = 1:symbs_per_worker
                symbols(:, :, 1) = symbols(:, :, 2);
                for jdx = 1:scaled_params.sim_sym_ratio
                    step = sim.step();
                    symbols(:, jdx, 2) = step;
                end
    
                curr_steps = symbols(:, :, 1);
                next_steps = symbols(:, :, 2);
                frame = [curr_steps(:); next_steps(:)];
                shifted_frame = delayseq(frame, -1 * delay);
                true_symbol = shifted_frame(1:scaled_params.symb_sz);
                demod_bits = sim.tags(1).demodulate(idx, true_symbol);
                start = ((idx - 1) * bits_per_symb) + 1;
                stop = start + bits_per_symb - 1;
                expected = sim.tags(1).bits(start:stop).';
                num_errors = num_errors + sum(expected ~= demod_bits);
            end
    
            bers(thread) = num_errors / (num_symbs * bits_per_symb);
        end
    
        savetotal = struct(replace(sprintf("totalber%dMEbNo%d", modulation_order, int32(snr_db)), ...
            "-", "NEG"), sum(bers, "omitnan"));
        save_out(filename, savetotal);
    end

catch errors
end

% delete(gcp);
if ~isempty(errors)
    rethrow(errors);
end

%%
function save_out(filename, structure)
    save(filename, "-append", "-struct", "structure");
end
