import bherber.thesis.PPersistPacketSimulator ...
    bherber.thesis.SimulationConstants

id_len = 16;
data_len = 0;

make_params = @(num_channels) bherber.thesis.SimulationConstants( ...
                2.4e9, ...   "Fc"
                1e6, ...   "Fb_base"
                1e6, ...   "Fb_step"
                0.5 / (1 / (100e3)), ...   "Fb_channel_spacing"
                num_channels, ...       "num_channels"
                4, ...       "fs_granularity"
                (id_len + data_len), ... "num_symbs"
                1e6, ...     "symb_freq"
                1, ...       "sim_sym_ratio"
                sqrt(8), ... "amplitude"
                4,  ...      "num_elements"
                3, ...        "pattern_len"
                2, ... "m_ary_modulation"
                id_len ... "preamble_len"
            );

make_sim = @(num_packets, num_tags, type, ppersist, params) ...
    PPersistPacketSimulator(num_packets, repelem([1;0;0], 1, num_tags), type, ppersist, params);

if isempty(gcp("nocreate"))
  pc = parcluster("local");
  pc.JobStorageLocation = "/tmp/";
  if ~isempty(getenv("SLURM_ARRAY_TASK_ID"))
      parpool(pc, 12);
  else
      parpool(pc, 6)
  end
end
pool = gcp;
num_threads = pool.NumWorkers;