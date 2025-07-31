% Quasi-manual runs of the inverse instrument model (reading from
% Orbit sim data) and processing the simulated L0-payload prepped
% netCDF files with the L1A algorithm.
% We do this by manually assigning the pd structure elements before
% calling produce_L1A_granule.
%
% Note that the associated prepped Bus data NC files were manually
% copied into the inputs directory (though they could be referenced
% directly to the /data/datasim locations)

addpath('../functions');
addpath('../PREFIRE_tools/functions');
addpath('../PREFIRE_L0/functions');

pd = configure_toplevel_IO;

calseq_ct = 662800000;
calseq_cadence = 675.0;
calseq_frame_count = 7;
noise_scale = 2.2;
constT = 290.0;

% TIRS1
rad_fpath = '/data/datasim/orbitsim_02/by_kind/ANC-SimRad-FullRes/PREFIRE_SAT1_ANC-SimRad-FullRes_B02_R00_20210101081841_00010.nc';
output_fpath = fullfile(pd.top_path, 'test', 'inputs', 'prefire_01_simpayload_tlm_ConstTwithnoiseNscale_20210101081841_00010.nc');
dummy = TIRS_inverse_instrument_model( ...
    pd, rad_fpath, output_fpath, calseq_ct, calseq_cadence, ...
    calseq_frame_count, true, false, noise_scale, constT);

rad_fpath = '/data/datasim/orbitsim_02/by_kind/ANC-SimRad-FullRes/PREFIRE_SAT1_ANC-SimRad-FullRes_B02_R00_20210101095410_00011.nc';
output_fpath = fullfile(pd.top_path, 'test', 'inputs', 'prefire_01_simpayload_tlm_ConstTwithnoiseNscale_20210101095410_00011.nc');
dummy = TIRS_inverse_instrument_model( ...
    pd, rad_fpath, output_fpath, calseq_ct, calseq_cadence, ...
    calseq_frame_count, true, false, noise_scale, constT);

rad_fpath = '/data/datasim/orbitsim_02/by_kind/ANC-SimRad-FullRes/PREFIRE_SAT1_ANC-SimRad-FullRes_B02_R00_20210101112938_00012.nc';
output_fpath = fullfile(pd.top_path, 'test', 'inputs', 'prefire_01_simpayload_tlm_ConstTwithnoiseNscale_20210101112938_00012.nc');
dummy = TIRS_inverse_instrument_model( ...
    pd, rad_fpath, output_fpath, calseq_ct, calseq_cadence, ...
    calseq_frame_count, true, false, noise_scale, constT);

pd.bus_L0_fpaths = {
    fullfile(pd.top_path, 'test', 'inputs', 'prefire_01_bus_tlm_20210101064635_20210101125654_20230902040053.nc'), ...
    fullfile(pd.top_path, 'test', 'inputs', 'prefire_01_bus_tlm_20210101125655_20210101192543_20230902040100.nc'), ...
                   };
pd.payload_L0_fpaths = { ...
    fullfile(pd.top_path, 'test', 'inputs', 'prefire_01_simpayload_tlm_ConstTwithnoiseNscale_20210101081841_00010.nc'), ...
    fullfile(pd.top_path, 'test', 'inputs', 'prefire_01_simpayload_tlm_ConstTwithnoiseNscale_20210101095410_00011.nc'), ...
    fullfile(pd.top_path, 'test', 'inputs', 'prefire_01_simpayload_tlm_ConstTwithnoiseNscale_20210101112938_00012.nc') ...
                   };
% these are the actual ctimes from granule 11, expanded by 0.25 second
pd.granule_beg_ts = 662810054.82965;
pd.granule_end_ts = 662815782.82065;
pd.make_L1A_aux_data = 1;
pd.granule_ID = '00011';
pd.output_dir = fullfile(pd.top_path, 'test', 'outputs', 'm3');

produce_L1A_granule(pd);

