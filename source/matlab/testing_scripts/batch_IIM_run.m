% Quasi-manual runs of the inverse instrument model (reading from
% Orbit sim data) and processing the simulated L0-payload prepped
% netCDF files with the L1A algorithm.
% We do this by manually assigning the pd structure elements before
% calling produce_L1A_granule.
%
% Note that the associated prepped Bus data NC files were manually
% copied into the inputs directory (though they could be referenced
% directly to the /data/datasim locations)

addpath('../PREFIRE_tools/functions');
addpath('../PREFIRE_L0/functions');
addpath('../functions');

pd = configure_toplevel_IO;
% sorta hack for testing
pd.L0_anc_data_dir = '/home/mmm/projects/PREFIRE_L0/dist/ancillary';
pd.tools_anc_data_dir = '/home/mmm/projects/PREFIRE_tools/dist/ancillary';

% nominal cal sequences (12-minute cadence)
calseq_ctime = 662800000;
calseq_cadence = 720.0;
calseq_length = 6 * 0.7;

% 5-minute long downlinks, every 12 h.
downlink_ctime = calseq_ctime + 7200;
downlink_cadence = 12 * 3600;
downlink_length = 600;
downlink_tperturb_size = 0.25;
downlink_tperturb_length = 60.0;

add_noise = true;
rng_seed = 1237;
noise_scale = 2.2;
write_matfile = false;

input_dir = '/data/datasim/orbitsim_03/by_kind/ANC-SimRad-FullRes/';
output_dir = '/data/users/mmm/TIRS_IIM_output/';

bus_tlm_dir = '/data/datasim/orbitsim_03/by_kind/0-PreppedBusTlm/';
prefire1_bus_tlm_files = dir([bus_tlm_dir 'prefire_01_bus_tlm*.nc']);
prefire2_bus_tlm_files = dir([bus_tlm_dir 'prefire_02_bus_tlm*.nc']);
% convert from structs to files
prefire1_bus_tlm_files = {prefire1_bus_tlm_files.name};
prefire2_bus_tlm_files = {prefire2_bus_tlm_files.name};

%%%%%%%%%%%%%%%%%%%%%%%
% Strings to construct for eventual PROC_MODE=3 L1A runs
%
% GRANULE_START_ID_END:
% ctime1_granuleID_ctime2
% where ctime1, stime2 are the start/end ctimes
%
% L0_PAYLOAD_FPATHS - one string containing N Payload files,
% delimited by ||
% L0_BUS_FPATHS - same as PAYLOAD bus for the bus files.

%%%%%%%
% TIRS1
TIRS1_simfiles = sort_file_structs(dir([input_dir 'PREFIRE_SAT1_*_000*.nc']));
TIRS1_L0_files = cell(0);

% create L0 payload granules ...
for n=1:length(TIRS1_simfiles)

    rad_f = TIRS1_simfiles(n).name;
    rad_fpath = [TIRS1_simfiles(n).folder filesep rad_f];
    L0_output_f = strrep(rad_f, 'PREFIRE_SAT1_ANC-SimRad-FullRes_B03_R00_', 'prefire_01_simpayload_tlm_');
    L0_output_fpath = [output_dir filesep L0_output_f];

    tic
    TIRS_inverse_instrument_model( ...
        pd, rad_fpath, L0_output_fpath, ...
        calseq_ctime, calseq_cadence, calseq_length, ...
        downlink_ctime, downlink_cadence, downlink_length, ...
        downlink_tperturb_size, downlink_tperturb_length, ...
        add_noise, rng_seed, noise_scale, write_matfile);

    TIRS1_L0_files{n} = L0_output_fpath;
    disp(['Created L0 payload file ' L0_output_f ' etime: ' num2str(toc)]);

end

% ... and collect the strings for driving L1a
TIRS1_granule_strings = cell(0);
TIRS1_payload_fpaths = cell(0);
TIRS1_bus_fpaths = cell(0);

for n=2:length(TIRS1_simfiles)-1
    L0_file = TIRS1_L0_files{n};
    rad_f = TIRS1_simfiles(n).name;
    rad_fpath = [TIRS1_simfiles(n).folder filesep rad_f];
    ctime = ncread(rad_fpath, '/Geometry/ctime');
    utc_tuples = ncread(rad_fpath, '/Geometry/time_UTC_values');

    granule_ID = L0_file(end-7:end-3);
    TIRS1_granule_strings{n-1} = ...
        sprintf('%15.5f_%s_%15.5f', ctime(1), granule_ID, ctime(end));
    TIRS1_payload_fpaths{n-1} = [ ...
        TIRS1_L0_files{n-1} '||' ...
        TIRS1_L0_files{n}   '||' ...
        TIRS1_L0_files{n+1} ];

    utc_ts1 = (    utc_tuples(6,1) + 100*utc_tuples(5,1) + ...
               1e4*utc_tuples(4,1) + 1e6*utc_tuples(3,1) + ...
               1e8*utc_tuples(2,1) +1e10*utc_tuples(1,1) );
    utc_ts2 = (    utc_tuples(6,end) + 100*utc_tuples(5,end) + ...
               1e4*utc_tuples(4,end) + 1e6*utc_tuples(3,end) + ...
               1e8*utc_tuples(2,end) +1e10*utc_tuples(1,end) );
    bus_tlm_files = select_bus_tlm_files(utc_ts1, utc_ts2, ...
                                         prefire2_bus_tlm_files);
    bus_fpaths = [bus_tlm_dir filesep bus_tlm_files{1}];
    for j=2:length(bus_tlm_files)
        bus_fpaths = [bus_fpaths '||' bus_tlm_dir filesep ...
                      bus_tlm_files{j}];
    end
    TIRS1_bus_fpaths{n-1} = bus_fpaths;

end

%%%%%%%
% repeat for TIRS2
TIRS2_simfiles = dir([input_dir 'PREFIRE_SAT2_*_000*.nc']);
TIRS2_L0_files = cell(0);

% create L0 payload granules ...
for n=1:length(TIRS2_simfiles)
        
    rad_f = TIRS2_simfiles(n).name;
    rad_fpath = [TIRS2_simfiles(n).folder filesep rad_f];
    L0_output_f = strrep(rad_f, 'PREFIRE_SAT2_ANC-SimRad-FullRes_B03_R00_', 'prefire_02_simpayload_tlm_');
    L0_output_fpath = [output_dir filesep L0_output_f];

    tic;
    TIRS_inverse_instrument_model( ...
        pd, rad_fpath, L0_output_fpath, ...
        calseq_ctime, calseq_cadence, calseq_length, ...
        downlink_ctime, downlink_cadence, downlink_length, ...
        downlink_tperturb_size, downlink_tperturb_length, ...
        add_noise, rng_seed, noise_scale, write_matfile);

    TIRS2_L0_files{n} = L0_output_fpath;
    disp(['Created L0 payload file ' L0_output_f ' etime: ' num2str(toc)]);

end

% ... and collect the strings for driving L1a
TIRS2_granule_strings = cell(0);
TIRS2_payload_fpaths = cell(0);
TIRS2_bus_fpaths = cell(0);

for n=2:length(TIRS2_simfiles)-1

    L0_file = TIRS2_L0_files{n};
    rad_f = TIRS2_simfiles(n).name;
    rad_fpath = [TIRS2_simfiles(n).folder filesep rad_f];
    ctime = ncread(rad_fpath, '/Geometry/ctime');
    utc_tuples = ncread(rad_fpath, '/Geometry/time_UTC_values');

    granule_ID = L0_file(end-7:end-3);
    TIRS2_granule_strings{n-1} = ...
        sprintf('%15.5f_%s_%15.5f', ctime(1), granule_ID, ctime(end));
    TIRS2_payload_fpaths{n-1} = [ ...
        TIRS2_L0_files{n-1} '||' ...
        TIRS2_L0_files{n}   '||' ...
        TIRS2_L0_files{n+1} ];

    utc_ts1 = (    utc_tuples(6,1) + 100*utc_tuples(5,1) + ...
               1e4*utc_tuples(4,1) + 1e6*utc_tuples(3,1) + ...
               1e8*utc_tuples(2,1) +1e10*utc_tuples(1,1) );
    utc_ts2 = (    utc_tuples(6,end) + 100*utc_tuples(5,end) + ...
               1e4*utc_tuples(4,end) + 1e6*utc_tuples(3,end) + ...
               1e8*utc_tuples(2,end) +1e10*utc_tuples(1,end) );
    bus_tlm_files = select_bus_tlm_files(utc_ts1, utc_ts2, ...
                                         prefire2_bus_tlm_files);
    bus_fpaths = [bus_tlm_dir filesep bus_tlm_files{1}];
    for j=2:length(bus_tlm_files)
        bus_fpaths = [bus_fpaths '||' bus_tlm_dir filesep ...
                      bus_tlm_files{j}];
    end
    TIRS2_bus_fpaths{n-1} = bus_fpaths;

end

save('batch_IIM_run_L1A_input_strings.mat', ...
     'TIRS1_granule_strings', 'TIRS1_payload_fpaths', 'TIRS1_bus_fpaths',...
     'TIRS2_granule_strings', 'TIRS2_payload_fpaths', 'TIRS2_bus_fpaths');


function sorted_dir_results = sort_file_structs(dir_results)
% I don't think the return from dir() is sorted;
% this helper makes sure it is sorted.
filenames = {dir_results.name};
[~,sort_idx] = sort(filenames);
sorted_dir_results = dir_results(sort_idx);
end


function bus_tlm_files = select_bus_tlm_files(utc1, utc2, available_bus_tlm_files)
% given a cell array of bus tlm files, assumed to have name pattern:
% prefire_0N_bus_tlm_startUTC_stopUTC_procUTC.nc
% return the list of files that covers the time range specified
% by utc1, utc2. The return is also a cell array, and can have zero
% elements, if no bus tlm files overlap.

nfiles = length(available_bus_tlm_files);
selected_n1 = 0;
selected_n2 = 0;

for n=1:nfiles
    tokens = strsplit(available_bus_tlm_files{n}, '_');
    tlm_utc_end = str2num(tokens{6});
    if tlm_utc_end > utc1
        selected_n1 = n;
        break
    end
end

bus_tlm_files = cell(0);

if selected_n1 == 0
    return
end

selected_n2 = selected_n1;
for n = selected_n1+1:nfiles
    tokens = strsplit(available_bus_tlm_files{n},'_');
    tlm_utc_start = str2num(tokens{5});
    disp(tlm_utc_start);
    if tlm_utc_start > utc2
        break
    end
    selected_n2 = n;
end

bus_tlm_files = available_bus_tlm_files(selected_n1:selected_n2);

end
