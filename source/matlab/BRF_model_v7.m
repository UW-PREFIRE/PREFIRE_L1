%% PREFIRE BRF Model
% 1/6/23
% normalize SRFs to fix scaling issue with real data that
% manifests through scaling of pre-calculated SRFMS and SRFMC values
%
% also will deprecate summing over grading orders and use the pre-summed
% version already available (SRF instead of SRFbyOrder)
% 
% also include arb correction to SRF, this will require commensurate change
% in L1_v1
%
% increase fidelity of SRFMC 11/13/22
% radcal hotops up to 350 K, precision better than 0.1 
% new precalculated emission is 200-400 K 0.1 K steps
% first iteration of instrument 2 model 11/3/22
% have tweaked instrument 1 SRFs for stray light near 7 um 10/18/22
%
% need to update filters 10/2/2022, this will change forward model but,
% SRF can/is updated separately...
%
% 9/30/22 running with instrument model adjusted to match spectral cal data
%
% March 8, 2021, update SRF input to PREFIRE_SRF_v0.10.4_360_2021-02-22_bw53p16.mat
% version six will bring back background calculations for each filter...
%
% rerun 5/26/20 with pi factor applied to emission read from file
% B. debugged output of v5 on 5/5/2020, units into exitance were wrong and
% the grating orders were not sorting correctly
%
% v5 development started 4/16/2020, will incorporate 2-D SRFs and create
% output for array of instrument temperatures, in principle this array need
% not provide a correct 'baseline' and only need reproduce the space and
% calibration scenes, the 'baseline' information calculated in earlier BRF
% models %should% not be required for L1A processing, but it could be
% important for understanding instrument dynamic range and saturation or
% gain non-linearity
%
% v5 strip out the fine detail calculation and just use SRF values to
% create the new table
%
% v4 implemented by B. 11/28/2019, adapted from SRF code by A.

% dependent on exitance.m that was annoted by B. once units were confirmed

% thermal model data stored thermal_model_data(); 
% orbital simulation thermal data from Ian Mckinley 11/20/2019, worst case based on SDL spacecraft model
% files names rev4_hot_case_transient.xlsm and rev4_cold_case_transient.xlsm

% B. updated for new filter transmission and reflection data 2/4/2020

[script_path, ~, ~] = fileparts(which('BRF_model_v7'));
orig_path = path;
path(orig_path, fullfile(script_path, 'functions'));

pd = configure_toplevel_IO;  % Get/set some shared filepaths and directories

% set up wavelength arrays, removed in v4 to extract from SRF load...
%(1/100 of TIRS sampling)
%spectral_channels=63; %Number of total channels is 64, the dispersive order falls across 63 of these channels
%bandwidth=54;  %(microns)
%channel_sampling = bandwidth/spectral_channels;
%channel_center_wavelens = (1:spectral_channels)'*channel_sampling;
%load('PREFIRE_SRF_v0.10.4_360_2021-02-22_bw53p16.mat');
%load('PREFIRE_SRF_v0.10.1_360_2020-04-19_ideal.mat');

%last used model file
%load('C:\Users\yyy\Documents\My Work\PREFIRE\Model\PREFIRE_SRF_v0.11_360_2022-11-15.mat');TIRSNUM = 0; % ideal model
%last used instrument 1 file
%load('Z:\instrument 1\model\PREFIRE_SRF_v0.11_360_2022-10-03_mod.mat'); TIRSNUM = 1; % use adjusted model for instrument #1, including grating efficiency update 10/3/2022
% note PREFIRE_SRF_v0.11_360_2022-11-03.mat includes 'tweaks' and has not been rerun through BRF yet!
%last used instrument 2 file
%load('C:\Users\yyy\Documents\My Work\PREFIRE\Model\PREFIRE_SRF_v0.11_360_2022-11-15_instrument2.mat');

% this file must contain:
% 'T_cal', 'T_space', 'space_emission' and 'calibrator emission' arrays
% it has no instrument dependence. Note that it contains a wavelen
% grid, that must match what is in the SRF file. (we check this.)
edat = load([pd.instrument_model_dir, ...
             '/PREFIRE_Emission_preCalculated_v0.3.0_2023-04-24.mat']);
% These variables are eventuially saved in MATLAB mat files, so
% that means we need to copy them out of the structure.
T_cal = edat.T_cal;
T_space = edat.T_space;

TIRSNUM = 1; % use adjusted model for instrument #1, including grating efficiency update 10/3/2022

% this file must contain 'SRF'
load([pd.instrument_model_dir, ...
      sprintf('/instrument_%1d_NETD_NEDR_SRF.mat', TIRSNUM)]);

num_spectral = size(SRF,2);
num_spatial = size(SRF,3);

% TIRS structure loaded for the 'arb' field.
[TIRS] = prefire_TIRS(pd,2,1,TIRSNUM);

% these are idealized center and edge wavelengths (Does not appear
% to be used anymore.)
%channel_center_wavelens = (channel_wavelens(:,1)+channel_wavelens(:,2))/2;

num_temps = size(T_cal,2);
SRFMS = zeros([num_spectral,num_spatial]);
SRFMC = zeros([num_spectral,num_spatial,num_temps]);

SRFsums = squeeze(sum(SRF, 1)); %this results in a 64x8 array
for k=1:num_spectral
    for j=1:num_spatial
        if SRFsums(k,j) > 0
            localSRF = SRF(:,k,j)/TIRS.arb(k,j);
            SRFMS(k,j) = reshape(sum(localSRF.*edat.space_emission,1),[1, 1]) / SRFsums(k,j);
            SRFMC(k,j,:) = reshape(sum(localSRF.*edat.calibrator_emission,1),[num_temps,1]) / SRFsums(k,j);
        end
    end % for j
end %for k

%Dumping cal_emission and space_emission means that total signal calculator will
%never need to call exitance.m

% Now, we need wavelength information. It seems like we
% only need the 'dd' parameter to get the high res wavelen grid.
% This assumes that the wavelen grid we recompute here (as part of
% 'load_efficiency_data' is equal to the one used in the
% precalculated Emission arrays. Assert this is true.
[~, wl] = geometry_v3(TIRSNUM);

%filter emission and reflection calculations start with high-res grid
[delta_wave, start_wave, stop_wave] = get_high_res(wl.dd);
[wavelen, grating, filter, filter_ir, filter_er] = ...
    load_efficiency_data_v4(pd, start_wave, delta_wave, stop_wave); %B. changed to v4 11/15/22...

assert(all(wavelen == edat.wavelen));

num_wave = size(wavelen, 1);
num_filters = size(filter,2);
[materials] = get_materials_properties();
BRFMD = zeros([num_filters, num_temps]); %filter reflections and detector emissions
BRFMF = zeros([num_filters, num_temps]); %filter emissivities
BRFMT = zeros([num_filters, num_temps]);
for i=1:num_wave
    for n=1:num_temps
        for f=1:num_filters
            %filter exit reflection (at detector temperatures), will include detector temperature by subtraction
            BRFMD(f,n) = (filter_er(i,f)*edat.calibrator_emission(i,n)/num_wave)+BRFMD(f,n); %precalculated emission is for any blackbody
            %filter emissivity (at filter temperature)
            BRFMF(f,n) = ((1-filter(i,f)-filter_er(i,f))*edat.calibrator_emission(i,n)/num_wave)+BRFMF(f,n);
            %instrument emission (at toroid temperature) is transmitted through filter
            BRFMT(f,n) = (materials.e_MB*filter(i,f)*edat.calibrator_emission(i,n)/num_wave)+BRFMT(f,n); %precalculated emission is for any blackbody
        end
        BRFMD(:,n) = BRFMD(:,n) - (materials.e_GB*edat.calibrator_emission(i,n)/num_wave); % this effectively wraps up the detector temperature into filter reflection
    end
end

%% get the git hash from command line; this will be stored in the
%% MATFILE; we set the filename here with a manual version number &
%% date.

output_name = [ ...
    pd.ancillary_data_dir, ...
    '/BRF/PREFIRE_BRF_v0.07_', ...
    datestr(datetime('now'), 'yyyy-mm-dd'), ...
    sprintf('_instrument%1d.mat', TIRSNUM)];

[stat, git_short_hash] = system('git rev-parse --short HEAD');
% removes the trailing return character
git_short_hash = strip(git_short_hash);

save(output_name, ...
         'wl', 'T_space', 'T_cal', 'SRFMS', 'SRFMC', ...
         'BRFMD', 'BRFMF', 'BRFMT', 'git_short_hash');
