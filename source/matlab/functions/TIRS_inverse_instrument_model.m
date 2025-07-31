function [srd_counts, sigDN, thermDN, encoder] = TIRS_inverse_instrument_model( ...
    pd, rad_fpaths, output_fpath, output_ctime, ...
    calseq_ctime_origin, calseq_cadence, ...
    calseq_lengths, downlink_ctime_origin, downlink_cadence, downlink_length, ...
    downlink_tperturb_size, downlink_tperturb_length, ...
    varargin )
%
% Inverse instrument model for TIRS. This function takes simulated,
% high spectral resolution data (from PCRTM) and runs it
% 'backward' through the instrument model to produce simulated DN
% (Digital Number) arrays that would be measured by TIRS.
%
% IMPORTANT: the model can take multiple input granules and
% concatenates these together for ease of processing. Because of
% the size of the high spectral resolution radiance arrays, and the
% various array copies that need to be made, the memory footprint
% is on order of 20 GB per single granule-sized input rad file.
% Thus we probably can reasonable process up to ~4
% granules together on the test machine, but no more than that.
%
% This is based on 'PREFIRE_Forward_Model_6' from B.; I decided
% to break with that naming convention, since describing this as a
% "Forward" model is arguably a misnomer (and we have been
% internally calling this the 'inverse instrument model' for some
% time.)
%
% required inputs:
% pd - config data and paths from configure_toplevel_IO.
% rad_fpaths - string filepaths for high spectral resolution simulated
%     data in netCDF format. These should also include a geometry group.
%     should be a cell array of string file paths.
% output_fpath - output file path
% output_ctime - a 2 element array with the ctime range for the
%    output file. An error is raised if the requested ctime does
%    not contain any of the input data.
% calseq_ctime_origin - timestamp (in seconds ctime) for a
%    calibration sequence start. All calseq will be relative to
%    this time, separated by multiples of the cadence time.
% calseq_cadence - the repetition cycle of calseq in seconds.
% calseq_lengths - length (in seconds) of each dwell in a
%    calseq (the space view and internal cal target)
%    This is a 2 element array, so that the space and cal dwells
%    can be different lengths.
% downlink_ctime_origin - timestamp (in seconds ctime) for a
%    downlink sequence start. All downlink will be relative
%    to this time, separated by multiples of the cadence time.
% downlink_cadence - the repetition cycle of downlink in seconds.
% downlink_length - length (in seconds) of downlink; this will
%    cause a data gap of this length.
% downlink_tperturb_size - size, in deg C, of internal temperature
%    perturbation that occurs during downlink.
% downlink_tperturb_length - time, in seconds, of temperature
%    perturbation from downlink. This length is used as an
%    e-folding time; the temperture perturbation's time variation
%    is basically exp(-ctime/tperturb_length)*tperturb_size, which
%    occurs at the end of the downlink gap.
%
% Options
% add_noise - whether to add noise (1 or true) or not (any other value)
% rng_seed - to force a rng seed (any integer) or not (False)
% noise_scale - apply an overall multiplicative scaling to the
%     simulated detector noise. In comparing to lab test data, the
%     internal noise values appear to be underestimates, so this
%     input allows for an ad hoc scaling to be applied to the
%     noise. Default is 1.0, which is effectively no scaling.
% write_matfile - whether to write out an additional mat format
%     file, with more internal variables. useful for debugging/testing.
% constant_temp - to force internal temperatures to the given
%     constant value, in K. useful for debugging/testing.
%
% returns:
% srd_counts: DN array, converted to telemetry order.
%     array size is (n_spectral x n_spatial, n_frames)
% sigDN: DN array, in science order.
%     array size is (n_spectral, n_spatial, n_frames)
% thermDN: DN array for thermistor measurements.
%     array size is (8, n_frame)
% encoder: encoder position array for instrument mirror, size (nframe)

ti = TIRS_info;  % Some TIRS-related parameters

% defaults
add_noise = false;
rng_seed = 'shuffle';
noise_scale = 1.0;
write_matfile = false;
apply_constant_temp = false;

if length(varargin) >= 1
    add_noise = varargin{1} == 1;
end

if length(varargin) >= 2
    if isnumeric(varargin{2})
        rng_seed = varargin{2};
    end
end

if length(varargin) >= 3
    if isnumeric(varargin{3})
        noise_scale = varargin{3};
    end
end

if length(varargin) >= 4
    write_matfile = varargin{4} == 1;
end

if length(varargin) >= 5
    if isnumeric(varargin{5})
        apply_constant_temp = true;
        constant_temp = varargin{5};
    end
end

% determine sensor num from TIRS instrument model file (has to be
% from filename at this time, should be fixed.)
[~,rad_file,~] = fileparts(rad_fpaths{1});
TIRS_num = str2num(rad_file(12));

[TIRS] = prefire_TIRS(ti.ROIC_tau, pd, 2, 1, TIRS_num);

% load the HyperSpectral radiance per wavenumber (hsrad) and
% wavenum from the input Fullspectra file.
% concatenate into single array.
sc_latitude = zeros(0);
ctime = zeros(0);
for f=1:length(rad_fpaths)
    rad_fpath = rad_fpaths{f};
    if f == 1
        hswavenum = ncread(rad_fpath, '/SimRad-FullRes/wavenum');
        hsrad_wavenum = zeros([length(hswavenum),8,0]);
    end
    hsrad_wavenum = cat(3, hsrad_wavenum, ...
                        ncread(rad_fpath, '/SimRad-FullRes/radiance'));
    sc_latitude = cat(1, sc_latitude, ...
                      ncread(rad_fpath, '/Geometry/subsat_latitude'));
    ctime = cat(1, ctime, ...
                ncread(rad_fpath, '/Geometry/ctime'));
end

% "unapply" the half integration time offset
ctime = ctime + (0.5*TIRS.tau);

% determine frame range for output.
[~,output_start_frame] = min(abs(ctime - output_ctime(1)));
[~,output_end_frame] = min(abs(ctime - output_ctime(2)));

if (output_end_frame == output_start_frame)
    disp(sprintf('ctime range: %12d %12d',round(ctime(1)), round(ctime(end))))
    disp(sprintf('requested:   %12d %12d',round(output_ctime(1)), round(output_ctime(2))))
    error('input ctime does not overlap data in input rad files')
end

% Set up instrument thermal data, need only one latitude 'strip' choose scene 4,
% this is the instrument thermal model
num_frames = length(ctime);
if apply_constant_temp
    thermistor_temps = zeros([8,num_frames]) + constant_temp;
else
    % 0.001 = 1mK here is the additional temperature variation
    thermistor_temps = TIRS_thermal_model(pd, sc_latitude, ctime, 0.001);
    thermistor_temps = thermistor_temps + 273.15;
end

% determines where the downlinks are located, and inserts a
% temperature perturbation in slot 1 (Toroid) which will impact
% the instrument thermal radiation background.
downlink_gaps = TIRS_IIM_calc_downlink_gaps( ...
    ctime, downlink_ctime_origin, downlink_cadence, downlink_length);
thermistor_temps = TIRS_IIM_downlink_tperturb( ...
    ctime, thermistor_temps, downlink_gaps, ...
    downlink_tperturb_size, downlink_tperturb_length);

fpath_to_load = fullfile( ...
    pd.ancillary_data_dir, 'BRF', ...
    sprintf('PREFIRE_BRFM_v0.12_TIRSN%1d_2023-10-05.mat', TIRS_num) );
BRF = load(fpath_to_load);

instr_model_fpath = fullfile( ...
    pd.ancillary_data_dir, 'SRF', ...
    sprintf('PREFIRE_TIRS%1d_SRF_v12_2024-04-18.nc', TIRS_num) );
SRF_wavelen = ncread(instr_model_fpath, '/wavelen');
SRF = ncread(instr_model_fpath, '/srf');
SRF_ncid = netcdf.open(instr_model_fpath, 'NC_NOWRITE');
SRF_im_version = netcdf.getAtt(SRF_ncid, netcdf.getConstant('NC_GLOBAL'), ...
                               'instrument_model_version');
netcdf.close(SRF_ncid);
SRF_delta_wl = SRF_wavelen(2) - SRF_wavelen(1);

% apply the wavelength delta here, so that dot products are
% effectively the discrete integrals.
SRF = SRF * SRF_delta_wl;

% prepare local variables -
% computing the channel spectrally-integrated radiance
ch_irad = zeros(TIRS.spectral_channels,TIRS.spatial_scenes,num_frames);

% compute per-wavenum conversion factor to convert the simulated
% full spectra to per-wavelength.
% conversion includes mW->W, and cm->um, resulting in factor of 10.
hswavelen = 1e4 ./ hswavenum;
wave_conversion = 10.0 ./ hswavelen.^2;

% using broadcasting here. L1b.full_spec_rad is (n_wave, n_xtrack, n_atrack)
hsrad_wavelen = hsrad_wavenum .* wave_conversion;

% now, interpolate the full_spec to the wavelen grid for SRF.
% this seems to do the right thing, broadcasting over axis 2 and 3
hsrad = interp1(hswavelen, hsrad_wavelen, SRF_wavelen, 'linear', 0.0);

%ideally this could be completely vectorized
%this results in a 63x8 array - the SRF array does not have the
%lead (undispersed) channel.
SRFsums = squeeze(sum(SRF, 3))';
for j=1:TIRS.spatial_scenes
    localSpec = squeeze(hsrad(:,j,:));
    for k=2:TIRS.spectral_channels
        % TBD - replace SRFsums check with detector mask
        if SRFsums(k-1,j) > 0
            localSRF = squeeze(SRF(j,k-1,:));
            ch_irad(k,j,:) = localSRF' * localSpec;
        end
    end % for k, loop over TIRS channels
end %for j, loop over spatial scenes


%B. add zero order sum 11/22/21
%overwrite zero order values
%zeroth order space and calibration values are sums of all channels
%(not quite true with SRF gaps at blocking between filters)
%there is no grating efficiency or reflection coefficient applied!
%AM's comment: since we are summing the "sig" computed above, this
%really just means we added those spectral radiances together
%(which is probably wrong) and implicitly assumed the order=-1
%grating efficiency (which is also probably wrong)
% TBD: should we just remove this?
for j=1:TIRS.spatial_scenes
    ch_irad(1,j,:) = sum(ch_irad(2:TIRS.spectral_channels,j,:));
end

% for now, assuming 'instrument housing' is "telescope" from the
% thermistor eng data.
T_cal = mean(thermistor_temps(4:6,:),1);
T_instr = thermistor_temps(2,:);

% Note, I implemented the helper functions to work on row vectors;
% so, need to transpose ctime.
calseqs = TIRS_IIM_prep_calseqs( ...
    ctime', calseq_ctime_origin, calseq_cadence, calseq_lengths);
calseqs = TIRS_IIM_create_calseq_chunks( ...
    calseqs, ch_irad, ctime', BRF, T_cal, T_instr);

% transfer calseq.irad_chunks into ch_irad here.
for n = 1:calseqs.num
    a = calseqs.start_frames(n) + 1;
    b = calseqs.end_frames(n);
    ch_irad(:,:,a:b) = calseqs.irad_chunks{n};
end

% for downlink gaps: overwrite the scene irradiance (ch_irad) with
% view of the internal calibrator target.
[num_downlinks, ~] = size(downlink_gaps);
for n=1:num_downlinks
    for i = downlink_gaps(n,1):downlink_gaps(n,2)
        dT = T_cal(i) - BRF.T_cal(1);
        fraci = mod(dT, 1);
        i_cT = fix(dT)+2;
        ch_irad(:,:,i) = ...
            (1-fraci) * BRF.scene_irad_c(:,:,i_cT-1) + ...
               fraci  * BRF.scene_irad_c(:,:,i_cT);
    end
end


encoder = TIRS_IIM_compute_calseq_encoder( ...
    ctime, calseqs, TIRS.space_enc, TIRS.cal_enc, TIRS.nadir_enc);
view_flag = zeros(num_frames,1);
view_flag(encoder == TIRS.nadir_enc) = 3;
view_flag(encoder == TIRS.space_enc) = 1;
view_flag(encoder == TIRS.cal_enc) = 2;

% Now, derive the background irradiance terms.

% linear interpolations for each temperature
dT_det = thermistor_temps(7,:) - BRF.T_cal(1); % FPA center/back
dT_flt = thermistor_temps(3,:) - BRF.T_cal(1); % Filter assembly
dT_tor = thermistor_temps(1,:) - BRF.T_cal(1); % toroidal mirror (optics)
fraci_det = mod(dT_det,1);
fraci_flt = mod(dT_flt,1);
fraci_tor = mod(dT_tor,1);
i_det = fix(dT_det)+2;
i_flt = fix(dT_flt)+2;
i_tor = fix(dT_tor)+2;

backgr_det = (1-fraci_det) .* BRF.BRFMD(:,i_det-1) + fraci_det .* BRF.BRFMD(:,i_det);
backgr_flt = (1-fraci_flt) .* BRF.BRFMF(:,i_flt-1) + fraci_flt .* BRF.BRFMF(:,i_flt);
backgr_tor = (1-fraci_tor) .* BRF.BRFMT(:,i_tor-1) + fraci_tor .* BRF.BRFMT(:,i_tor);
% "scene" background also uses toroid temperature
backgr_sce = (1-fraci_tor) .* BRF.BRFMS(:,i_tor-1) + fraci_tor .* BRF.BRFMS(:,i_tor);

% now compute the net irradiance time series per channel/scene, by
% multiplying the observed radiance (ch_irad, [W/(m2 sr)]) by the geometric factor
% (in sr)....
ch_irrad = ch_irad * BRF.scene_geom_fac;
% and build the per-channel backgrounds from each term.
% these are filter dependent so we need to map these from filter to
% channel. Leave them split into components first, for
% vis/debugging purposes.
ch_backgr_det = zeros(size(ch_irrad));
ch_backgr_flt = zeros(size(ch_irrad));
ch_backgr_tor = zeros(size(ch_irrad));
ch_backgr_sce = zeros(size(ch_irrad));

for k=1:TIRS.spectral_channels
    f = TIRS.filter_channel(k);
    if (f > 0)
        for j=1:TIRS.spatial_scenes
            ch_backgr_det(k,j,:) = backgr_det(f,:);
            ch_backgr_flt(k,j,:) = backgr_flt(f,:);
            ch_backgr_tor(k,j,:) = backgr_tor(f,:);
            ch_backgr_sce(k,j,:) = backgr_sce(f,:);
        end
    end
end

% combine into total net irradiance
% (could be negative if the detector is warm and the scene is cold)
ch_background_irrad = ch_backgr_det + ch_backgr_flt + ch_backgr_tor + ch_backgr_sce;
total_irrad = ch_irrad + ch_background_irrad;

% convert from irrad to volts by multiplying the detector area (to
% get power in watts) and the responsivity (V/W)
d = TIRS.d/1e6; % convert from microns to meters
for m=1:num_frames;
    sigV(:,:,m) = TIRS.R.*total_irrad(:,:,m)*d^2;
end

% create random noise vector for signal, in units of Volts.
if add_noise
    if isnumeric(rng_seed)
        rng(rng_seed)
    elseif ~strcmp(rng_seed,'none')
        rng(rng_seed)
    end
    r = random(makedist('Normal'),[TIRS.spectral_channels,TIRS.spatial_scenes,num_frames]);
    sigmaV = r*TIRS.noise;
else
    sigmaV = zeros([TIRS.spectral_channels, TIRS.spatial_scenes, num_frames]);
end

% apply noise scaling
sigmaV = sigmaV * noise_scale;

% gain reverses polarity on opposite channel; 'phase shift' between
% scenes 1-4 and 5-8.
% Note this looks different than the similar ("inverse")
% calculation in L1A because this is done on the science ordering,
% but the polarity shift is applied in L1A to the ROIC ordered data.
sigDN = zeros(size(sigV));
s1 = 1;
s2 = TIRS.spatial_scenes/2;
for k=1:TIRS.spectral_channels
    sigDN(k,s1:s2,:) = TIRS.o + ((-1)^k)*TIRS.g * (sigV(k,s1:s2,:)+sigmaV(k,s1:s2,:));
end
s1 = TIRS.spatial_scenes/2 + 1;
s2 = TIRS.spatial_scenes;
for k=1:TIRS.spectral_channels
    sigDN(k,s1:s2,:) = TIRS.o + ((-1)^(k+1))*TIRS.g * (sigV(k,s1:s2,:)+sigmaV(k,s1:s2,:));
end

sigDNfloat = sigDN;
sigDN = round(sigDN,0);

min_sigDN = min(sigDN,[],'all');
max_sigDN = max(sigDN,[],'all');
if min_sigDN < 0
    disp('signal DN has values below zero');
end
if max_sigDN >= 2^16
    disp('signal DN has values at or above 2^16')
end

% force to 0 - 2^16 (uint 16 range) via modulo - I assume this is
% what would happen in the payload firmware.
sigDN = mod(sigDN, 2^16);

% Re-order data to assumed readout format
sigDN_ROIC = ROIC_map_oi([TIRS.spatial_scenes, num_frames],TIRS.spectral_channels,sigDN);

% flattening to (nframes, 512)
srd_counts = zeros([TIRS.spectral_channels * TIRS.spatial_scenes, num_frames]);
for j=1:TIRS.spatial_scenes
    for k=1:TIRS.spectral_channels
        srd_counts(8*(k-1)+j,:) = sigDN_ROIC(k,j,:);
    end
end

% convert thermistor temps to DN
thermDN = inverse_tempest_temperature(thermistor_temps);

% new transformations for updated L0 format. This looks a bit
% strange, because we are undoing some steps just done immediately
% above; however, we want to output prepped L0 payload data, and in
% the current processing flow, this file includes some
% transformations of the data (conversion of engineering DN to
% thermistor temps; and translation of "SRD" raw counts to "SCI"
% raw counts.
adj_sci_counts = srd_to_sci_DN(TIRS, srd_counts);

% clip final data according to output ctime range (we checked it at
% the start to make sure it covered the available times.)
a = output_start_frame;
b = output_end_frame;
ctime = ctime(a:b);
encoder = encoder(a:b);
adj_sci_counts = adj_sci_counts(:,:,a:b);
thermDN = thermDN(:,a:b);
thermistor_temps = thermistor_temps(:,a:b);

% load needed global attributes from input file
% only creating the minimal attributes here - right now this is
% just the input file, orbit_sim_version, SRF version, sensor_ID.
[rad_path, rad_file, rad_file_ext] = fileparts(rad_fpath);
rad_file = strcat(rad_file, rad_file_ext);

global_atts = struct;
global_atts.input_product_files = rad_file;

ncid = netcdf.open(rad_fpath, "NC_NOWRITE");
global_atts.sensor_ID = netcdf.getAtt( ...
    ncid, netcdf.getConstant('NC_GLOBAL'), 'sensor_ID');
global_atts.orbit_sim_version = netcdf.getAtt( ...
    ncid, netcdf.getConstant('NC_GLOBAL'), 'orbit_sim_version');
netcdf.close(ncid)

global_atts.SRF_NEdR_version = SRF_im_version;

write_prepped_L0_payload( ...
    pd.L0_anc_data_dir, adj_sci_counts, encoder, ...
    thermDN, thermistor_temps, ctime, ...
    output_fpath, global_atts);

% for debugging
if write_matfile
    save(strcat(output_fpath, '.mat'), ...
         'thermistor_temps', 'calseqs', ...
         'srd_counts', 'view_flag', 'encoder', 'thermDN', 'ctime', ...
         'ch_irrad', 'ch_background_irrad', ...
         'ch_backgr_det', 'ch_backgr_flt', 'ch_backgr_tor', 'ch_backgr_sce', ...
         'sigV', 'sigmaV', 'sigDN', 'sigDNfloat');
end

end
