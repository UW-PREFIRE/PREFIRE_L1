function [error_ca] = produce_L1A_granule(pd)

error_ca = {'#NONE#', 0};  % Default

uc = constants_unit_conversion;  % Unit conversion constants
tc = constants_time(pd.tools_anc_data_dir);  % Time-related constants
ti = TIRS_info;  % Some TIRS-related parameters

  % If false, this is for 1A-GEOM -- otherwise it is for 1A-RAD:
process_payload_tlm = ~contains(pd.ref_payload_L0_fpaths{1}, '_prototlm_');

% Read in product file data specifications (done here to have access to
%  fill_value parameters for the calculations below):
if process_payload_tlm
   JSON_fspecs_L1A_fpath = fullfile(pd.ancillary_data_dir, ...
                                    'L1A_product_filespecs.json');
else
   JSON_fspecs_L1A_fpath = fullfile(pd.ancillary_data_dir, ...
                                    'L1A_GEOM_product_filespecs.json');
end
[L1A_schema, L1A_jdat] = init_nc_schema_from_JSON(JSON_fspecs_L1A_fpath, false);
n_UTC_parts = 7;

% Define/initialize output data structure array:
odat = struct;
odat.global_atts = struct;
odat.ProtoGeometry = struct;
odat.ProtoGeometry.group_atts = struct;
odat.ProtoGeometry.group_dims = struct;
if process_payload_tlm
   odat.NonGeoLoc_Radiance = struct;
%   odat.NonGeoLoc_Radiance.group_atts = struct;
   odat.NonGeoLoc_Radiance.group_dims = struct;
   odat.NonGeoLoc_BT = struct;
%   odat.NonGeoLoc_BT.group_atts = struct;
   odat.NonGeoLoc_BT.group_dims = struct;
   odat.NonGeoLoc_Channel_0 = struct;
%   odat.NonGeoLoc_Channel_0.group_atts = struct;
   odat.NonGeoLoc_Channel_0.group_dims = struct;
   odat.Cal_Artifacts = struct;
   odat.Cal_Artifacts.group_atts = struct;
   odat.Cal_Artifacts.group_dims = struct;
else
   odat.NonGeoLoc_FauxRad = struct;
   odat.NonGeoLoc_FauxRad.group_dims = struct;
end

% Parameters for L0-payload processing:

if process_payload_tlm
   % Read in relevant L0-payload fields:
   frame_types = {'obstgt', 'space', 'caltgt'};
   n_frame_types = length(frame_types);
   [error_ca, L0, TIRS_ID] = read_PREFIRE_payload_L0(pd, ti, frame_types);
   if (error_ca{2} ~= 0)
      fprintf(2, '%s\n', error_ca{1});
      return
   end
else
   % Read in relevant L0-protopayload fields:
   frame_types = {'obstgt'};
   n_frame_types = length(frame_types);
   [error_ca, L0, TIRS_ID] = read_PREFIRE_protopayload_L0(pd, frame_types);
   if (error_ca{2} ~= 0)
      fprintf(2, '%s\n', error_ca{1});
      return
   end
end
odat.global_atts.sensor_ID = sprintf('TIRS%02d', TIRS_ID);

% Load/calculate wavelength, effective pixel shapes, et cetera for this sensor:
[SRFMS, SRFMC, NEDR, dtbf, T_cal, gm, TIRS, ARF, sensor_d, factor] = ...
                                        load_geom_SRFM(pd, TIRS_ID, uc, ti, L0);

% load information for flagging portions of orbit due to eclipse transitions
info_file = fullfile(pd.ancillary_data_dir, ...
                     sprintf('TIRS%1d_orbitphase_qflaginfo.csv', TIRS_ID));
[error_ca, orbphase_qinfo] = ...
    read_PREFIRE_orbitphase_qflaginfo(pd.granule_beg_ts, info_file);
if (error_ca{2} ~= 0)
    fprintf(2, '%s\n', error_ca{1});
    return
end

if process_payload_tlm
   % Convert from detection units (count, DN) to uncalibrated spectral radiance
   %  (i.e., science units), with masking (of the masked channels common
   %  to all scenes) applied before calibration
   this_mask = TIRS.mask_blocked';
   for i_t=1:n_frame_types
      frametype = frame_types{i_t};
      sig.(frametype) = this_mask .* L0.(frametype).srd_DN;
   end

   % Determine instrument temperatures for all packets:
   %  Index of first dimension:
   %  1: 'Toroidal Mirror', 2: 'FPM Mounting Bracket', 3: 'Filter Assembly',
   %  4: 'Cal Target C', 5: 'Cal Target B', 6: 'Cal Target A',
   %  7: 'Moly Block B', 8: 'Moly Block A'
   for i_t=1:n_frame_types
      frametype = frame_types{i_t};
      L1.(frametype).T_val = L0.(frametype).TIRS_T;
   end

   odat.Cal_Artifacts.group_dims.atrack = L0.pdims.atrack_obstgt;
   odat.Cal_Artifacts.group_dims.calp_input_atr = length(L0.scipkt_frame_type);
   odat.Cal_Artifacts.group_dims.etemp = L0.pdims.n_TIRS_therm;
   odat.Cal_Artifacts.group_dims.mangoffs_atr = L0.pdims.mangoffs_atr;

   odat.Cal_Artifacts.group_atts.scipkt_ctime_start_s = L0.scipkt_ctime_start_s;
   odat.Cal_Artifacts.group_atts.scipkt_ctime_end_s = L0.scipkt_ctime_end_s;
   odat.Cal_Artifacts.group_atts.mirror_ang_offset0_deg = ...
                                                      L0.mirror_ang_offset0_deg;

   odat.Cal_Artifacts.scipkt_rel_ctime = L0.scipkt_rel_ctime;
   odat.Cal_Artifacts.scipkt_frame_type = L0.scipkt_frame_type;
   odat.Cal_Artifacts.engineering_temp = L1.obstgt.T_val;
   odat.Cal_Artifacts.mirror_aostart_rel_ctime = L0.mirror_aostart_rel_ctime;
   odat.Cal_Artifacts.mirror_ang_offset = L0.mirror_ang_offset;

   if L0.has_bus_bitflags
      odat.Cal_Artifacts.scipkt_bus_status_bitflags = ...
                                                  L0.scipkt_bus_status_bitflags;
      odat.Cal_Artifacts.scipkt_bus_events_bitflags = ...
                                                  L0.scipkt_bus_events_bitflags;
   end

   % Determine which calibration sequences are valid, and return relevant info
   %  from those:
   [error_ca, calseq] = process_TIRS_calseqs(L0, L1, sig);
   if (error_ca{2} ~= 0)
      fprintf(2, '%s\n', error_ca{1});
      return
   end
end

% Calculate approximate ctime of the midpoint of each TIRS image
%  integration ('L0.obstgt.ctime' is the *end* of each TIRS image integration):
L1.ctime = L0.obstgt.ctime-(0.5*TIRS.tau);  % [s]

odat.ProtoGeometry.orbit_phase_metric = (L1.ctime-pd.granule_beg_ts)/ ...
                            (pd.granule_end_ts-pd.granule_beg_ts)*360.;  % [deg]

if process_payload_tlm
   % Calculate calibration fields interpolated to the 'obstgt' frames' ctime:
   [error_ca, calseq_gain, gain_at_obs, space_at_obs, valid_cal] = ...
                    orbit_cal_interp(L0, calseq, L1, SRFMC, SRFMS, T_cal, TIRS);
   if (error_ca{2} ~= 0)
      fprintf(2, '%s\n', error_ca{1});
      return
   end

   % Calibrate image radiances, masking out disconnected/bad channels common to
   %  all scenes:
   tmp_rad = (sig.obstgt-space_at_obs)./gain_at_obs;  % [W/m^2/sr/um]
   tmp_rad = TIRS.mask_blocked'.*tmp_rad;             %

   % Set calibration "failures" and masked pixels to FillValue.  Any remaining
   %  negative radiances are clipped to be -5*NEdR (to allow for
   %  statistically-reasonable/relevant negative radiances, all > FillValue).
   for i=1:L0.pdims.xtrack
      for j=1:L0.pdims.allspectral
         if (TIRS.mask_blocked(i,j)==0) | (valid_cal(j,i)==0)
            tmp_rad(j,i,:) = -9.999e3;
         else
            threshold = -5.*NEDR.data(j,i);
            bad_ind = tmp_rad(j,i,:) < threshold;
            tmp_rad(j,i,bad_ind) = threshold;
         end
      end
   end

   % Replace any NaN values ("missing" data) with a fill value:
   tmp_rad = fillmissing(tmp_rad, 'constant', -9.999e3);

   odat.NonGeoLoc_Radiance.spectral_radiance = ...
                                            tmp_rad(2:end,:,:);  % [W/m^2/sr/um]
   odat.NonGeoLoc_Radiance.spectral_radiance_unc = ...
      repmat(NEDR.data(2:end,:), 1, 1, L0.pdims.atrack_obstgt);  % [W/m^2/sr/um]
   bad_ind = find(odat.NonGeoLoc_Radiance.spectral_radiance < -9.999e3*0.99);
   odat.NonGeoLoc_Radiance.spectral_radiance_unc(bad_ind) = -9.999e3;

   odat.NonGeoLoc_Radiance.detector_bitflags = dtbf(2:end,:);      % [-]

   odat.NonGeoLoc_Channel_0.channel_0_radiance = tmp_rad(1,:,:); % [W/m^2/sr]
   odat.NonGeoLoc_Channel_0.channel_0_radiance_unc = ...
      repmat(NEDR.data(1,:)', 1, L0.pdims.atrack_obstgt);         % [W/m^2/sr]
   bad_ind = find(odat.NonGeoLoc_Channel_0.channel_0_radiance < -9.999e3*0.99);
   odat.NonGeoLoc_Channel_0.channel_0_radiance_unc(bad_ind) = -9.999e3;
   odat.NonGeoLoc_Channel_0.channel_0_detector_bitflags = dtbf(1,:);

   calseq_mid_time = 0.5 * (calseq.cal_seqmean_ts + calseq.space_seqmean_ts);
   calibrator_T = mean(L1.obstgt.T_val(4:6,:), 1);
   odat.NonGeoLoc_Radiance.observation_bitflags = ...
           create_observation_bitflags(pd.granule_beg_ts, pd.granule_end_ts, ...
           L1.ctime, calseq_mid_time, odat.ProtoGeometry.orbit_phase_metric, ...
           orbphase_qinfo, calibrator_T, L0.scipkt_ctime, ti, TIRS_ID, ...
           L0.has_bus_bitflags, odat.Cal_Artifacts);

   % currently a calibration failure will flag the scene+channel
   % for the whole granule
   % Note that we skip channel 0 here, hence (j-1) index on bitflags.
   cal_bitflags = ...
       zeros(L0.pdims.cspectral, L0.pdims.xtrack, L0.pdims.atrack_obstgt, 'uint8');
   for i=1:L0.pdims.xtrack
      for j=2:L0.pdims.allspectral
         if valid_cal(j,i) == 0
            cal_bitflags(j-1,i,:) = 1;
         end
         if TIRS.mask_blocked(i,j) == 0
            cal_bitflags(j-1,i,:) = 2;
         end
      end
   end
   odat.NonGeoLoc_Radiance.calibration_bitflags = cal_bitflags;

   qflag = create_L1_quality_flags( ...
          odat.NonGeoLoc_Radiance.calibration_bitflags, ...
          dtbf(2:end,:), dtbf(1,:), ...
          odat.NonGeoLoc_Radiance.observation_bitflags );
   odat.NonGeoLoc_Radiance.radiance_quality_flag = qflag.rad;
   odat.NonGeoLoc_Radiance.calibration_quality_flag = qflag.calib;
   odat.NonGeoLoc_Radiance.detector_quality_flag = qflag.detector;
   odat.NonGeoLoc_Radiance.observation_quality_flag = qflag.obs;

   odat.NonGeoLoc_Channel_0.channel_0_radiance_quality_flag = qflag.channel_0_rad;
   odat.NonGeoLoc_Channel_0.channel_0_detector_quality_flag = qflag.channel_0_detector;

   tmp_gain = fillmissing(calseq_gain, 'constant', -9.999e3);
   tmp_gain(isinf(tmp_gain)) = -9.999e3;

   odat.Cal_Artifacts.framecount_caltgt = calseq.cal_framecount;
   odat.Cal_Artifacts.framecount_space = calseq.space_framecount;
   odat.Cal_Artifacts.DN_caltgt = permute(calseq.cal_seqmean_v(:,2:end,:),[2,3,1]);
   odat.Cal_Artifacts.DN_space = permute(calseq.space_seqmean_v(:,2:end,:),[2,3,1]);
   odat.Cal_Artifacts.gain = permute(tmp_gain(:,2:end,:), [2,3,1]);

   odat.Cal_Artifacts.channel_0_DN_caltgt = squeeze(calseq.cal_seqmean_v(:,1,:))';
   odat.Cal_Artifacts.channel_0_DN_space = squeeze(calseq.space_seqmean_v(:,1,:))';
   odat.Cal_Artifacts.channel_0_gain = squeeze(tmp_gain(:,1,:))';
   odat.Cal_Artifacts.caltgt_ctime = calseq.cal_seqmean_ts;
   odat.Cal_Artifacts.space_ctime = calseq.space_seqmean_ts;
   odat.Cal_Artifacts.calseq_engineering_temp = calseq.cal_seqmean_TIRS_T;

   spectral_BT = NaN(size(odat.NonGeoLoc_Radiance.spectral_radiance), ...
                     'like', odat.NonGeoLoc_Radiance.spectral_radiance);
   spectral_BT_unc = NaN(size(odat.NonGeoLoc_Radiance.spectral_radiance), ...
                         'like', odat.NonGeoLoc_Radiance.spectral_radiance);
   BT_quality_flag = odat.NonGeoLoc_Radiance.radiance_quality_flag;
   T_y = T_cal';
   for i=1:L0.pdims.xtrack
      for j=1:L0.pdims.allspectral-1  % 1-63
         % Only compute BT where spectral_radiance is not set to
         % its fill_value; this occurs for masked channels, or any
         % scene+channel combination with failed calibration calculation
         if (odat.NonGeoLoc_Radiance.spectral_radiance(j,i,1) > -9.999e3)
            rad_x0 = squeeze(SRFMC(j+1,i,:));
              % Need to make `makima` x values "unique", so where rad_x0 values
              %  are zero, set them to a "unique" value slightly more than zero:
            m = (rad_x0 < 1.e-30);
            rad_x0(m) = T_y(m)*1.e-30;
            Srad = makima(rad_x0, T_y);
            tmp_BT = ppval(Srad, ...
                             odat.NonGeoLoc_Radiance.spectral_radiance(j,i,:));
            tmp_BT(isinf(tmp_BT)) = 1.e36;  % Remove Inf
            tmp_BT(tmp_BT > 1.e36) = 1.e4;  % Clip very large values
            tmp_BT(tmp_BT < 1.e-30) = 1.e-30;  % Clip negative/zero values
            spectral_BT(j,i,:) = tmp_BT;

              % Interpolator for derivative (could use drad_dT, this seems
              %  simpler):
            D = [3 0 0 0;0 2 0 0;0 0 1 0]';
            dSrad_dL = Srad;
            dSrad_dL.order = 3;
            dSrad_dL.coefs = dSrad_dL.coefs*D;
            tmp_unc = ppval(dSrad_dL, ...
                         odat.NonGeoLoc_Radiance.spectral_radiance(j,i,:))* ...
                         NEDR.data(j+1,i);
            tmp_unc(isinf(tmp_unc)) = 1.e36;  % Remove Inf
            tmp_unc(tmp_unc > 1.e36) = 1.e36;  % Clip very large values
            tmp_unc(tmp_unc < -8.e3) = -8.e3;  % Clip very-negative values
            spectral_BT_unc(j,i,:) = tmp_unc;

            % replace BT derived from negative radiance with fill
            % value (these tend to have divergent values due to
            % cubic spline extrapolation)
            neg_value_msk = ...
                odat.NonGeoLoc_Radiance.spectral_radiance(j,i,:)<0;
            if any(neg_value_msk)
                spectral_BT(j,i,neg_value_msk) = -9.999e3;
                spectral_BT_unc(j,i,neg_value_msk) = -9.999e3;
                BT_quality_flag(j,i,neg_value_msk) = 2;
            end
         end
      end
   end

   % Replace any mising (NaN) values with a fill value:
   odat.NonGeoLoc_BT.spectral_BT = fillmissing(spectral_BT, 'constant', ...
                                               -9.999e3);  % [K]
   odat.NonGeoLoc_BT.spectral_BT_unc = fillmissing(spectral_BT_unc, ...
                                                   'constant', -9.999e3);  % [K]
   odat.NonGeoLoc_BT.BT_quality_flag = BT_quality_flag;

end

% Read in relevant L0-bus telemetry fields:
bt_frametype = 'bustlm';
[error_ca, L0, selected_atts, bt_fldnames] = read_PREFIRE_bus_L0(pd, L0, ...
                                                                 bt_frametype);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end
odat.global_atts.spacecraft_ID = selected_atts{1};
odat.global_atts.orbit_sim_version = selected_atts{2};

% Read in relevant L0-orbit (reconstruction) fields:
rt_frametype = 'orbrec';
[error_ca, L0, selected_atts, rt_fldnames] = read_PREFIRE_orbit_L0(pd, L0, ...
                                                                  rt_frametype);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

% Revise 'sat_solar_illumination_flag' field based on orbit reconstruction info:
cr_pp = makima(L0.orbrec.ctime(:), ...
               L0.orbrec.SC_in_Earth_umbra(:));
tmp_e = ppval(cr_pp, L0.bustlm.ctime(:));
ind = find(tmp_e < 0.001);
L0.bustlm.sat_solar_illumination_flag(ind) = 1;  % not in shadow
ind = find((tmp_e >= 0.001) & (tmp_e < 0.999));
L0.bustlm.sat_solar_illumination_flag(ind) = 0.5;  % partial shadow
ind = find((tmp_e >= 0.999));
L0.bustlm.sat_solar_illumination_flag(ind) = 0;  % in shadow

% Interpolate relevant bus and orbit fields to 'obstgt' ctime:
fldnames_to_interp = {'sat_solar_illumination_flag', 'beta_angle', ...
                      'q_scbody_wrt_ECI_1', 'q_scbody_wrt_ECI_2', ...
                      'q_scbody_wrt_ECI_3', 'q_scbody_wrt_ECI_4'};
ic = 0;
for i_v=1:length(fldnames_to_interp)
   ic = ic+1;
   c1_pp{ic} = makima(L0.(bt_frametype).ctime(:), ...
                      L0.(bt_frametype).(fldnames_to_interp{i_v})(:));
   odat.ProtoGeometry.(fldnames_to_interp{i_v}) = ppval(c1_pp{ic}, L1.ctime);
end
fldnames_to_interp = {'position_wrt_ECI_1', 'position_wrt_ECI_2', ...
                      'position_wrt_ECI_3', 'velocity_wrt_ECI_1', ...
                      'velocity_wrt_ECI_2', 'velocity_wrt_ECI_3'};
ic = 0;
for i_v=1:length(fldnames_to_interp)
   ic = ic+1;
   c2_pp{ic} = makima(L0.(rt_frametype).ctime(:), ...
                      L0.(rt_frametype).(fldnames_to_interp{i_v})(:));
   odat.ProtoGeometry.(fldnames_to_interp{i_v}) = ppval(c2_pp{ic}, L1.ctime);
end

% Determine pass type:
odat.ProtoGeometry.satellite_pass_type = ...
                                   sign(odat.ProtoGeometry.velocity_wrt_ECI_3);
odat.ProtoGeometry.satellite_pass_type( ...
                               odat.ProtoGeometry.satellite_pass_type == 0) = 1;

% Convert interpolated 'sat_solar_illumination_flag' field to a new integer
%  representation, where 0=no, 1=partial, 2=full:
ii = find(odat.ProtoGeometry.sat_solar_illumination_flag < 0.01);
odat.ProtoGeometry.sat_solar_illumination_flag(ii) = 0;
ii = find(odat.ProtoGeometry.sat_solar_illumination_flag > 0.99);
odat.ProtoGeometry.sat_solar_illumination_flag(ii) = 2;
ii = find(odat.ProtoGeometry.sat_solar_illumination_flag >= 0.01 & ...
          odat.ProtoGeometry.sat_solar_illumination_flag <= 0.99);
odat.ProtoGeometry.sat_solar_illumination_flag(ii) = 1;

% Calculate the mean value of the 'beta_angle' field and store it as a group
%  attribute:
odat.ProtoGeometry.group_atts.solar_beta_angle_deg = ...
                    mean(odat.ProtoGeometry.beta_angle, 'omitnan');  % [deg]
odat.ProtoGeometry = rmfield(odat.ProtoGeometry, 'beta_angle');  % do not write

% Calculate output ctime and UTC fields(leap-second-aware):
odat.ProtoGeometry.ctime = L1.ctime;  % [s since 2000-01-01T00:00:00 UTC]
[error_ca, UTC_DT, ctime_minus_UTC] = ctime_to_UTC_DT( ...
                                  odat.ProtoGeometry.ctime, 'seconds', tc, uc);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end
odat.ProtoGeometry.ctime_minus_UTC = ctime_minus_UTC;  % [s]
odat.ProtoGeometry.group_atts.TAI_minus_ctime_at_epoch_s = ...
                                       tc.ref_ctimeOffsetFromUTC_atEp_s;  % [s]
UTC_6splitvals = datevec(UTC_DT);
odat.ProtoGeometry.time_UTC_values = zeros([n_UTC_parts, ...
                                            length(ctime_minus_UTC)], 'int16');
for i=1:5
   odat.ProtoGeometry.time_UTC_values(i,:) = ...
                               UTC_6splitvals(:,i);  % Fill year through minute
end
odat.ProtoGeometry.time_UTC_values(6,:) = fix(UTC_6splitvals(:,6));  % [s]
s_frac = mod(UTC_6splitvals(:,6), 1);  % [s] frac part of seconds component
odat.ProtoGeometry.time_UTC_values(7,:) = s_frac*uc.s_to_ms;  % [ms]
ds_field = s_frac*uc.s_to_ds;  % [ds] used by obs_ID creation below

odat.ProtoGeometry.group_atts.start_granule_edge_ctime_s = pd.granule_beg_ts;
odat.ProtoGeometry.group_atts.end_granule_edge_ctime_s = pd.granule_end_ts;

odat.global_atts.ctime_coverage_start_s = ...
                                       L1.ctime(1);  % [s since 2000-01-01 UTC]
odat.global_atts.ctime_coverage_end_s = ...
                 L1.ctime(length(ctime_minus_UTC));  % [s since 2000-01-01 UTC]
odat.global_atts.UTC_coverage_start = char(UTC_DT(1), ...
                                           'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');
odat.global_atts.UTC_coverage_end = char(UTC_DT(length(ctime_minus_UTC)), ...
                                           'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');

% Construct 'obs_ID' integer representations:
a_sum = zeros(L0.pdims.atrack_obstgt, 1, 'int64');
reppwr = [13, 11, 9, 7, 5, 3];
reppwr_v = cast((10).^reppwr, 'int64');
for ip=1:length(reppwr)
    itmp = cast(odat.ProtoGeometry.time_UTC_values(ip,:), 'int64');
    a_sum = a_sum+itmp'.*reppwr_v(ip);
end
a_sum = a_sum+cast(fix(ds_field).*100, 'int64');
a_sum = a_sum+cast(TIRS_ID*10, 'int64');
odat.ProtoGeometry.obs_ID = zeros(L0.pdims.xtrack, L0.pdims.atrack_obstgt, ...
                                  'int64');
for i=1:L0.pdims.xtrack
   odat.ProtoGeometry.obs_ID(i,:) = a_sum+cast(i, 'int64');
end

chan = 1:63;
if process_payload_tlm
   % Channel reference info (ignoring first "broadband" channel):
   odat.NonGeoLoc_Radiance.detector_ID = zeros(L0.pdims.cspectral, ...
                                               L0.pdims.xtrack, 'int16');
   odat.NonGeoLoc_Radiance.wavelength = zeros(L0.pdims.cspectral, ...
                                              L0.pdims.xtrack, 'single');
   odat.NonGeoLoc_Radiance.idealized_wavelength = zeros(L0.pdims.cspectral, ...
                                                    L0.pdims.xtrack, 'single');
   for i=1:L0.pdims.xtrack
      odat.NonGeoLoc_Radiance.detector_ID(:,i) = cast(i*100, 'int64')+ ...
                                                 cast(chan, 'int64');
      odat.NonGeoLoc_Radiance.wavelength(:,i) = ...
         fillmissing(NEDR.channel_mean_wavelen(:,i), 'constant', -9.999e3);  % [um]
      odat.NonGeoLoc_Radiance.idealized_wavelength(:,i) = ...
         NEDR.channel_center_wavelen(:,i); % [um]
   end
else
   % Channel reference info (ignoring first "broadband" channel):
   odat.NonGeoLoc_FauxRad.detector_ID = zeros(L0.pdims.cspectral, ...
                                              L0.pdims.xtrack, 'int16');
   for i=1:L0.pdims.xtrack
      odat.NonGeoLoc_FauxRad.detector_ID(:,i) = cast(i*100, 'int64')+ ...
                                                cast(chan, 'int64');
   end
end

%% Write ("raw") L1A product file:

% Construct output filepath:
if process_payload_tlm
   fmtstr = 'raw-PREFIRE_SAT%1d_1A-RAD_%s_%04d%02d%02d%02d%02d%02d_%s.nc';
else
   fmtstr = 'raw-PREFIRE_SAT%1d_1A-GEOM_%s_%04d%02d%02d%02d%02d%02d_%s.nc';
end
output_fn = sprintf(fmtstr, TIRS_ID, pd.product_fullversion, ...
                    odat.ProtoGeometry.time_UTC_values(1:6,1), pd.granule_ID);
output_fpath = fullfile(pd.output_dir, output_fn);

% create aux data for testing.
% this outputs many of the MATLAB arrays directly, and also makes
% some diagnostic plots for one channel from all scenes
if (pd.make_L1A_aux_data == 1) & process_payload_tlm

    % copy out radiance and into into single variables to be able
    % to be written into the aux file (I don't want the whole odat
    % or L0 included.)
    cal_radiance = odat.NonGeoLoc_Radiance.spectral_radiance;
    obs_time = L1.ctime;

    % additional intermediate file:
    aux_fmtstr = 'raw-aux-PREFIRE_SAT%1d_1A-RAD_%s_%04d%02d%02d%02d%02d%02d_%s.mat';
    aux_output_fn = sprintf(aux_fmtstr, TIRS_ID, pd.product_fullversion, ...
                            odat.ProtoGeometry.time_UTC_values(1:6,1), pd.granule_ID);
    aux_output_fpath = fullfile(pd.output_dir, aux_output_fn);

    % write out intermediate file directly as Mat-format
    save(aux_output_fpath, ...
         'obs_time', 'sig', 'calseq', 'calseq_gain', ...
         'gain_at_obs', 'space_at_obs', 'cal_radiance');

    %quicklook plots: using a hardcoded channel ID (14) which
    %appears to be good for all 8 x 2 scenes.
    %TBD: make sure this works in a 'headless' batch process.
    channel_num = 14;
    aux_fmtstr = 'QL_L1A_PREFIRE_SAT%1d_1A-RAD_%s_%04d%02d%02d%02d%02d%02d_%s_ch%02dsc%d.png';
    for scene_num = 1:8
        plot_output_fn = sprintf(aux_fmtstr, TIRS_ID, pd.product_fullversion, ...
                                odat.ProtoGeometry.time_UTC_values(1:6,1), ...
                                pd.granule_ID, ...
                                channel_num, scene_num);
        plot_output_fpath = fullfile(pd.output_dir, plot_output_fn);
        fig = quicklook_L1A_data( ...
            sig.obstgt, obs_time, calseq, gain_at_obs, ...
            space_at_obs, cal_radiance, channel_num, scene_num);
        saveas(fig, plot_output_fpath);
        close(fig);
    end
end

% Load any additional data (including dimension lengths) into output structure
%  array:

[~, name, sfx] = fileparts(output_fpath);
odat.global_atts.file_name = [name sfx];

odat.global_atts.full_versionID = pd.product_fullversion;
odat.global_atts.granule_ID = pd.granule_ID;

line_parts = strsplit(pd.product_fullversion, '_');
if (line_parts{1}(1) == 'R')
   odat.global_atts.archival_versionID = strrep(line_parts{1}, 'R', '');
else  % Old nomenclature
   odat.global_atts.archival_versionID = strrep(line_parts{2}, 'R', '');
end

if process_payload_tlm
    odat.global_atts.SRF_NEdR_version = NEDR.im_version;  % for 1A-RAD
else
    odat.global_atts.SRF_NEdR_version = ' ';  % for 1A-GEOM
end

tmp_fpath = fullfile(pd.top_path, 'dist', ...
                     sprintf('prdgit_version_m%d.txt', pd.proc_mode));
fid = fopen(tmp_fpath, 'rt');
tmp_line = fgetl(fid);
fclose(fid);
idx = strfind(tmp_line, '(');
odat.global_atts.provenance = sprintf('%s%s (%s', ...
                   extractBefore(tmp_line, idx(1)), pd.product_fullversion, ...
                   extractAfter(tmp_line, idx(1)));

icf = 0;
for i_f=1:length(pd.ref_payload_L0_fpaths)
   icf = icf+1;
   [~, name, ext] = fileparts(pd.ref_payload_L0_fpaths{i_f});
   in_fnames{icf} = [name ext];
end
for i_f=1:length(pd.ref_bus_L0_fpaths)
   icf = icf+1;
   [~, name, ext] = fileparts(pd.ref_bus_L0_fpaths{i_f});
   in_fnames{icf} = [name ext];
end
odat.global_atts.input_product_files = strjoin(in_fnames, ', ');

tmp_fpath = fullfile(pd.top_path, 'VERSION.txt');
fid = fopen(tmp_fpath, 'rt');
odat.global_atts.processing_algorithmID = fgetl(fid);
fclose(fid);

odat.ProtoGeometry.group_atts.image_integration_duration_ms = ...
                                                   TIRS.tau*uc.s_to_ms;  % [ms]
odat.ProtoGeometry.group_dims.atrack = L0.pdims.atrack_obstgt;
odat.ProtoGeometry.group_dims.xtrack = L0.pdims.xtrack;
odat.ProtoGeometry.group_dims.UTC_parts = n_UTC_parts;

if process_payload_tlm
   odat.NonGeoLoc_Radiance.group_dims.atrack = L0.pdims.atrack_obstgt;
   odat.NonGeoLoc_Radiance.group_dims.xtrack = L0.pdims.xtrack;
   odat.NonGeoLoc_Radiance.group_dims.spectral = L0.pdims.cspectral;
   odat.NonGeoLoc_BT.group_dims.atrack = L0.pdims.atrack_obstgt;
   odat.NonGeoLoc_BT.group_dims.xtrack = L0.pdims.xtrack;
   odat.NonGeoLoc_BT.group_dims.spectral = L0.pdims.cspectral;
   odat.NonGeoLoc_Channel_0.group_dims.atrack = L0.pdims.atrack_obstgt;
   odat.NonGeoLoc_Channel_0.group_dims.xtrack = L0.pdims.xtrack;
   odat.Cal_Artifacts.group_dims.calseq = calseq.n;
   odat.Cal_Artifacts.group_dims.xtrack = L0.pdims.xtrack;
   odat.Cal_Artifacts.group_dims.spectral = L0.pdims.cspectral;
else
   odat.NonGeoLoc_FauxRad.group_dims.atrack = L0.pdims.atrack_obstgt;
   odat.NonGeoLoc_FauxRad.group_dims.xtrack = L0.pdims.xtrack;
   odat.NonGeoLoc_FauxRad.group_dims.spectral = L0.pdims.cspectral;
end

now_UTCn_DT = datetime('now', 'TimeZone', 'UTC', 'Format', ...
                      'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');  % faux-UTC (no leap-s)
odat.global_atts.UTC_of_file_creation = char(now_UTCn_DT);

% Prepare for NetCDF-format file write, then define and write the file:
[m_ncs, vars_to_write] = modify_nc_schema(L1A_schema, odat);
if isfile(output_fpath)
   delete(output_fpath);
end
repaired_ncwriteschema(output_fpath, m_ncs);
write_predef_nc_vars(output_fpath, vars_to_write);

end
