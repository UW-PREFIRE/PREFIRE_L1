function [error_ca] = produce_L1B_granulepart(pd)
% Produce part of an L1B granule, based on the input parameters
%  specified by a 'pd' struct returned from configure_toplevel_IO.
%
% Returns a cell array error code.

error_ca = {'#NONE#', 0};  % Default

uc = constants_unit_conversion;  % Unit conversion constants
tc = constants_time(pd.tools_anc_data_dir);  % Time-related constants
ti = TIRS_info;  % Some TIRS-related parameters

% Open input L1A granule data file, read in some useful parameters:
L1A.ncid = netcdf.open(pd.L1A_rad_fpath, 'NC_NOWRITE');
[L1A.n_dims, ~, L1A.n_globalatts, ~] = netcdf.inq(L1A.ncid);
gg_att_C = netcdf.getConstant('NC_GLOBAL');

  % If false, this is for 1B-*GEOM -- otherwise it is for 1B-RAD:
tmp_s = netcdf.getAtt(L1A.ncid, gg_att_C, 'input_product_files');
not_geom_only = ~contains(tmp_s, '_prototlm_');
use_DEM = ~strcmp(pd.DEM_root_dir, 'NONE');

% Read in product file data specifications (done here to have access to
%  fill_value parameters for the calculations below):
if not_geom_only
   JSON_fspecs_L1B_fpath = fullfile(pd.ancillary_data_dir, ...
                                     'L1B_product_filespecs.json');
else
   JSON_fspecs_L1B_fpath = fullfile(pd.ancillary_data_dir, ...
                                     'L1B_GEOM_product_filespecs.json');
end
[L1B_schema, L1B_jdat] = init_nc_schema_from_JSON(JSON_fspecs_L1B_fpath, false);

% Define/initialize output data structure array:
odat = struct;
odat.global_atts = struct;
odat.Geometry = struct;
odat.Geometry.group_atts = struct;
odat.Geometry.group_dims = struct;
if not_geom_only
   odat.Radiance = struct;
%   odat.Radiance.group_atts = struct;
   odat.Radiance.group_dims = struct;
   odat.BT = struct;
%   odat.BT.group_atts = struct;
   odat.BT.group_dims = struct;
   odat.Channel_0 = struct;
%   odat.Channel_0.group_atts = struct;
   odat.Channel_0.group_dims = struct;
else
   odat.FauxRad = struct;
   odat.FauxRad.group_dims = struct;
end

L1A.PG_gid = netcdf.inqNcid(L1A.ncid, 'ProtoGeometry');
if not_geom_only
   L1A.NR_gid = netcdf.inqNcid(L1A.ncid, 'NonGeoLoc_Radiance');
   L1A.NB_gid = netcdf.inqNcid(L1A.ncid, 'NonGeoLoc_BT');
   L1A.N0_gid = netcdf.inqNcid(L1A.ncid, 'NonGeoLoc_Channel_0');
   L1A.CA_gid = netcdf.inqNcid(L1A.ncid, 'Cal_Artifacts');
else
   L1A.NR_gid = netcdf.inqNcid(L1A.ncid, 'NonGeoLoc_FauxRad');
end
xid1 = L1A.ncid;
xid2 = L1A.ncid;
xid3 = L1A.ncid;
if (strcmp(netcdf.getAtt(L1A.ncid, gg_att_C, ...
           'is_this_a_final_product_file'), 'no'))
   xid1 = L1A.PG_gid;
   xid2 = L1A.NR_gid;
   xid3 = L1A.NB_gid;
end
dimid = netcdf.inqDimID(xid1, 'atrack');
[~, max_nframes] = netcdf.inqDim(xid1, dimid);
dimid = netcdf.inqDimID(xid1, 'xtrack');
[~, nxtrack] = netcdf.inqDim(xid1, dimid);
nUTCparts = 7;
dimid = netcdf.inqDimID(xid2, 'spectral');
[~, nspectral] = netcdf.inqDim(xid2, dimid);

% Copy some L1A global and group attribute values to the output structure array:
g_atts_to_copy = {'granule_ID', 'spacecraft_ID', 'sensor_ID', ...
                  'ctime_coverage_start_s', 'ctime_coverage_end_s', ...
                  'orbit_sim_version', 'SRF_NEdR_version'};
for ia=1:length(g_atts_to_copy)
   g_att_name = g_atts_to_copy{ia};
   g_att_value = netcdf.getAtt(L1A.ncid, gg_att_C, g_att_name);
   odat.global_atts.(g_att_name) = g_att_value;
end

gp_atts_to_copy = {'image_integration_duration_ms', 'solar_beta_angle_deg', ...
                   'start_granule_edge_ctime_s', 'end_granule_edge_ctime_s'};
for ia=1:length(gp_atts_to_copy)
   gp_att_name = gp_atts_to_copy{ia};
   odat.Geometry.group_atts.(gp_att_name) = ...
                              netcdf.getAtt(L1A.PG_gid, gg_att_C, gp_att_name);
end

% Read in instrument/sensor parameters and geometry:
sensor_num = str2num(odat.global_atts.sensor_ID(5:6));
TIRS = prefire_TIRS(ti.ROIC_tau, pd, 2, 1, sensor_num);
gm_fpath = fullfile(pd.instrument_model_dir, ...
                    ['geometry_instrument' num2str(sensor_num) '.mat']);
gm = load(gm_fpath);

% (pd.idxbeg_atrack:pd.idxend_atrack) is the inclusive range/subset to process.
% NOTE that when pd.idxend_atrack == -1, the actual atrack dimension length
%  (determined from the L1A file, since it varies across files) should be used
%  for the end index of the subset.
%
% This means that to process the first N frames, we expect idxbeg_atrack to be 1
% and idxend_atrack to be N.
frame_beg = pd.idxbeg_atrack;
frame_end = pd.idxend_atrack;
if (frame_end == -1)
    % if idxend is -1, that means do the whole file.
    % the inclusive end frame is then the max N frame in the L1A file.
    frame_end = max_nframes;
end
nframes = frame_end - frame_beg + 1;

% Read some ProtoGeometry data from L1A file:
[error_ca, L1A] = read_PREFIRE_L1A(L1A, L1A.PG_gid, {'obs_ID'}, ...
                                   [1,frame_beg], [nxtrack,nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

fields_to_read = {'position_wrt_ECI_1', 'position_wrt_ECI_2', ...
                  'position_wrt_ECI_3', 'velocity_wrt_ECI_1', ...
                  'velocity_wrt_ECI_2', 'velocity_wrt_ECI_3', ...
                  'q_scbody_wrt_ECI_1', 'q_scbody_wrt_ECI_2', ...
                  'q_scbody_wrt_ECI_3', 'q_scbody_wrt_ECI_4', ...
                  'ctime', 'sat_solar_illumination_flag', ...
                  'satellite_pass_type', 'orbit_phase_metric'};
[error_ca, L1A] = read_PREFIRE_L1A(L1A, L1A.PG_gid, fields_to_read, ...
                                   frame_beg, nframes);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

if not_geom_only
   fields_to_read = {'detector_ID', 'detector_bitflags', ...
                     'detector_quality_flag', ...
                     'wavelength', 'idealized_wavelength'};
   [error_ca, L1A] = read_PREFIRE_L1A(L1A, L1A.NR_gid, fields_to_read, ...
                                      [1, 1], [nspectral,nxtrack]);
else
   fields_to_read = {'detector_ID'};
   [error_ca, L1A] = read_PREFIRE_L1A(L1A, L1A.NR_gid, fields_to_read, ...
                                      [1, 1], [nspectral,nxtrack]);
end
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

if not_geom_only
   fields_to_read = {'spectral_radiance', 'spectral_radiance_unc', ...
                     'radiance_quality_flag', 'calibration_quality_flag', ...
                     'calibration_bitflags'};
   [error_ca, L1A] = read_PREFIRE_L1A(L1A, L1A.NR_gid, fields_to_read, ...
                               [1, 1, frame_beg], [nspectral,nxtrack,nframes]);
   if (error_ca{2} ~= 0)
      fprintf(2, '%s\n', error_ca{1});
      return
   end

   fields_to_read = {'observation_quality_flag', 'observation_bitflags'};
   [error_ca, L1A] = read_PREFIRE_L1A(L1A, L1A.NR_gid, fields_to_read, ...
                                      frame_beg, nframes);
   if (error_ca{2} ~= 0)
      fprintf(2, '%s\n', error_ca{1});
      return
   end

   fields_to_read = {'spectral_BT', 'spectral_BT_unc', 'BT_quality_flag'};
   [error_ca, L1A] = read_PREFIRE_L1A(L1A, L1A.NB_gid, fields_to_read, ...
                               [1, 1, frame_beg], [nspectral,nxtrack,nframes]);
   if (error_ca{2} ~= 0)
      fprintf(2, '%s\n', error_ca{1});
      return
   end

   fields_to_read = {'channel_0_radiance', 'channel_0_radiance_unc', ...
                     'channel_0_radiance_quality_flag', };
   [error_ca, L1A] = read_PREFIRE_L1A(L1A, L1A.N0_gid, fields_to_read, ...
                                            [1, frame_beg], [nxtrack,nframes]);
   if (error_ca{2} ~= 0)
      fprintf(2, '%s\n', error_ca{1});
      return
   end

   fields_to_read = {'channel_0_detector_bitflags', ...
                     'channel_0_detector_quality_flag'};
   [error_ca, L1A] = read_PREFIRE_L1A(L1A, L1A.N0_gid, fields_to_read, ...
                                      1, nxtrack);
   if (error_ca{2} ~= 0)
      fprintf(2, '%s\n', error_ca{1});
      return
   end

   varid = netcdf.inqVarID(L1A.CA_gid, 'mirror_aostart_rel_ctime');
   L1A.mirror_aostart_rel_ctime = netcdf.getVar(L1A.CA_gid, varid);
   varid = netcdf.inqVarID(L1A.CA_gid, 'mirror_ang_offset');
   L1A.mirror_ang_offset = netcdf.getVar(L1A.CA_gid, varid);
   L1A.mirror_ang_offset0_deg = netcdf.getAtt(L1A.CA_gid, ...
                     netcdf.getConstant('NC_GLOBAL'), 'mirror_ang_offset0_deg');
   L1A.scipkt_ctime_start_s = netcdf.getAtt(L1A.CA_gid, ...
                       netcdf.getConstant('NC_GLOBAL'), 'scipkt_ctime_start_s');
   L1A.scipkt_ctime_end_s = netcdf.getAtt(L1A.CA_gid, ...
                       netcdf.getConstant('NC_GLOBAL'), 'scipkt_ctime_end_s');
   L1A.mirror_aostart_ctime = double(L1A.mirror_aostart_rel_ctime)+ ...
                              double(L1A.scipkt_ctime_start_s);

else
   fields_to_read = {'radiance_obsq_bitflags'};
   [error_ca, L1A] = read_PREFIRE_L1A(L1A, L1A.NR_gid, fields_to_read, ...
                                      [1, frame_beg], [nxtrack,nframes]);
   if (error_ca{2} ~= 0)
      fprintf(2, '%s\n', error_ca{1});
      return
   end
end

% Now finished with input file reads, so close that file:
netcdf.close(L1A.ncid);

% Now reorganize into MATLAB structs similar to the old style.
sat = struct;
sat.R_eci = [L1A.position_wrt_ECI_1'; L1A.position_wrt_ECI_2'; ...
             L1A.position_wrt_ECI_3'];
sat.V_eci = [L1A.velocity_wrt_ECI_1'; L1A.velocity_wrt_ECI_2'; ...
             L1A.velocity_wrt_ECI_3'];
sat.quaternion = [L1A.q_scbody_wrt_ECI_1'; L1A.q_scbody_wrt_ECI_2'; ...
                  L1A.q_scbody_wrt_ECI_3'; L1A.q_scbody_wrt_ECI_4'];

  % Full counts of leap seconds (equivalent to TAI - UTC); needed for
  %  `process_sat_rv` below:
sat.ctime = L1A.ctime;  % [s]
[error_ca, UTC_DT, ctime_minus_UTC] = ctime_to_UTC_DT(sat.ctime, ...
                                                      'seconds', tc, uc);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end
sat.ctimeN_minus_UTC = double(ctime_minus_UTC)+ ...
                                       tc.ref_ctimeOffsetFromUTC_atEp_s;  % [s]
sat.TIRS_tau = TIRS.tau;  % [s]
sat.UTC_splitvals = datevec(UTC_DT);

time_UTC_values = zeros([nUTCparts, length(ctime_minus_UTC)], 'int16');
for i=1:5
   time_UTC_values(i,:) = sat.UTC_splitvals(:,i);  % Fill year through minute
end
time_UTC_values(6,:) = fix(sat.UTC_splitvals(:,6));  % [s]
s_frac = mod(sat.UTC_splitvals(:,6), 1);  % [s] frac part of seconds component
time_UTC_values(7,:) = s_frac*uc.s_to_ms;  % [ms]

UTC_granule_start = time_UTC_values(1:6,1);  % int16 array

%obtain basic parameters for Earth from the 1984 Geodesy data
earth_param_flag = 84;
[E_km, ~, ~, ~,] = def_Earth(earth_param_flag, 'kilometer');
[E_m, ~, ~, ~,] = def_Earth(earth_param_flag, 'meter');

% Create geolocation quality flag var:
geo_qf = zeros([nxtrack, nframes], 'uint16');

% takes satellite R and V vectors, computes ellipsoidal sat
% position (lat, lon, alt, to ellipsoid)
sat = process_sat_rv(sat, E_m, uc, tc);

% these are rotation angles, assuming my pointing vector sat.quaternion is an
%  'Euler-Rodrigues vector', yaw, pitch, roll
% based on TIRS geometry, the instrument will be positioned in this axis
% system yaw = local xy rotation, pitch = local xz rotation, roll = local yz
%  rotation, where x is the RAM direction
% thus the along track local dimension is x and the cross track is y
% it looks as though the BCT quaternion data (read into sat.quaternion) is
%  already attitude corrected compared to (unread) velocity field vectors, so
%  these values will be necessary only as offsets from sim velocities
ypr = zeros([3,nframes]);
ypr(3,:) = 0;  % roll adjustment [deg]

% The 180 deg yaw adjustment is calculated from
%   * the bus telemetry on SAT1 and SAT2 indicating that the power-positive
%      spacecraft attitude was with the ram direction along +X (spacecraft body
%      reference frame)
%   * pre-launch L1 algorithm development assumed that the ram direction would
%      be -X (spacecraft body reference frame)
% The -1.02 deg pitch adjustment is calculated from
%   * the nadir boresight vector provided by BCT
%   * in the spacecraft body reference frame, the ram direction is along +X
% Note that since during the slow mirror scan the mirror travels (w.r.t. the
%  spacecraft body reference frame) from about {a bit +X, +Y} to
%  {a lot -X, a bit +Y}, a positive encoder unit offset corresponds to a
%  negative pitch adjustment.
if (sensor_num == 1)  %TIRS1
   ypr(2,:) = -1.02;  % pitch adjustment [deg]
   ypr(1,:) = 180;  % yaw adjustment [deg]
else  % TIRS2
   ypr(2,:) = -1.02;  % pitch adjustment [deg]
   ypr(1,:) = 180;  % yaw adjustment [deg]
end

% Additional pitch adjustment(s):
[ypr] = get_mirror_ang_offset(ypr, L1A, nframes);

%satazel takes the satellite positions, the yaw/pitch/roll vector, and the
% instrument optical and geometric parameters to determine the azimuth and
% elevation angles or the 8 scenes and their corners, the centers are currently
% calculated redundantly

% TODO: redo ordering of the vertices. satazel is where the vertex
% line of sight angles are created. Changing the order (and going
% from 5 to 4) will impact process_ground_footprints and
% apply_terrain_correction.

% core processing steps:
% 1) compute relative line of sight angles the vertices and center
% of each scene.
sat = satazel(sat, ypr, TIRS, gm);

% 2) geolocate the vertices and centers to the Ellipsoidal surface
[ground, geo_qf] = process_ground_footprints(TIRS, sat, E_km, geo_qf, ypr);
nvertices = size(ground.Ps, 3)-1;

if use_DEM
   % 3) do a DEM lookup for the scene footprint surface elevations,
   % and correct the geolocated lat/lon positions
   % this is by far the most computationally expensive step.
   terrain_correction_verbose = false;

   [ground_corrected, geo_qf] = apply_terrain_correction(...
                 pd, ground, sat, terrain_correction_verbose, E_m, uc, geo_qf);
else  % mode intended to produce FOV estimates faster
   ground_corrected.lat = ground.P';
   ground_corrected.lon = ground.Q';
   tmp1 = squeeze(ground.Ps(1,:,:,:));  % (xtrack, n_vert, atrack)
   ground_corrected.lat_vertices = ...
                            permute(tmp1, [2,1,3]);  % (n_vert, xtrack, atrack)
   tmp1 = squeeze(ground.Qs(1,:,:,:));  % (xtrack, n_vert, atrack)
   ground_corrected.lon_vertices = ...
                            permute(tmp1, [2,1,3]);  % (n_vert, xtrack, atrack)

   tmp1 = squeeze(ground.Pi(1,:,:,:));  % (xtrack, n_vert, atrack)
   ground_corrected.lat_maxintgz_verts = ...
                            permute(tmp1, [2,1,3]);  % (n_vert, xtrack, atrack)
   tmp1 = squeeze(ground.Qi(1,:,:,:));  % (xtrack, n_vert, atrack)
   ground_corrected.lon_maxintgz_verts = ...
                            permute(tmp1, [2,1,3]);  % (n_vert, xtrack, atrack)
     % `compute_view_angles` below needs these fields to be set, so just set
     %  them to zero:
   ground_corrected.alt = zeros(size(ground_corrected.lat));
   ground_corrected.alt_std = zeros(size(ground_corrected.lat));
   % Note that 'ground_corrected' fields 'land_fraction', 'alt_min', 'alt_max'
   %  are left unset in this mode.
end

% 4) compute the view angles (sensor and solar azimuth/zenith)
ground_corrected = compute_view_angles(ground_corrected, sat, E_m, uc, tc);

% Load any additional data (including dimension lengths) into output structure
%  array:
odat.global_atts.full_versionID = pd.product_fullversion;

line_parts = strsplit(pd.product_fullversion, '_');
if (line_parts{1}(1) == 'R')
   odat.global_atts.archival_versionID = strrep(line_parts{1}, 'R', '');
else  % Old nomenclature
   odat.global_atts.archival_versionID = strrep(line_parts{2}, 'R', '');
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

odat.global_atts.UTC_coverage_start= char(UTC_DT(1), ...
                                            'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');
odat.global_atts.UTC_coverage_end = char(UTC_DT(length(ctime_minus_UTC)), ...
                                           'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');

[~, name, ext] = fileparts(pd.L1A_rad_fpath);
odat.global_atts.input_product_files = [name ext];

tmp_fpath = fullfile(pd.top_path, 'VERSION.txt');
fid = fopen(tmp_fpath, 'rt');
odat.global_atts.processing_algorithmID = fgetl(fid);
fclose(fid);

odat.Geometry.group_dims.atrack = nframes;
odat.Geometry.group_dims.xtrack = nxtrack;
odat.Geometry.group_dims.UTC_parts = nUTCparts;
odat.Geometry.group_dims.FOV_vertices = nvertices;

odat.Geometry.group_atts.TAI_minus_ctime_at_epoch_s = ...
                                       tc.ref_ctimeOffsetFromUTC_atEp_s;  % [s]

odat.Geometry.obs_ID = L1A.obs_ID;

odat.Geometry.time_UTC_values = time_UTC_values;
odat.Geometry.ctime = L1A.ctime;
odat.Geometry.ctime_minus_UTC = ctime_minus_UTC;

odat.Geometry.latitude = ground_corrected.lat;
odat.Geometry.longitude = ground_corrected.lon;

odat.Geometry.satellite_pass_type = L1A.satellite_pass_type;

  % Reorder/cull FOV vertices:
vll_fields = {'vertex_latitude', 'vertex_longitude', 'maxintgz_verts_lat', ...
              'maxintgz_verts_lon'};
for i=1:length(vll_fields)
   field = vll_fields{i};
   odat.Geometry.(field) = zeros([nvertices, nxtrack, nframes]);
end
iv_reorder = [4, 2, 3, 5];
for iv=1:nvertices
   odat.Geometry.vertex_latitude(iv,:,:) =  ...
                             ground_corrected.lat_vertices(iv_reorder(iv),:,:);
   odat.Geometry.vertex_longitude(iv,:,:) =  ...
                             ground_corrected.lon_vertices(iv_reorder(iv),:,:);
   odat.Geometry.maxintgz_verts_lat(iv,:,:) =  ...
                       ground_corrected.lat_maxintgz_verts(iv_reorder(iv),:,:);
   odat.Geometry.maxintgz_verts_lon(iv,:,:) =  ...
                       ground_corrected.lon_maxintgz_verts(iv_reorder(iv),:,:);
end

if use_DEM
   odat.Geometry.land_fraction = ground_corrected.land_fraction;
end
odat.Geometry.elevation = ground_corrected.alt;
odat.Geometry.elevation_stdev = ground_corrected.alt_std;

odat.Geometry.viewing_zenith_angle = ground_corrected.sensor_zenith;
odat.Geometry.viewing_azimuth_angle = ground_corrected.sensor_azimuth;

odat.Geometry.solar_zenith_angle = ground_corrected.solar_zenith;
odat.Geometry.solar_azimuth_angle = ground_corrected.solar_azimuth;
odat.Geometry.solar_distance = ground_corrected.solar_distance*uc.m_to_km; % [km]

odat.Geometry.subsat_latitude = sat.geod_lat;
odat.Geometry.subsat_longitude = sat.geod_lon;
odat.Geometry.sat_altitude = sat.geod_alt;

odat.Geometry.sat_solar_illumination_flag = L1A.sat_solar_illumination_flag;
odat.Geometry.orbit_phase_metric = L1A.orbit_phase_metric;

odat.Geometry.geoloc_quality_bitflags = geo_qf;

if not_geom_only
   odat.Radiance.group_dims.atrack = nframes;
   odat.Radiance.group_dims.xtrack = nxtrack;
   odat.Radiance.group_dims.spectral = nspectral;

   odat.Radiance.detector_ID = L1A.detector_ID;
   odat.Radiance.detector_bitflags = L1A.detector_bitflags;
   odat.Radiance.detector_quality_flag = L1A.detector_quality_flag;
   odat.Radiance.wavelength = L1A.wavelength;
   odat.Radiance.idealized_wavelength = L1A.idealized_wavelength;

   % flip all bad quality flag data to FillValue here.
   bad_ind = find(L1A.radiance_quality_flag == 2);
   odat.Radiance.spectral_radiance = L1A.spectral_radiance;
   odat.Radiance.spectral_radiance(bad_ind) = -9.999e3;
   odat.Radiance.spectral_radiance_unc = L1A.spectral_radiance_unc;
   odat.Radiance.spectral_radiance_unc(bad_ind) = -9.999e3;
   odat.Radiance.radiance_quality_flag = L1A.radiance_quality_flag;
     % Set top-level quality flag to FillValue where field is FillValue:
   odat.Radiance.radiance_quality_flag(bad_ind) = -99;
   odat.Radiance.observation_bitflags = L1A.observation_bitflags;
   odat.Radiance.calibration_bitflags = L1A.calibration_bitflags;
   odat.Radiance.observation_quality_flag = L1A.observation_quality_flag;
   odat.Radiance.calibration_quality_flag = L1A.calibration_quality_flag;

   odat.BT.group_dims.atrack = nframes;
   odat.BT.group_dims.xtrack = nxtrack;
   odat.BT.group_dims.spectral = nspectral;

   bad_ind = find(L1A.BT_quality_flag == 2);
   odat.BT.spectral_BT = L1A.spectral_BT;
   odat.BT.spectral_BT(bad_ind) = -9.999e3;
   odat.BT.spectral_BT_unc = L1A.spectral_BT_unc;
   odat.BT.spectral_BT_unc(bad_ind) = -9.999e3;
   odat.BT.BT_quality_flag = L1A.BT_quality_flag;
     % Set top-level quality flag to FillValue where field is FillValue:
   odat.BT.BT_quality_flag(bad_ind) = -99;

   odat.Channel_0.group_dims.atrack = nframes;
   odat.Channel_0.group_dims.xtrack = nxtrack;

   bad_ind = find(L1A.channel_0_radiance_quality_flag == 2);
   odat.Channel_0.channel_0_radiance = L1A.channel_0_radiance;
   odat.Channel_0.channel_0_radiance(bad_ind) = -9.999e3;
   odat.Channel_0.channel_0_radiance_unc = L1A.channel_0_radiance_unc;
   odat.Channel_0.channel_0_radiance_unc(bad_ind) = -9.999e3;
   odat.Channel_0.channel_0_detector_bitflags = ...
                                               L1A.channel_0_detector_bitflags;
   odat.Channel_0.channel_0_detector_quality_flag = ...
                                           L1A.channel_0_detector_quality_flag;
   odat.Channel_0.channel_0_radiance_quality_flag = ...
                                           L1A.channel_0_radiance_quality_flag;
     % Set top-level quality flag to FillValue where field is FillValue:
   odat.Channel_0.channel_0_radiance_quality_flag(bad_ind) = -99;

else
   odat.FauxRad.group_dims.atrack = nframes;
   odat.FauxRad.group_dims.xtrack = nxtrack;
   odat.FauxRad.group_dims.spectral = nspectral;

   odat.FauxRad.detector_ID = L1A.detector_ID;
   odat.FauxRad.radiance_obsq_bitflags = L1A.radiance_obsq_bitflags;
end

now_UTCn_DT = datetime('now', 'TimeZone', 'UTC', 'Format', ...
                      'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');  % faux-UTC (no leap-s)
odat.global_atts.UTC_of_file_creation = char(now_UTCn_DT);

% TODO:
% review array shapes in odat. It looks like things are not
% consistent about [atrack, xtrack] vs [xtrack, atrack]. I am not
% sure what it NEEDS to be here, so cannot fix is at the moment.

% Construct output filepath, inserting a suffix to describe the
% 'granule part' (as an inclusive atrack range.)
if not_geom_only
   fn_fmtstr = 'raw-PREFIRE_SAT%1d_1B-RAD_%s_%04d%02d%02d%02d%02d%02d_%s';
elseif use_DEM
   fn_fmtstr = 'raw-PREFIRE_SAT%1d_1B-GEOM_%s_%04d%02d%02d%02d%02d%02d_%s';
else
   fn_fmtstr = 'raw-PREFIRE_SAT%1d_1B-NTCG_%s_%04d%02d%02d%02d%02d%02d_%s';
end
output_fn = sprintf(fn_fmtstr, sensor_num, pd.product_fullversion, ...
                    UTC_granule_start(1:6,1), odat.global_atts.granule_ID);
suffix_fmtstr = '-%s_%05d_%05d_of_%05df.nc';

% Output frame range (in filename) should be the original inclusive
%  0-index values.
if pd.idxend_atrack == -1
   iend = frame_end-1;  % 0-based index
else
   iend = pd.idxend_atrack-1;  % 0-based index
end
part_suffix = sprintf(suffix_fmtstr, 'atrack', pd.idxbeg_atrack-1, iend, ...
                      max_nframes);
output_fpath = fullfile(pd.output_dir, [output_fn part_suffix]);

[~, name, ~] = fileparts(output_fpath);
odat.global_atts.file_name = name;

% Prepare for NetCDF-format file write, then define and write the file:
[m_ncs, vars_to_write] = modify_nc_schema(L1B_schema, odat);
if isfile(output_fpath)
   delete(output_fpath);
end

repaired_ncwriteschema(output_fpath, m_ncs);
write_predef_nc_vars(output_fpath, vars_to_write);

end
