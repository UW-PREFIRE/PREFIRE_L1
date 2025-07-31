function [error_ca, L0, TIRS_ID] = read_PREFIRE_payload_L0(pd, ti, frame_types)

n_frame_types = length(frame_types);
   
%% Extract/determine useful dim values, open 'prepped' input files for reading:

% Creates/sets the L0.pdims.* fields :  values of named L0-payload dimensions

fldnames1 = {'atrack_c_scipkt'};  % Cumulative along-track dim lengths
fldnames2 = {'atrack_scipkt'};  % Subset along-track dim lengths
fldnames3 = {'atrack_f_scipkt'};  % Per-file along-track dim len
fldnames4 = {'n_CDH_therm', 'n_ROIC_therm', 'n_TIRS_therm', 'n_engd', ...
             'n_misc', 'allspectral', 'xtrack'};  % Other dimension values

L0.pdims.(fldnames1{1}) = 0;  % Initialize sum

for i_f=1:length(pd.payload_L0_fpaths)
   ncfi.ncid{i_f} = netcdf.open(pd.payload_L0_fpaths{i_f}, 'NC_NOWRITE');
   ncfi.gid{i_f} = netcdf.inqNcid(ncfi.ncid{i_f}, 'TIRS_L0_Data');
   finfo = ncinfo(pd.payload_L0_fpaths{i_f});

   % MATLAB-generated L0 payload (from inverse instrument model)
   % has dimension metadata within the group:
   if isempty(finfo.Dimensions)
      dim_names = {finfo.Groups(1).Dimensions.Name};
      dim_lengths = [finfo.Groups(1).Dimensions.Length];
   else
      dim_names = {finfo.Dimensions.Name};
      dim_lengths = [finfo.Dimensions.Length];
   end

   i_atrack = strcmp(dim_names, fldnames2{1});
   ncfi.(fldnames3{1}){i_f} = dim_lengths(i_atrack);
   tmp = L0.pdims.(fldnames1{1})+dim_lengths(i_atrack);
   L0.pdims.(fldnames1{1}) = tmp;
end
for i_d=1:length(fldnames4)
   i_n = strcmp(dim_names, fldnames4{i_d});
   L0.pdims.(fldnames4{i_d}) = dim_lengths(i_n);
end
L0.pdims.cspectral = L0.pdims.allspectral-1;  % number of channels in L1A output
L0.pdims.spectralxtrack = L0.pdims.allspectral*L0.pdims.xtrack;

%% Extract/determine useful dim values and frame-type fields:

ncid = netcdf.open(pd.pld_cat_fpath, 'NC_NOWRITE');
gid = netcdf.inqNcid(ncid, 'TIRS_L0_Pld_Categorization');
finfo = ncinfo(pd.pld_cat_fpath);
dim_names = {finfo.Dimensions.Name};
dim_lengths = [finfo.Dimensions.Length];

i_n = strcmp(dim_names, 'atrack_scipkt');
if dim_lengths(i_n) > L0.pdims.atrack_c_scipkt
   error_ca = {['ERROR: the atrack_scipkt dim of the frame-categorization ' ...
                'file is greater than the cumulative atrack_scipkt dim of ' ...
                'all input prepped L0-payload files.'], 45};
   TIRS_ID = -999;  % Allow error code propagation (avoid "not assigned")
   return
end

varid = netcdf.inqVarID(gid, 'ctime');
L0.scipkt_ctime = netcdf.getVar(gid, varid);
L0.scipkt_ctime_start_s = L0.scipkt_ctime(1);
L0.scipkt_ctime_end_s = L0.scipkt_ctime(length(L0.scipkt_ctime));
L0.scipkt_rel_ctime = L0.scipkt_ctime-L0.scipkt_ctime_start_s;

L0.mirror_ang_offset0_deg = netcdf.getAtt(gid, ...
                     netcdf.getConstant('NC_GLOBAL'), 'mirror_ang_offset0_deg');

varid = netcdf.inqVarID(gid, 'mirror_aostart_ctime');
tmp_ct = netcdf.getVar(gid, varid);
L0.mirror_aostart_rel_ctime = tmp_ct-L0.scipkt_ctime_start_s;

varid = netcdf.inqVarID(gid, 'mirror_ang_offset');
L0.mirror_ang_offset = netcdf.getVar(gid, varid);
tmp_d = size(L0.mirror_ang_offset);
L0.pdims.mangoffs_atr = tmp_d(1);

try
   varid = netcdf.inqVarID(gid, 'scipkt_frame_type');
catch
   varid = netcdf.inqVarID(gid, 'scipkt_look_type');  % In pre-Oct2024 files
end
L0.scipkt_frame_type = netcdf.getVar(gid, varid);

try  % Pre-Oct2024 files do not have these fields
   varid = netcdf.inqVarID(gid, 'scipkt_bus_status_bitflags');
   L0.scipkt_bus_status_bitflags = netcdf.getVar(gid, varid);
   varid = netcdf.inqVarID(gid, 'scipkt_bus_events_bitflags');
   L0.scipkt_bus_events_bitflags = netcdf.getVar(gid, varid);
   L0.has_bus_bitflags = true;
catch
   L0.has_bus_bitflags = false;
end

netcdf.close(ncid);

ind.obstgt = find(L0.scipkt_frame_type == ti.OBSTGT_NOM);
ind.space = find(L0.scipkt_frame_type == ti.SPACE_NOM);
ind.caltgt = find(L0.scipkt_frame_type == ti.CALTGT_NOM);

%% Read in relevant L0 payload data fields, reorder spectroradiometer DNs:

% Fill structure array fields for each of the three main frame types:
%  'caltgt', 'space', and 'obstgt'.  For example, for the 'caltgt' frame type:
%     L0.caltgt.ctime :  timestamp of each frame's end, continuous_time
%                         [seconds since 2000-01-01T00:00:00 UTC]
%     L0.caltgt.TIRS_T :  TIRS thermistor temperatures
%     L0.caltgt.adj_srd_DN :  raw spectroradiometer counts (DNs)
%     L0.caltgt.srd_DN :  same as adj_srd_DN but with permuted dimension order

for i_t=1:n_frame_types
   frametype = frame_types{i_t};

     % t_subset = [ctime_framecat_start, ctime_subset_start, ctime_subset_end,
     %             ctime_frame_offset_for_subset_search]
   if strcmp(frametype, 'obstgt')
      t_subset = [L0.scipkt_ctime_start_s, pd.granule_beg_ts, ...
                  pd.granule_end_ts, -0.5*ti.ROIC_tau];  % [s]
   else
      t_subset = [L0.scipkt_ctime_start_s, L0.scipkt_ctime_start_s, ...
                  L0.scipkt_ctime_end_s, 0];  % [s]
   end
   [error_ca, L0, sensor_IDs] = read_frame_type(frametype, ind.(frametype), ...
                                                L0, pd, ncfi, t_subset);
   if (error_ca{2} ~= 0)
      fprintf(2, '%s\n', error_ca{1});
      TIRS_ID = -999;  % Allow error code propagation (avoid "not assigned")
      return
   end

   if (i_t == 1)
      if (all(strcmp(sensor_IDs, sensor_IDs{1})))
         tmp = sensor_IDs{1};
         TIRS_ID = str2num(tmp(5:6));
      else
         error_ca = {['ERROR: all input L0-payload files must be from the ' ...
                      'same TIRS sensor.'], 44};
         TIRS_ID = -999;  % Allow error code propagation (avoid "not assigned")
         return
      end
   end

   % Reorder 3-D array
   L0.(frametype).srd_DN = permute(L0.(frametype).adj_srd_DN, [3 2 1]);
end

% Close NetCDF input file(s):
for i_f=1:length(pd.payload_L0_fpaths)
   netcdf.close(ncfi.ncid{i_f});
end

if length(L0.obstgt.ctime) < 2
   error_ca = {['ERROR: less than 2 obstgt frames within the given ' ...
                'time range.'], 46};
   TIRS_ID = -999;  % Allow error code propagation (avoid "not assigned")
   return
elseif length(L0.caltgt.ctime) < 6
   error_ca = {['ERROR: less than 6 caltgt frames within the given ' ...
                'time range (i.e., calibration not possible).'], 47};
   TIRS_ID = -999;  % Allow error code propagation (avoid "not assigned")
   return
elseif length(L0.space.ctime) < 6
   error_ca = {['ERROR: less than 6 space frames within the given ' ...
                'time range (i.e., calibration not possible).'], 48};
   TIRS_ID = -999;  % Allow error code propagation (avoid "not assigned")
   return
end

end


function [error_ca, L0, sensor_IDs] = read_frame_type(frametype, ...
                                          frametype_ind, L0, pd, ncfi, t_subset)
% Extract data fields for a given TIRS frame type.

error_ca = {'#NONE#', 0};  % Default

% Extract sensor ID info (for further use in the calling routine):
for i_f=1:length(pd.payload_L0_fpaths)
   ncid = ncfi.ncid{i_f};
   sensor_IDs{i_f} = netcdf.getAtt(ncid, ...
                                 netcdf.getConstant('NC_GLOBAL'), 'sensor_ID');
end

varnm_list = {'TIRS_T', 'adj_srd_DN', 'ctime'};

% Read (full) ctime field, without subsetting by frame type, and check whether
%  there is any data:
max_num_frames = sum([ncfi.atrack_f_scipkt{:}]);
dummy_ind = (1:max_num_frames)';
[error, ts_full] = read_L0p_var('ctime', dummy_ind, ncfi, pd, {[]});
if (error)
   error_ca = {['ERROR: No relevant payload telemetry field data is ' ...
                'available.'], 40};
   return
end

atr_dnm = strcat('atrack_', frametype);
% Search ctime field for the (relative) packet index range that is
%  needed to produce this granule:
ip_beg = find(ts_full+t_subset(4) > t_subset(2)-0.001, 1);
ip_end = find(ts_full+t_subset(4) < t_subset(3)+0.001, 1, 'last');
if (isempty(ip_beg) | isempty(ip_end))
   error_ca = {['ERROR: The given payload telemetry data does not ' ...
                   'contain the requested granule buffer bounds.'], 41};
   return
end
if (ts_full(1) > t_subset(2) | ts_full(end) < t_subset(3))
   fprintf(1, 'WARNING: %s\n', ...
              'The input payload telemetry files given do not encompass the full time range requested.');
end
ip_range = [ip_beg, ip_end];  % <ts_full index>

% * frametype_ind = 1 is a <relative index> valid at t_subset(1)
% * t_subset(1) should be >= ts_full(1)]
%
% * ip_range(1) is the <ts_full index> with a ctime > t_subset(2)
% * ip_range(2) is the <ts_full index> with a ctime < t_subset(3)
%
% * ip_refbeg is the <ts_full index> with a ctime >= t_subset(1)
% * mod_frametype_ind is the <ts_full index> for each frametype_ind entry

ip_refbeg = find(ts_full >= t_subset(1), 1);  % <ts_full index>
mod_frametype_ind = frametype_ind+ip_refbeg-1;  % <ts_full index>

frametype_rabsind = mod_frametype_ind(mod_frametype_ind >= ip_range(1) & ...
                                      mod_frametype_ind <= ip_range(2));
L0.pdims.(atr_dnm) = length(frametype_rabsind);

% Read time-subsets of any other relevant fields (including ctime)
for i_v=1:length(varnm_list)
   varnm = varnm_list{i_v};
   [error, subset_var] = read_L0p_var(varnm, mod_frametype_ind, ncfi, pd, ...
                                      {ip_range, frametype_rabsind});
   if (error)
      error_ca = {['ERROR: Some requested (subset) payload telemetry ' ...
                   'field data is unavailable.'], 42};
      return
   end
   L0.(frametype).(varnm) = subset_var;
end

end


function [error, var_data] = read_L0p_var(varnm, frametype_absind, ncfi, pd, ...
                                          subset_in_t_cfg)
% Extract data for a given field name and frame type indices, drawing from up to
%  multiple contiguous input files.  ctime-subsetting of the output field is
%  enabled if the function argument 'subset_in_t_cfg' does not contain an empty
%  array.

% Determine input variable presence, shape, and length(s):
i_beg = 1;
for i_f=1:length(pd.payload_L0_fpaths)
   ncid = ncfi.ncid{i_f};
   gid = ncfi.gid{i_f};

   % Determine varid; if it is not found in this file, move on to the next
   %  file (if any):
   try
      varid = netcdf.inqVarID(gid, varnm);
   catch
      v_ind{i_f} = [0, 0];
      continue
   end

   % Determine some info about this variable:
   [varn, vtype, vdimids, vnumatts] = netcdf.inqVar(gid, varid);
   n_dims = length(vdimids);
   if (n_dims == 1)
      odimlen1 = 0;
      [dname, dlen] = netcdf.inqDim(ncid, vdimids(1));
   elseif (n_dims == 2)
      [dname, odimlen1] = netcdf.inqDim(ncid, vdimids(1));
      [dname, dlen] = netcdf.inqDim(ncid, vdimids(2));
   else  % n_dims == 3, only for 'adj_srd_DN'
      [dname, dlen] = netcdf.inqDim(ncid, vdimids(1));
      [dname, odimlen1] = netcdf.inqDim(ncid, vdimids(2));
      [dname, odimlen2] = netcdf.inqDim(ncid, vdimids(3));
   end

   % Multifile concatenation indices:
   i_end = i_beg+dlen-1;
   v_ind{i_f} = [i_beg, i_end];

   i_beg = i_end+1;
end
n_ip = i_end;

error = (i_beg == 1);
if (error)
   var_data = -999;  % Allow error code propagation (avoid "not assigned")
   return  % No field data is available
end

ip_range = subset_in_t_cfg{1};

if (isempty(ip_range))  % Read in full file contributions to this field
   % Allocate array to hold all full file contributions to this field:
   if (n_dims == 1)
      var_data0 = zeros(n_ip, 1);  % a "column vector"
   elseif (n_dims == 2)
      var_data0 = zeros(odimlen1, n_ip);
   else  % n_dims == 3, only for 'adj_srd_DN'
      var_data0 = zeros(n_ip, odimlen1, odimlen2);
   end

   for i_f=1:length(pd.payload_L0_fpaths)
      gid = ncfi.gid{i_f};
      ind = v_ind{i_f};

      v = netcdf.getVar(gid, varid);
      if (n_dims == 1)
         var_data0(ind(1):ind(2)) = v(:);
      elseif (n_dims == 2)
         var_data0(:,ind(1):ind(2)) = v(:,:);
      else  % n_dims == 3, only for 'adj_srd_DN'
         var_data0(ind(1):ind(2),:,:) = v(:,:,:);
      end
   end

   % Retrieve the frame type subset:
   if (n_dims == 1)
      var_data = var_data0(frametype_absind);
   elseif (n_dims == 2)
      var_data = var_data0(:,frametype_absind);
   else  % n_dims == 3, only for 'adj_srd_DN'
      var_data = var_data0(frametype_absind,:,:);
   end

else  % ctime-subsetting is desired -- allocate and fill subset of field
   % Allocate array to hold a subset of this field:
   n_ip = ip_range(2)-ip_range(1)+1;
   if (n_dims == 1)
      var_data0 = zeros(n_ip, 1);  % a "column vector"
   elseif (n_dims == 2)
      var_data0 = zeros(odimlen1, n_ip);
   else  % n_dims == 3, only for 'adj_srd_DN'
      var_data0 = zeros(n_ip, odimlen1, odimlen2);
   end

   vi = [v_ind{:}];
   vi_b = find(vi >= ip_range(1), 1);
   vi_e = find(vi <= ip_range(2), 1, 'last');
   i_f_i = [idivide(vi_b-1, int8(2))+1, idivide(vi_e-1, int8(2))+1];
   for i_f=1:length(v_ind)
      ind = v_ind{i_f};
      r_ind{i_f} = [ind(1)-ip_range(1), ind(2)-ip_range(2)];
   end

   ib_o = 1;
   for i_f=i_f_i(1):i_f_i(2)
      ind = r_ind{i_f};
      if ind(1) >= 0
         ib = 1;  % Start at the beginning of this file's data
      else
         ib = -ind(1)+1;  % Start partway through this file's data
      end
      if ind(2) > 0
         ie = ncfi.atrack_f_scipkt{i_f}-ind(2);  % End partway through this file
      else
         ie = ncfi.atrack_f_scipkt{i_f};  % End of this file's data
      end
      ie_o = ib_o+ie-ib;

      gid = ncfi.gid{i_f};

      v = netcdf.getVar(gid, varid);
      if (n_dims == 1)
         var_data0(ib_o:ie_o) = v(ib:ie);
      elseif (n_dims == 2)
         var_data0(:,ib_o:ie_o) = v(:,ib:ie);
      else  % n_dims == 3, only for 'adj_srd_DN'
         var_data0(ib_o:ie_o,:,:) = v(ib:ie,:,:);
      end
      ib_o = ie_o+1;
   end

   % Retrieve the frame type subset:
   frametype_relind = subset_in_t_cfg{2}-ip_range(1)+1;
   if (n_dims == 1)
      var_data = var_data0(frametype_relind);
   elseif (n_dims == 2)
      var_data = var_data0(:,frametype_relind);
   else  % n_dims == 3, only for 'adj_srd_DN'
      var_data = var_data0(frametype_relind,:,:);
   end

end

end
