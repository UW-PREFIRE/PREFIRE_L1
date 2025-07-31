function [error_ca, L0, selected_atts, varnm_list] = ...
                                        read_PREFIRE_orbit_L0(pd, L0, frametype)
   
%% Extract/determine useful dimension values, open input files for reading:

% Creates/sets the L0.pdims.* fields :  values of named L0-orbit dimensions
fldname1 = strcat('atrack_', frametype);  % Along-track dim length
L0.pdims.(fldname1) = 0;  % Initialize sum 
for i_f=1:length(pd.orbit_L0_fpaths)
   ncfi.ncid{i_f} = netcdf.open(pd.orbit_L0_fpaths{i_f}, 'NC_NOWRITE');
   ncfi.gid{i_f} = netcdf.inqNcid(ncfi.ncid{i_f}, 'Orbit_L0_Data');
   finfo = ncinfo(pd.orbit_L0_fpaths{i_f});

   dim_names = {finfo.Dimensions.Name};
   dim_lengths = [finfo.Dimensions.Length];

   i_atrack = strcmp(dim_names, 'atrack');
   ncfi.(fldname1){i_f} = dim_lengths(i_atrack);
   tmp = L0.pdims.(fldname1)+dim_lengths(i_atrack);
   L0.pdims.(fldname1) = tmp;
end

%% Read in relevant L0 orbit reconstruction data fields:

% Fill structure array fields for the 'orbrec' frame type. For example:
%     L0.orbrec.timestamp :  timestamp of each entry, continuous_time
%                             [seconds since 2000-01-01T00:00:00 UTC]

[error_ca, L0, selected_atts_ca, varnm_list] = read_orbrec(frametype, L0, pd, ...
                                                           ncfi, true);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

if (all(strcmp(selected_atts_ca{1}, selected_atts_ca{1}{1})))
   spacecraft_ID = selected_atts_ca{1}{1};
   orbit_sim_v = selected_atts_ca{2}{1};
else
   error_ca = {['ERROR: all input L0-orbit files must be from the same ' ...
                'PREFIRE spacecraft.'], 14};
   return
end
selected_atts = {spacecraft_ID, orbit_sim_v};

% Close NetCDF input file(s):
for i_f=1:length(pd.orbit_L0_fpaths)
   netcdf.close(ncfi.ncid{i_f});
end

end


function [error_ca, L0, selected_atts_ca, varnm_list] = ...
                               read_orbrec(frametype, L0, pd, ncfi, subset_in_t)
% Extract data fields for the 'orbrec' frame type.

error_ca = {'#NONE#', 0};  % Default

% Extract spacecraft ID and orbit sim info (for further use in the calling
%  routine):
global_att_C = netcdf.getConstant('NC_GLOBAL');
for i_f=1:length(pd.orbit_L0_fpaths)
   ncid = ncfi.ncid{i_f};
   spacecraft_IDs{i_f} = netcdf.getAtt(ncid, global_att_C, 'spacecraft_ID');
   orbit_sim_v{i_f} = netcdf.getAtt(ncid, global_att_C, 'orbit_sim_version');
end
selected_atts_ca = {spacecraft_IDs, orbit_sim_v};

varnm_local0_l = {'SC_in_Earth_umbra'};
varnm_l = {'position_wrt_ECI_1', 'position_wrt_ECI_2', 'position_wrt_ECI_3', ...
           'velocity_wrt_ECI_1', 'velocity_wrt_ECI_2', 'velocity_wrt_ECI_3'};
varnm_list = [varnm_local0_l, varnm_l];

% Read (full) relevant ctime field; check whether there is any data:
[error, ts_full] = read_L0b_var('ctime', ncfi, pd, []);
if (error)
   error_ca = {'ERROR: No relevant orbit reconstruction field data is available.', 10};
   return
end

atr_dnm = strcat('atrack_', frametype);
if (subset_in_t)
   % Search ctime field for the (relative) packet index range that is
   %  needed to produce this granule:
   ip_beg = find(ts_full >= pd.granule_beg_ts, 1);
   ip_end = find(ts_full < pd.granule_end_ts, 1, 'last');
   if (isempty(ip_beg) | isempty(ip_end))
      error_ca = {['ERROR: The given orbit reconstruction data does not ' ...
                   'contain the requested granule buffer bounds.'], 11};
      return
   end
   if (ts_full(1) > pd.granule_beg_ts | ts_full(end) < pd.granule_end_ts)
      fprintf(1, 'WARNING: %s\n', ...
              'The input orbit reconstruction files given do not encompass the full time range requested.')
   end
   ip_range = [ip_beg, ip_end];

   % Add up to 5 entries on either end to help with higher-order interpolation
   %  later:
   ip_range(1) = max([ip_range(1)-5, 1]);
   ip_range(2) = min([ip_range(2)+5, length(ts_full)]);

   L0.pdims.(atr_dnm) = ip_range(2)-ip_range(1)+1;

   % Extract ctime field subset:
   L0.(frametype).ctime = ts_full(ip_range(1):ip_range(2));

   % Read time-subsets of any other relevant fields:
   for i_v=1:length(varnm_list)
      varnm = varnm_list{i_v};
      [error, subset_var] = read_L0b_var(varnm, ncfi, pd, ip_range);
      if (error)
         error_ca = {['ERROR: Some requested (subset) orbit reconstruction ' ...
                      'field data is unavailable.'], 12};
         return
      end
      L0.(frametype).(varnm) = subset_var;
   end

else

   L0.pdims.(atr_dnm) = length(ts_full);

   L0.(frametype).ctime = ts_full;

   % Read any other relevant fields:
   for i_v=1:length(varnm_list)
      varnm = varnm_list{i_v};
      [error, full_var] = read_L0b_var(varnm, ncfi, pd, []);
      if (error)
         error_ca = {['ERROR: Some requested orbit reconstruction field ' ...
                      'data is unavailable.'], 13};
         return
      end
      L0.(frametype).(varnm) = full_var;
   end
end

end


function [error, var_data] = read_L0b_var(varnm, ncfi, pd, ip_range)
% Extract data for a given field name, drawing from up to multiple
%  contiguous input files.  ctime-subsetting of the output field is enabled
%  if the function argument 'ip_range' is not an empty array.

% Determine input variable presence, shape, and length(s):
i_beg = 1;
for i_f=1:length(pd.orbit_L0_fpaths)
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
      odimlen = 0;
      [dname, dlen] = netcdf.inqDim(ncid, vdimids(1));
   else
      [dname, odimlen] = netcdf.inqDim(ncid, vdimids(1));
      [dname, dlen] = netcdf.inqDim(ncid, vdimids(2));
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

if (isempty(ip_range))  % Read in full file contributions to this field
   % Allocate array to hold all full file contributions to this field:
   if (n_dims == 1)
      var_data = zeros(n_ip, 1);  % a "column vector"
   else  % n_dims = 2
      var_data = zeros(odimlen, n_ip);
   end

   for i_f=1:length(pd.orbit_L0_fpaths)
      gid = ncfi.gid{i_f};
      ind = v_ind{i_f};

      v = netcdf.getVar(gid, varid);
      if (n_dims == 1)
         var_data(ind(1):ind(2)) = v(:);
      else  % n_dims = 2
         var_data(:,ind(1):ind(2)) = v(:,:);
      end
   end

else  % ctime-subsetting is desired -- allocate and fill subset of field
   % Allocate array to hold a subset of this field:
   n_ip = ip_range(2)-ip_range(1)+1;
   if (n_dims == 1)
      var_data = zeros(n_ip, 1);  % a "column vector"
   else  % n_dims = 2
      var_data = zeros(odimlen, n_ip);
   end

   vi = [v_ind{:}];
   vi_b = find(vi >= ip_range(1), 1);
   vi_e = find(vi <= ip_range(2), 1, 'last');
   i_f_i = [idivide(vi_b-1, int8(2))+1, idivide(vi_e-1, int8(2))+1];
   for i_f=1:length(v_ind)
      ind = v_ind{i_f};
      r_ind{i_f} = [ind(1)-ip_range(1), ind(2)-ip_range(2)];
   end

   dlnm = 'atrack_orbrec';
   ib_o = 1;
   for i_f=i_f_i(1):i_f_i(2)
      ind = r_ind{i_f};
      if ind(1) >= 0
         ib = 1;  % Start at the beginning of this file's data
      else
         ib = -ind(1)+1;  % Start partway through this file's data
      end
      if ind(2) > 0
         ie = ncfi.(dlnm){i_f}-ind(2);  % End partway through this file's data
      else
         ie = ncfi.(dlnm){i_f};  % End of this file's data
      end
      ie_o = ib_o+ie-ib;

      gid = ncfi.gid{i_f};

      v = netcdf.getVar(gid, varid);
      if (n_dims == 1)
         var_data(ib_o:ie_o) = v(ib:ie);
      else  % n_dims = 2
         var_data(:,ib_o:ie_o) = v(:,ib:ie);
      end
      ib_o = ie_o+1;
   end

end

end
