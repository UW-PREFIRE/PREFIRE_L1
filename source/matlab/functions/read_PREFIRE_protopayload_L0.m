function [error_ca, L0, TIRS_ID] = read_PREFIRE_protopayload_L0(pd, frame_types)

n_frame_types = length(frame_types);
   
%% Extract/determine useful dimension values, open input files for reading:

% Creates/sets the L0.pdims.* fields :  values of named L0-payload dimensions

fldnames1 = {'atrack_c_scipkt'};  % Cumulative along-track dim lengths
fldnames2 = {'atrack_scipkt'};  % Subset along-track dim lengths
fldnames3 = {'atrack_f_scipkt'};  % Per-file along-track dim len
fldnames4 = {'allspectral', 'xtrack'};  % Other dimension values

L0.pdims.(fldnames1{1}) = 0;  % Initialize sum

for i_f=1:length(pd.payload_L0_fpaths)
   ncfi.ncid{i_f} = netcdf.open(pd.payload_L0_fpaths{i_f}, 'NC_NOWRITE');
   ncfi.gid{i_f} = netcdf.inqNcid(ncfi.ncid{i_f}, 'TIRS_L0_Data');
   finfo = ncinfo(pd.payload_L0_fpaths{i_f});

   dim_names = {finfo.Dimensions.Name};
   dim_lengths = [finfo.Dimensions.Length];

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

%% Read in relevant L0 protopayload data fields:

% Fill structure array fields for each of the three main frame types: 'caltgt',
%  'space', and 'obstgt'.  For example, for the 'caltgt' frame type:
%     L0.caltgt.ctime :  timestamp of each frame's end, continuous_time
%                         [seconds since 2000-01-01T00:00:00 UTC]

for i_t=1:n_frame_types
   frametype = frame_types{i_t};

   subset_in_t = strcmp(frametype, 'obstgt');
   [error_ca, L0, sensor_IDs] = read_protoframe_type(frametype, L0, pd, ncfi, ...
                                                    subset_in_t);
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
end

end


function [error_ca, L0, sensor_IDs] = read_protoframe_type(frametype, L0, pd, ...
                                                          ncfi, subset_in_t)
% Extract data fields for a given TIRS frame type.

error_ca = {'#NONE#', 0};  % Default

% Extract sensor ID info (for further use in the calling routine):
for i_f=1:length(pd.payload_L0_fpaths)
   ncid = ncfi.ncid{i_f};
   sensor_IDs{i_f} = netcdf.getAtt(ncid, ...
                                 netcdf.getConstant('NC_GLOBAL'), 'sensor_ID');
end

% Read (full) ctime field, without subsetting by frametype
% check whether there is any data:
[error, ts_full] = read_L0proto_var('ctime', frametype, ncfi, pd, {[]});
if (error)
   error_ca = {['ERROR: No relevant protopayload telemetry field data is ' ...
                'available.'], 40};
   return
end

atr_dnm = strcat('atrack_', frametype);
if (subset_in_t)
   % Search ctime field for the (relative) packet index range that is
   %  included in this granule:
   ip_beg = find(ts_full >= pd.granule_beg_ts, 1);
   ip_end = find(ts_full < pd.granule_end_ts, 1, 'last');
   if (isempty(ip_beg) | isempty(ip_end))
      error_ca = {['ERROR: The given protopayload telemetry data does not (or ' ...
                   'not fully) contain the requested granule bounds.'], 41};
      return
   end
   if (ts_full(1) > pd.granule_beg_ts | ts_full(end) < pd.granule_end_ts)
      fprintf(1, 'WARNING: %s\n', ...
              'The input protopayload telemetry files given do not encompass the full time range requested.');
   end
   ip_range = [ip_beg, ip_end];

   L0.pdims.(atr_dnm) = ip_range(2)-ip_range(1)+1;

   % Extract ctime field subset:
   L0.(frametype).ctime = ts_full(ip_range(1):ip_range(2));

else

   L0.pdims.(atr_dnm) = length(ts_full);

   L0.(frametype).ctime = ts_full;

end

end


function [error, var_data] = read_L0proto_var(varnm, frametype, ncfi, pd, ...
                                              subset_in_t_cfg)
% Extract data for a given field name and frame type, drawing from up to multiple
%  contiguous input files.  ctime-subsetting of the output field is enabled
%  if the function argument 'subset_in_t_cfg' does not contain an empty
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
   else
      [dname, odimlen1] = netcdf.inqDim(ncid, vdimids(1));
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

ip_range = subset_in_t_cfg{1};

if (isempty(ip_range))  % Read in full file contributions to this field
   % Allocate array to hold all full file contributions to this field:
   if (n_dims == 1)
      var_data = zeros(n_ip, 1);  % a "column vector"
   else  % n_dims = 2
      var_data = zeros(odimlen1, n_ip);
   end

   for i_f=1:length(pd.payload_L0_fpaths)
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
      var_data = zeros(odimlen1, n_ip);
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
         var_data(ib_o:ie_o) = v(ib:ie);
      else  % n_dims = 2
         var_data(:,ib_o:ie_o) = v(:,ib:ie);
      end
      ib_o = ie_o+1;
   end

end

end
