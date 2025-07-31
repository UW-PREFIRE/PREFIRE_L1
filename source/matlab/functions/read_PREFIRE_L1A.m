function [error_ca, L1A] = read_PREFIRE_L1A(L1A, gid, fields_to_read, ...
                                            start, count)
% Reads in one or more variables from an already-open L1A data file.
%
%-- Needed input parameters:
%   L1A :  (structure array)
%      .ncid :  ID of already-open L1A data file (NetCDF-format)
%   gid :  NetCDF-file group ID value
%   fields_to_read :  cell array of names of variables to read from the L1A file
%                     (NOTE: per call of this routine, all requested variables
%                      should have the same dimensions)
%   start :  scalar start index, or array of start indices (the number of values
%            provided should match the number of dimensions of the requested
%            variables)
%   count :  scalar element count, or array of element counts (the number of
%            values provided should match the number of dimensions of the
%            requested variables)
%
%-- Output:
%   Various new or updated fields of the structure array L1A, with those fields
%    named the same as given in 'fields_to_read'

error_ca = {'#NONE#', 0};  % Default

enotvar_C = netcdf.getConstant('NC_ENOTVAR');
for iv=1:length(fields_to_read)
   v_name = fields_to_read{iv};

   % Determine varid:
   try
      varid = netcdf.inqVarID(gid, v_name);
   catch
      msg = sprintf('ERROR: Variable %s does not exist in L1A file.', v_name);
      error_ca = {msg, 80};
      return
   end

   L1A.(v_name) = netcdf.getVar(gid, varid, start-1, count);  % C-indexing
end

end
