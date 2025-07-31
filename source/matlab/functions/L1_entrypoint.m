function [successful_completion, error_ca] = L1_entrypoint(varargin)

successful_completion = false;

pd = configure_toplevel_IO;  % Get/set shared filepaths, dirs, control params

if pd.there_is_a_pre_error
   error_ca = {'there is a pre-error', 2};
   fprintf(2, '%s\n',  error_ca{1});
   return
end

switch pd.proc_mode
   case 3
      error_ca = produce_L1A_granule(pd);
   case 4
      error_ca = produce_L1B_granulepart(pd);
   otherwise
      tmp_str = sprintf('ERROR: unknown processing mode (%d)', pd.proc_mode);
      error_ca = {tmp_str, 2};
end

if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

% Output which MatLab toolboxes are used by this program
%license('inuse')

successful_completion = true;

end
