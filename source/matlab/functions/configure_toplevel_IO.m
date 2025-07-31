function [pd] = configure_toplevel_IO
% Assigns filepaths, directories/paths, and other control parameters to the
%  fields of a new structure array 'pd'
pd.there_is_a_pre_error = false;

if isunix

   [func_path, ~, ~] = fileparts(which('configure_toplevel_IO'));
   top_path = fullfile(func_path, '..', '..', '..'); % has 'dist','source', etc.
   begDir = cd(top_path);
   pd.top_path = pwd;  % Resolve '.' and '..'
   cd(begDir);
   top_path = pd.top_path;

   tmp = getenv('PROC_MODE');
   if strlength(tmp) < 1  % relevant environment variables are not set
      pd.proc_mode = 3;

      pd.ancillary_data_dir = fullfile(top_path, 'dist/ancillary');
      pd.tools_anc_data_dir = fullfile(top_path, 'PREFIRE_tools/dist/ancillary');
      pd.instrument_model_dir = fullfile(top_path, ...
                                         'dist/ancillary/instrument_model');
      pd.output_dir = fullfile(top_path, 'test/outputs');

      pd.product_fullversion = 'R01_P00';

      if (pd.proc_mode == 3)
         pd.SRF_disambig_str = 'SRF_v13_2024-09-15';

         [pd.payload_L0_fpaths{1:2}] = deal( ...
              fullfile(top_path, ...
              'test/inputs/prefire_02_payload_tlm_20221028193650_20221029032453_20240324215847.nc'), ...
              fullfile(top_path, ...
              'test/inputs/prefire_02_payload_tlm_20221029032450_20221029135736_20240324215903.nc'));
         pd.ref_payload_L0_fpaths = pd.payload_L0_fpaths;

         [pd.bus_L0_fpaths{1}] = deal( ...
              fullfile(top_path, ...
                      'test/inputs/prefire_02_bus_tlm_20221029145913_20221029164727_20230330143233.nc'));
         pd.ref_bus_L0_fpaths = pd.bus_L0_fpaths;

         [pd.orbit_L0_fpaths{1}] = deal( ...
              fullfile(top_path, ...
                      'test/inputs/prefire_02_orbit_reconst_20221029145913_20221029164727_20230330143233.nc'));
         pd.ref_orbit_L0_fpaths = pd.orbit_L0_fpaths;

         pd.pld_cat_fpath = fullfile(top_path, ...
                 'test/inputs/prefire_02_pld_cat_20221028193650_20221029135736_20240324214512.nc');

         [pd.granule_beg_ts, pd.granule_end_ts] = deal({717502532095.872, ...
                                                        717508255640.243});
         pd.granule_ID = '00009';

      elseif (pd.proc_mode == 4)

         pd.L1A_rad_fpath = 'test/outputs/m3/raw-PREFIRE_SAT1_1A-RAD_R01_P00_20220926150106_00009.nc';
         pd.DEM_root_dir = '/data/users/mmm/DEM/copernicus-dem-90m/tiles';

         atrack_idx_range_0bi = 'ATRACK_IDXRANGE_0BASED_INCLUSIVE:0:END';
         tmp = strsplit(atrack_idx_range_0bi, ':');
         pd.idxbeg_atrack = str2num(tmp{2})+1;  % 1-based index
         if strcmp(tmp{3}, 'END')
            pd.idxend_atrack = -1;  % Sentinel value
         else
            pd.idxend_atrack = str2num(tmp{3})+1;  % 1-based index
         end

      end

   else
      pd = obtain_from_env_vars(pd, tmp);
   end

else  % MS-Windows

   [func_path, ~, ~] = fileparts(which('configure_toplevel_IO'));
   top_path = fullfile(func_path, '..', '..', '..'); % has 'dist','source', etc.
   begDir = cd(top_path);
   pd.top_path = pwd;  % Resolve '.' and '..'
   cd(begDir);
   top_path = pd.top_path;

   tmp = getenv('PROC_MODE');
   if strlength(tmp) < 1  % relevant environment variables are not set
      pd.proc_mode = 3;  % Change this to run different modes (3 and 4 supported)

      pd.ancillary_data_dir = fullfile(top_path, 'dist\ancillary');
      pd.tools_anc_data_dir = fullfile(top_path, 'PREFIRE_tools\dist\ancillary');
      pd.instrument_model_dir = fullfile(top_path, ...
                                         'dist\ancillary\instrument_model');
      pd.output_dir = fullfile(top_path, 'test\outputs');

      pd.product_fullversion = 'R01_P00';

      if (pd.proc_mode == 3)
         pd.SRF_disambig_str = 'SRF_v13_2024-09-15';

         [pd.payload_L0_fpaths{1:2}] = deal( ...
                fullfile(top_path, ...
                 'test\inputs\prefire_02_payload_tlm_20221028193650_20221029032453_20240324215847.nc'), ...
                fullfile(top_path, ...
                 'test\inputs\prefire_02_payload_tlm_20221029032450_20221029135736_20240324215903.nc'));
         pd.ref_payload_L0_fpaths = pd.payload_L0_fpaths;

         [pd.bus_L0_fpaths{1}] = deal( ...
                fullfile(top_path, ...
                 'test\inputs\prefire_02_bus_tlm_20221029071721_20221029145912_20230330143233.nc'));
         pd.ref_bus_L0_fpaths = pd.bus_L0_fpaths;

         [pd.orbit_L0_fpaths{1}] = deal( ...
                fullfile(top_path, ...
                 'test\inputs\prefire_02_orbit_reconst_20221029071721_20221029145912_20230330143233.nc'));
         pd.ref_orbit_L0_fpaths = pd.orbit_L0_fpaths;

         pd.pld_cat_fpath = fullfile(top_path, ...
                'test\inputs\prefire_02_pld_cat_20221028193650_20221029135736_20240324214512.nc');

         [pd.granule_beg_ts, pd.granule_end_ts] = deal({717502532095.872, ...
                                                     717508255640.243});
         pd.granule_ID = '00009';

      elseif (pd.proc_mode == 4)

         pd.L1A_rad_fpath = 'test\outputs\m3\raw-PREFIRE_SAT1_1A-RAD_R01_P00_20220926150106_00009.nc';
         pd.DEM_root_dir = 'C:\data\users\mmm\DEM\copernicus-dem-90m\tiles';
   
         atrack_idx_range_0bi = 'ATRACK_IDXRANGE_0BASED_INCLUSIVE:0:END';
         tmp = strsplit(atrack_idx_range_0bi, ':');
         pd.idxbeg_atrack = str2num(tmp{2})+1;  % 1-based index
         if strcmp(tmp{3}, 'END')
            pd.idxend_atrack = -1;  % Sentinel value
         else
            pd.idxend_atrack = str2num(tmp{3})+1;  % 1-based index
         end

      end

   else
      pd = obtain_from_env_vars(pd, tmp);
   end

end

 function pd = obtain_from_env_vars(pd, procmode)
    pd.top_path = getenv('PACKAGE_TOP_DIR');
    pd.proc_mode = str2double(procmode);

    pd.ancillary_data_dir = getenv('ANCILLARY_DATA_DIR');
    pd.tools_anc_data_dir = getenv('TOOLS_ANC_DIR');
    pd.instrument_model_dir = getenv('INSTRUMENT_MODEL_DIR');
    pd.output_dir = getenv('OUTPUT_DIR');

    pd.product_fullversion = getenv('PRODUCT_FULLVER');

    pd.there_is_a_pre_error = ~isempty(getenv('THERE_IS_A_PRE_ERROR'));

    if (pd.proc_mode == 3)
       pd.SRF_disambig_str = 'SRF_v13_2024-09-15';
       tmp = getenv('SRF_DISAMBIG_STR');
       if ~isempty(tmp)
          pd.SRF_disambig_str = getenv('SRF_DISAMBIG_STR');
       end

       tmp = getenv('L0_PAYLOAD_FPATHS');
       pd.ref_payload_L0_fpaths = strsplit(tmp, '||');
       tmp = getenv('PLD_TMPCUR_FPATH');
       if isempty(tmp)
          pd.payload_L0_fpaths = pd.ref_payload_L0_fpaths;
       else
          pd.payload_L0_fpaths{1} = tmp;
       end

       tmp = getenv('L0_BUS_FPATHS');
       pd.ref_bus_L0_fpaths = strsplit(tmp, '||');
       tmp = getenv('BUS_TMPCUR_FPATH');
       if isempty(tmp)
          pd.bus_L0_fpaths = pd.ref_bus_L0_fpaths;
       else
          pd.bus_L0_fpaths{1} = tmp;
       end

       tmp = getenv('L0_ORBIT_FPATHS');
       pd.ref_orbit_L0_fpaths = strsplit(tmp, '||');
       tmp = getenv('ORBIT_TMPCUR_FPATH');
       if isempty(tmp)
          pd.orbit_L0_fpaths = pd.ref_orbit_L0_fpaths;
       else
          pd.orbit_L0_fpaths{1} = tmp;
       end

       pd.pld_cat_fpath = getenv('PLD_CAT_FPATH');

       tmp = strsplit(getenv('GRANULE_START_ID_END'), '_');
       pd.granule_beg_ts = str2double(tmp{1});
       pd.granule_end_ts = str2double(tmp{3});
       pd.granule_ID = tmp{2};

       tmp = getenv('MAKE_L1A_AUX_DATA');
       if isempty(tmp)
           pd.make_L1A_aux_data = 0;
       else
           pd.make_L1A_aux_data = str2double(tmp);
       end

    elseif (pd.proc_mode == 4)
       pd.L1A_rad_fpath = getenv('L1A_RAD_FILE');
       pd.DEM_root_dir = getenv('DEM_ROOT_DIR');

       tmp = strsplit(getenv('ATRACK_IDX_RANGE_0BI'), ':');
       pd.idxbeg_atrack = str2num(tmp{2})+1;  % 1-based index
       if strcmp(tmp{3}, 'END')
          pd.idxend_atrack = -1;  % Sentinel value
       else
          pd.idxend_atrack = str2num(tmp{3})+1;  % 1-based index
       end

    end
 end

end
