function [SRFMS, SRFMC, NEDR, dtbf, T_cal, gm, TIRS, ARF, sensor_d, ...
          factor] = load_geom_SRFM(pd, TIRS_ID, uc, ti, L0)

gg_att_C = netcdf.getConstant('NC_GLOBAL');

if (TIRS_ID == 0)
   fpath_to_load = fullfile(pd.ancillary_data_dir, 'BRF', ...
                            'PREFIRE_BRF_v0.06_2022-11-15_ideal.mat');
   load(fpath_to_load);

   fpath_to_load = fullfile(pd.instrument_model_dir, 'geometry_ideal.mat');
   [gm] = load(fpath_to_load);
elseif (TIRS_ID == 1) | (TIRS_ID == 2)
    instr_model_fpath = fullfile(pd.ancillary_data_dir, 'SRF', ...
               sprintf('PREFIRE_TIRS%1d_%s.nc', TIRS_ID, pd.SRF_disambig_str));

    ncid = netcdf.open(instr_model_fpath, 'NC_NOWRITE');
    NEDR.im_version = netcdf.getAtt(ncid, gg_att_C, 'instrument_model_version');
    netcdf.close(ncid)  % Done with input SRF file

    NEDR.channel_mean_wavelen = ncread(instr_model_fpath, ...
                                       '/channel_mean_wavelen')';
    NEDR.channel_center_wavelen = ncread(instr_model_fpath, ...
                                         '/channel_center_wavelen')';

    % implementation note: in the instrument model file, the
    % channel 0 information is stored separately. for various
    % reasons, it is simpler to combine the precomputed radiance
    % together to a single  array the same size as the full FPA.
    % (the SRFMC array) radiances will be "split" back into the
    % channel 0 and spectral channels for output.
    % we don't do this to the detector flats

    T_cal = ncread(instr_model_fpath, '/T_grid')';
    tmp = permute(ncread(instr_model_fpath, '/rad'), [2,1,3]);
    tmp0 = ncread(instr_model_fpath, '/channel_0_rad');

    tmp0 = reshape(tmp0, [1,1,length(T_cal)]);
    SRFMC = [repmat(tmp0, [1,L0.pdims.xtrack,1]); tmp];  % Add channel "0"

    % Note that MATLAB converts FillValue to NaN - so we need to
    % change that back to the fill value.
    tmp = ncread(instr_model_fpath, '/NEDR')';
    tmp = fillmissing(tmp, 'constant', -9.999e3);
    tmp0 = ncread(instr_model_fpath, '/channel_0_NEDR')';
    NEDR.data = [tmp0; tmp];  % Add channel "0"

    tmp = ncread(instr_model_fpath, '/detector_bitflags')';
    tmp0 = ncread(instr_model_fpath, '/channel_0_bitflags')';
    dtbf = [tmp0; tmp];  % Add channel "0"

    SRFMS = zeros(L0.pdims.allspectral, L0.pdims.xtrack);  % Has channel "0"

    fpath_to_load = fullfile(pd.instrument_model_dir, ...
        sprintf('geometry_instrument%1d.mat', TIRS_ID));
    [gm] = load(fpath_to_load);
end

[TIRS] = prefire_TIRS(ti.ROIC_tau, pd, 2, 1, TIRS_ID);

sensor_d = TIRS.d*uc.um_to_m;  % [m] width and length of a single TIRS detector

% Left edges of each projected pixel, averaged
edge_left = (gm.xy_p(1,:,:,2)+gm.xy_p(1,:,:,4))*0.5;
% Right edges of each projected pixel, averaged
edge_right = (gm.xy_p(1,:,:,3)+gm.xy_p(1,:,:,5))*0.5;

% Calculate the radiometric correction due to defocused pixels (this was
%  explicitly excluded in the SRF output after 2020-04-24).  This ARF matrix is
%  unity for an ideal alignment.  This is currently only one component of the
%  total ARF, which also will contain slit broadening and motion blur.
ARF = reshape(abs(edge_left-edge_right)/gm.pix_width, ...
              [L0.pdims.allspectral, L0.pdims.xtrack]);

% AJM's note: the application of ARF and factor is likely incorrect,
% as this scales only the earth view data, not the calibration view
% data. Therefore it directly biases the earth view data relative
% to the calibrator. In theory, ARF might have some role in
% correcting per-detector biases - it is shaped (spectral x xtrack),
% but right now ARF is equal to 1.
% For now, simply skip the application of the seeming arbitrary
% 'factor' equal to 2, and let the ARF proceed through L1B, as it
% does not modify the radiance right now. (it is equal to 1.0)

factor = 1.0;
%factor = 2.0;  % This scales the emission and skews the brightness temperature
ARF = factor*ARF;  % This doubling is consistent with the slit width

end
