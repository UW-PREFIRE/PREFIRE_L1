function [outp_wl_fpath] = output_wavelengths_2(slit_width, pix_width, ...
                   num_spectral, num_spatial, x_slit, y_slit, tau_slit, ...
                   x_fp, y_fp, z_fp, tau_fp, tip_fp, tilt_fp, dd, ...
                   wavelengths, TIRS_number);
% angles are in degrees

wl.slit_width = slit_width;
wl.pix_width = pix_width;
wl.num_spectral = num_spectral;
wl.num_spatial = num_spatial;
wl.x_slit = x_slit;
wl.y_slit = y_slit;
wl.tau_slit = tau_slit;
wl.x_fp = x_fp;
wl.y_fp = y_fp;
wl.z_fp = z_fp;
wl.tau_fp = tau_fp;
wl.tip_fp = tip_fp;
wl.tilt_fp = tilt_fp;
wl.dd = dd;
wl.wavelengths = wavelengths;

if TIRS_number > 0
   outp_wl_fpath = sprintf('wavelengths_instrument%1d.mat', TIRS_number);
else
   outp_wl_fpath = 'wavelengths_ideal.mat';
end
save(outp_wl_fpath, '-struct', 'wl');

end
