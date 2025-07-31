function [outp_gm_fpath] = output_geometry_2(slit_width, pix_width, ...
                   num_spectral, num_spatial, x_slit, y_slit, tau_slit, ...
                   x_fp, y_fp, z_fp, tau_fp, tip_fp, tilt_fp, pixels, slit, ...
                   IFOVpixels, TIRS_number);
%% 10/17/22 rewrite to use matlab output in structure form, remove offsets
%% angles are in degrees !!
geometry.x_slit = x_slit;
geometry.y_slit = y_slit;
geometry.tau_slit  = tau_slit;
geometry.x_fp = x_fp;
geometry.y_fp = y_fp;
geometry.z_fp = z_fp;
geometry.tau_fp = tau_fp;
geometry.tip_fp = tip_fp;
geometry.tilt_fp = tilt_fp;
%% four dimensional pixel data 
% index 1 axes label x=1, y = 2
% index 2 spectral channel 1 through 64
% index 3 spatial channel 1 through 8
% index 4 center/corners -+, ++, --, +-


geometry.slit = slit;
geometry.slit_width = slit_width;
geometry.pix_width = pix_width;
geometry.num_spectral = num_spectral;
geometry.num_spatial = num_spatial;
geometry.xy_p = pixels;
geometry.IFOVxy_p = IFOVpixels;
geometry.slit = slit;

if TIRS_number > 0
   outp_gm_fpath = sprintf('geometry_instrument%1d.mat', TIRS_number);
else
   outp_gm_fpath = 'geometry_ideal.mat';
end
save(outp_gm_fpath, '-struct', 'geometry');

end
