function [NEDR, bandwidth,detector_noise,responsivity,Dstar] = calculate_NEDR2(wl,TIRS,wavelen,SRF)
%% NED(S)R calculation (Noise-Equiv Delta Spectral Radiance)
%
% these are instrument design parameters
%
detector_noise = TIRS.noise; 
% originally
% 60 nV noise level for high gain
% 90 nV noise level for low gain
responsivity = TIRS.R; % this is an array now 
% nominally 4200 V/W, responsivity of thermopile detectors
% also requires fno, pix_size (set above)
% convert pix_size (width in um) to det_size (area, in m)
det_size = (TIRS.d*1e-6)^2;
% Noise equiv Power (W)
NEP = detector_noise ./ responsivity;
% Noise Equivalent Exitance from Scene in W/m2
NEE=(NEP.*4*TIRS.fno^2)/det_size;
% derive Dstar, just for comparison.
% bandwidth in Hz.
bandwidth = 1/TIRS.tau;
% Dstar is conventially [cm sqrt(Hz) / W], so convert pix_size to cm
Dstar = (TIRS.d*1e-4) * sqrt(bandwidth) ./ NEP;

% A. intuitive  guess for radiance noise level (supported by
% Dereniak & Boreman 14.12) - just the NEE divided by pi,
% similar to converting from the irradiance to radiance.
% I am not 100% sure why that makes sense, becaise E = pi L
% only when you integrate a L-field that is constant with angle
% and you integrate over the hemisphere.
NEDR_constant = NEE / pi;
wl_delta = wavelen(2) - wavelen(1);
SRF_norm = reshape(sum(SRF, 1) .* wl_delta, [TIRS.spectral_channels,TIRS.spatial_scenes]); % adjustment for 8 scenes in version 2
NEDR = NEDR_constant ./ SRF_norm;

end
