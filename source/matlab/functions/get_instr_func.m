function [instr_func] = get_instr_func( ...
    wavelen, center_wave, slit_width, pix_size, pix_wave, fno);
%
% function instr_func = get_instr_func( ...
%    wavelen, center_wave, slit_width, pix_size, pix_wave, fno);
%
% computes an instrument function:
% convolution of a sinc^2 function (diffraction) with a rectangle
% (for the slit).
%
% inputs:
% wavelen: the array of wavelengths [um] where the instrument function
%   needs to be defined and calculated
% center_wave: center wavelength for the instrument function [um].
%   could consider this to be the (monochromatic) wavelength
%   incident on the detector.
% slit_width: width of slit [um]
% pix_size: size of pixel [um] in spectral domain
% pix_wave: pixel wavelength range [um], a 2-element array
% fno: f-number of optical beam
%
% proceedure:
% I assume a linear mapping between the wavelength of the dispersed
% light incident on the detector. This means that:
% (pix_wave(2)-pix_wave(1))/pix_size gives the linear relationship
% between the dispersed wavelength and the linear coordinate on
% the detector pixel.
% this allows me to map the diffraction pattern (which is computed
% in terms of the linear x-coordinate on the detector plane) to the
% spectral wavelength domain.


% convert the input wavelen grid (in spectral um)
% to linear positions at the detector array
pix_wave_delta = pix_wave(2) - pix_wave(1);
scaling = pix_size / pix_wave_delta;

x = wavelen * scaling;
center_x = center_wave * scaling;
xr = x - center_x;

% get diffraction pattern for the central wavelength.
gamma = 1.22 * fno * center_wave;
diff_pattern = sinc(xr/gamma).^2 / gamma;

% normalize the diff pattern, since this will be used as a
% convolution kernel over the slit function.
diff_pattern = diff_pattern ./ sum(diff_pattern);

% make a rectangle, with the width of the slit, to convolve with
% the diff pattern.
dx = x(2) - x(1);
slit_x = (-ceil(slit_width/dx)*dx:dx:ceil(slit_width/dx)*dx)';
slit_func = rectpuls(slit_x / slit_width);

instr_func = conv(diff_pattern, slit_func, 'same');

% Note: I now think it is not appropriate to renormalize this to
% have a unit integral. Instead, I think it is more correct to
% think about this term as a "throughput", so the peak amplitude
% matters. So note above that the slit_func (the rect function) is
% not renormalized (it has a peak value of 1), but we normalize the
% diff_pattern since that is a convolution kernel.
% So, for now, skip this normalization step at the end.

% this is normalized for the x-coord, but we need it to be norm in
% the wavelength coord. 
%instr_func = instr_func / trapz(wavelen, instr_func);

end
