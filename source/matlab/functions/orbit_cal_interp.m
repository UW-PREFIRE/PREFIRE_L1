function [error_ca, calseq_gain, gain_at_obs, space_at_obs, valid_cal] = ...
        orbit_cal_interp(L0, calseq, L1, SRFMC, SRFMS, T_cal, TIRS)
% Calculate calibration fields interpolated to the ctime of the 'obstgt' frames.

error_ca = {'#NONE#', 0};  % Default

% Note on methods:
% we compute the gain at the cal sequence times, using the raw
% signal from the calibrator and space views. The temporal
% interpolation to the earth observation times (the obstgt frames)
% is done on the gain and the offset (which is equal to the space
% view). This now also returns a (allspectral, xtrack) array,
% valid_cal, which specifies which channels had valid cal data

%  Determine temperature lookup-table indices and associated interpolation
%  fractions for the calibrator target values
calseq_T = calseq.cal_seqmean_temp;
[cminT, cmaxT] = bounds(calseq_T);
if (cmaxT > T_cal(end))
   tmp_str = sprintf(['ERROR: the internal calibration target temperature '...
                      '(%.2fK) is > %.2fK (the maximum lookup table ' ...
                      'temperature)'], cmaxT, T_cal(end));
   gain_at_obs = -999;  % Allow error code propagation (avoid "not assigned")
   space_at_obs = -999;  % Allow error code propagation (avoid "not assigned")
   error_ca = {tmp_str, 30};
   return
elseif (cminT < T_cal(1))
   tmp_str = sprintf(['ERROR: the internal calibration target temperature ' ...
                      '(%.2fK) is < %.2fK (the minimum lookup table ' ...
                      'temperature)'], cminT, T_cal(1));
   error_ca = {tmp_str, 31};
   gain_at_obs = -999;  % Allow error code propagation (avoid "not assigned")
   space_at_obs = -999;  % Allow error code propagation (avoid "not assigned")
   return
end

% find linear interpolation to precomputed T grid (T_cal)
% We assume the grid spacing is 1K here.
calseq_dT = calseq_T - T_cal(1);  % [K]
fraci = mod(calseq_dT, 1);
i_cT = fix(calseq_dT)+2;
calseq_gain = zeros(calseq.n, L0.pdims.allspectral, L0.pdims.xtrack);

% perform the linear interp on the fly, to compute the gain for
% each cal sequence.
for i=1:L0.pdims.xtrack
   for j=1:L0.pdims.allspectral
      srfmc_intrp = (1-fraci') .* squeeze(SRFMC(j,i,i_cT-1)) + ...
                       fraci'  .* squeeze(SRFMC(j,i,i_cT));
      calseq_gain(:,j,i) = calseq.contrast_seqmean_v(:,j,i) ./ (srfmc_intrp - SRFMS(j,i));
   end
end

% using Modified Akima piecewise cubic Hermite interpolation (makima)
% interpolate the gain and offset (space mean) to the earth obs times.
ic = 0;
valid_cal = logical(zeros(L0.pdims.allspectral, L0.pdims.xtrack));
for i=1:L0.pdims.xtrack
   for j=1:L0.pdims.allspectral
      ic = ic+1;
      % catch problem values: if the calseq gain was NaN (occurs
      % for the shortest wavelength channels where there is no
      % significant radiance from the calibrator), or zero (occurs
      % when the detector is dead, and reports the same DN for
      % space and cal views), then the cal radiance will end up
      % being NaN.
      if any(isnan(calseq_gain(:,j,i))) | any(calseq_gain(:,j,i)==0)
          valid_cal(j,i) = false;
      else
          g_pp{ic} = makima(calseq.contrast_seqmean_ts(:), ...
                            calseq_gain(:,j,i));
          s_pp{ic} = makima(calseq.space_seqmean_ts(:), ...
                            calseq.space_seqmean_v(:,j,i));
          valid_cal(j,i) = true;
      end
   end
end

gain_at_obs = zeros(L0.pdims.allspectral, L0.pdims.xtrack, ...
                    length(L0.obstgt.ctime));
space_at_obs =  zeros(L0.pdims.allspectral, L0.pdims.xtrack, ...
                      length(L0.obstgt.ctime));
ic = 0;
for i=1:L0.pdims.xtrack
   for j=1:L0.pdims.allspectral
      ic = ic+1;
      x_ts = L1.ctime;  % [s] ctime of obstgt images
      if valid_cal(j,i)
          gain_at_obs(j,i,:) = ppval(g_pp{ic}, x_ts);
          space_at_obs(j,i,:) = ppval(s_pp{ic}, x_ts);
      end
   end
end

end
