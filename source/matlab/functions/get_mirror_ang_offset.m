function [ypr] = get_mirror_ang_offset(ypr, L1A, nframes)

% A movable-mirror angular offset is incorporated into the geolocation
%  calculations as a term added to the pitch angle.

% The movable-mirror "static" angular offset ('L1A.mirror_ang_offset0_deg';
%  with respect to the nadir aperture boresight) is determined in
%  `categorize_L0_payload_tlm`, communicated via the 1A-RAD granule file, and
%  then used here.
ypr(2,:) = ypr(2,:)+L1A.mirror_ang_offset0_deg;  % static pitch adjustment [deg]

% Similarly, the movable-mirror time-variable angular offset array
%  ('L1A.mirror_ang_offset'; with respect to the best available estimate of
%  the nadir aperture boresight), is communicated via the 1A-RAD granule file,
%  and then used here to determine the appropriate mirror angular offset value
%  for each frame.

mangoffs_atr = length(L1A.mirror_aostart_ctime);
i_maoa = 1;
if (L1A.mirror_aostart_ctime(i_maoa) > L1A.ctime(1)) & ...
        (L1A.mirror_aostart_ctime(mangoffs_atr) > L1A.ctime(nframes))
%@   ib = 1;
%@   ie = nframes;
%@   ypr(2,ib:ie) = ypr(2,ib:ie)+0.;
   disp('WARNING: no mirror info for any frame(s), using default value');
   return  % Nothing more to do
end

  % Determine first frame index with mirror info: 
if L1A.mirror_aostart_ctime(i_maoa) > L1A.ctime(1) % No mirror info @ first frame
   ib = find(L1A.ctime >= L1A.mirror_aostart_ctime(i_maoa), 1);
%@   ie = max(1, ib-1);
%@   ypr(2,1:ie) = ypr(2,1:ie)+0.;
   disp('WARNING: no mirror info for first frame(s), using default value for these'); 
else
   ib = 1;
end

  % There is mirror info for the first frame (at least).  Find which mirror info
  %  element to start using:
tmp_i_maoa = find(L1A.mirror_aostart_ctime > L1A.ctime(1), 1);
if isempty(tmp_i_maoa)
   ie = nframes;
   i_maoa = mangoffs_atr;
   ypr(2,ib:ie) = ypr(2,ib:ie)+L1A.mirror_ang_offset(i_maoa);
   disp('WARNING: assuming last mirror info element is valid for remaining frames');
   return  % Nothing more to do
else
   i_maoa = max(1, tmp_i_maoa-1);
end

  % Determine the end ctime index for this chunk of frames, then apply the
  %  appropriate mirror angle to that chunk:
if mangoffs_atr == 1
   ie = nframes;
else
   for i=tmp_i_maoa:mangoffs_atr
      ii = find(L1A.ctime >= L1A.mirror_aostart_ctime(i), 1);
      if isempty(ii)
         ie = nframes;
         break;
      else
         ie = ii-1;
      end
      ypr(2,ib:ie) = ypr(2,ib:ie)+L1A.mirror_ang_offset(i_maoa);
      i_maoa = i;  % For next iteration
      ib = ie+1;   %
   end
end
if ib > ie
   ie = nframes;  % need to finish setting values for this chunk's frames
end
ypr(2,ib:ie) = ypr(2,ib:ie)+L1A.mirror_ang_offset(i_maoa);

end
