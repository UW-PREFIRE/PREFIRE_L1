function sci_DN = srd_to_sci_DN(TIRS, srd_DN)
%
% function sci_DN = srd_to_sci_DN(srd_DN)
%
% convert srd_DN (raw DN in "SRD" convention, as they are stored in
% raw L0 packet data) into sci_DN (raw DN in "SCI" convention).
%
% SRD convention is expected to be an array with size (512,nframe);
% SCI convention will be an array with size (64,8,nframe).
%
% SCI convention includes alternating polarity shift, reshaping
% from 512 to 64,8, and "ROIC unmapping".

nframes = size(srd_DN,2);
srd_tmp = zeros([TIRS.spectral_channels, ...
                 TIRS.spatial_scenes, ...
                 nframes]);

for i=1:TIRS.spatial_scenes
    for j=1:TIRS.spectral_channels
        k = TIRS.spatial_scenes*(j-1) + i;
        o = TIRS.o * 2 * mod(j,2);
        srd_tmp(j,i,:) = ((-1)^(j)) .* srd_DN(k,:) + o;
    end
end

sci_DN = ROIC_map_io( ...
    TIRS.spectral_channels, ...
    TIRS.spatial_scenes, ...
    nframes, ...
    srd_tmp);

end
