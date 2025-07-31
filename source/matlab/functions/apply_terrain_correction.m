function [ground_corrected, geo_qf] = apply_terrain_correction( ...
                                       pd, ground, sat, verbose, E, uc, geo_qf)
%apply_terrain_correction - main function to apply terrain model to geolocation
%
% Inputs
%
% pd: the MATLAB structure containing data paths (see configure_toplevel_IO)
% ground: MATLAB structure containing geolocation to WGS84
%     ellipsoid. (e.g., what we get from MATLAB built-in lookAtSpheroid)
%     This contains 4 fields:
%     P: scene latitude in degrees, shape (nscene, nframe)
%     Q: scene longitude in degrees, shape (nscene, nframe)
%     Ps latitude vertices, shape (nchannel, nscene, 5, nframe)
%     Qs longitude vertices, shape (nchannel, nscene, 5, nframe)
% Note that the vertices contain 5 elements, currently ordered as
% follows: (center, forward left, rear left, forward right, rear right.)
% The left/right sense is from facing forward long track and
% looking down at Earth. From an ascending (northward) part an
% orbit track, the ordering is: (Center, NW, SW, NE, SE)
% this is a backward "N" shape for an ascending scene, when
% plotting lat/lon in a convetional cartesian sense (+x = east, +y
% = north)
% the center vertex is identical to the scene latitude.
%
% The spectral points in Ps/Qs are identical, e.g. the scene
% footprint vertices are not wavelength dependent.
%
% sat: MATLAB structure with satellite positions. This has many
% fields, but only the following are used in this function:
%     geod_lat, geod_lon: geodetic lat/lon of the satellite nadir
%         point (also called ground track).
%     geod_alt: the geodetic altitude of the satellite (in km).
%         This is referenced to the WGS84 ellipsoid (I think - confirm.)
%
% verbose, a boolean specifying if progress info should be printed
%     to console.
% E, the referenceEllipsoid object, in units of meters.
% uc, the return from 'constants_unit_conversion'
%
% outputs:
% ground_corrected, a MATLAB structure with the following fields:
%     lat: terrain-corrected scene latitude (nscene, nframe)
%     lon: terrain-corrected scene longitude (nscene, nframe)
%     lat_vertices: terrain-corrected scene vertex latitudes
%         shape (nscene, 5, nframe)
%     lon_vertices: terrain-corrected scene vertex longitudes
%         shape (nscene, 5, nframe)
%     alt: mean altitude of scene in meters.
%         This value (and the other 3 altitude stats, and the land
%         fraction) are derived from the DEM pixels within the
%         scene polygon defined by the vertices.
%     alt_min: min altitude within scene
%     alt_max: max altitude within scene
%     alt_std: altitude std deviation within scene
%     land_fraction: fraction of land cover in the scene.


% parameters for the elevation model class.
% buffer size is the number of DEM tiles that can be loaded into
% memory at once.
dem_buffer_size = 20;
dem_debug_flag = false;

% vertex reordering
% This should probably be moved upstream, so save time (e.g., just
% create the data in this order)
% adds 1 to skip over the first (central position)
vv = [1, 2, 4, 3] + 1;

% collect input variables from the initial L1B structures.
% shapes are: (spectral, spatial, [5], along_track)
lat = ground.P;
lon = ground.Q;
vlat = ground.Ps;
vlon = ground.Qs;

sat_lat = sat.geod_lat;
sat_lon = sat.geod_lon;
sat_h = sat.geod_alt;

sat_h = sat_h * uc.km_to_m;

geoid_fpath = fullfile(pd.DEM_root_dir, 'us_nga_egm96_15.tif');
geoid = geoid_model(geoid_fpath);
dem = elevation_model( ...
    pd.DEM_root_dir, dem_buffer_size, verbose, dem_debug_flag);

[nspectral, nscenes, nvertex, nframes] = size(vlat);

% note that we are reshaping the vertex LL here, from [nscene, 5, nframe]
% in the input, to [5, nscene, nframe] in the output.
% TODO: better to just fix this upstream.
vlat_corr = zeros([5, nscenes, nframes]);
vlon_corr = zeros([5, nscenes, nframes]);
vlat0_corr = zeros([5, nscenes, nframes]);
vlon0_corr = zeros([5, nscenes, nframes]);

alt = zeros([nscenes, nframes]);
alt_min = zeros([nscenes, nframes]);
alt_max = zeros([nscenes, nframes]);
alt_std = zeros([nscenes, nframes]);
lfrac = zeros([nscenes, nframes]);

for n = 1:nframes

    if verbose
        disp(['Frame: ' num2str(n) ' of ' num2str(nframes)]);
    end
    for s = 1:8

        start_h = 12500.0;  % [m] Maximum possible geoid-referenced topo value

        % leading dimension: spectral, is somewhat redundant since it is
        % not being used (all values are constant in this axis)
        % So, just use the (1, ...) values.
        % note also the reshaping - (s,v) swapped to (v,s)
        for v = 1:5
           [vlon_corr(v,s,n), vlat_corr(v,s,n), h_corr, iter_ct, ...
               geo_qf(s,n)] = terrain_ll_hcorrection(sat_lon(n), sat_lat(n), ...
                                     sat_h(n), vlon(1,s,v,n), vlat(1,s,v,n), ...
                                          start_h, E, geoid, dem, geo_qf(s,n));
        end

        for v=2:5  % Do not waste time computing center point
           [vlon0_corr(v,s,n), vlat0_corr(v,s,n), h_corr, iter_ct, ...
               geo_qf(s,n)] = terrain_ll_hcorrection(sat_lon(n), sat_lat(n), ...
                           sat_h(n), ground.Qi(1,s,v,n), ground.Pi(1,s,v,n), ...
                                           start_h, E, geoid, dem, geo_qf(s,n));
        end

        P = polyshape(vlon_corr(vv,s,n), vlat_corr(vv,s,n), 'Simplify', false);
        [h_vals, m_vals, wts] = dem.get_pixel_group_h(P);
        alt(s,n) = sum(h_vals .* wts);
        alt_min(s,n) = min(h_vals);
        alt_max(s,n) = max(h_vals);
        alt_std(s,n) = std(h_vals);
        lfrac(s,n) = sum(m_vals .* wts);

    end
end

% remember that the upstream code stores 5 vertices, where
% the first is a copy of the center value (which was done above)
ground_corrected.lat = squeeze(vlat_corr(1,:,:));
ground_corrected.lon = squeeze(vlon_corr(1,:,:));
ground_corrected.lat_vertices = vlat_corr;
ground_corrected.lon_vertices = vlon_corr;

ground_corrected.lat_maxintgz_verts = vlat0_corr;
ground_corrected.lon_maxintgz_verts = vlon0_corr;

ground_corrected.alt = alt;
ground_corrected.alt_min = alt_min;
ground_corrected.alt_max = alt_max;
ground_corrected.alt_std = alt_std;
ground_corrected.land_fraction = lfrac;

end
