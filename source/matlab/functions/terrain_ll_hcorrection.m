function [lon_corr, lat_corr, h_corr, iter_ct, this_geo_qf] = ...
        terrain_ll_hcorrection(sat_lon, sat_lat, sat_h, ground_lon, ...
                               ground_lat, start_h, E, geoid, dem, this_geo_qf)
%terrain_ll_hcorrection apply terrain correction to a lat/lon position that is geolocated to the ellipsoid
%
% [lon_corr, lat_corr, h_corr, iter_ct] = terrain_ll_hcorrection( ...
%     sat_lon, sat_lat, sat_h, ground_lon, ground_lat, start_h, ...
%     geoid, dem)
%
% Inputs:
%
% sat_lon, sat_lat, sat_h describe the geodetic ECEF position of
% the satellite above the ellipsoid. Units are degrees (lat/lon)
% and meters (h)
%
% ground_lon, ground_lat are the geolocated positions of the
% viewing line of sight on the ellipsoid.
%
% start_h: a start altitude [m], which initializes the iterative loop
% to determine the new ground lon/lat/alt. This should be chosen
% higher than any possible terrain near to the starting
% ground_lon/ground_lat. "Near" is tricky to define: this would
% require a larger spatial domain for very large slant angles.
%
% E, the referenceEllipsoid object, in units of meters.
%
% geoid: an instance of the geoid_model object class (see
% geoid_model.m)
%
% dem: an instance of the elevation_model object class (see
% elevation_model.m).
%
%
% Outputs:
%
% lon_corr, lat_corr, h_corr: the terrain-corrected lon, lat and
% altitude values (degrees and meters, following the input
% inputs). The reported altitude is relative to the geoid.
%
% iter_ct: the number of iterations used to find the corrected
% values.
%
% Implementation details:
%
% The ellipsoid is assumed to be the WGS84.
%
% The algorithm largely follows the method
% in the VIIRS geolocation ATBD, with some modifications to switch
% to a fixed step size that seems more appropriate for near-nadir
% views.
% 
% Joint Polar Satellite System (JPSS) VIIRS Geolocation Algorithm
% Theoretical Basis Document (ATBD)
% Rev A, Effective Date March 8, 2017.
%
%

verbose = false;

[sat_x, sat_y, sat_z] = geodetic2ecef(E, sat_lat, sat_lon, sat_h);
[ground_x, ground_y, ground_z] = geodetic2ecef(E, ground_lat, ground_lon, 0.0);

% grouping into column vector, [3,1].
sat_r = [sat_x; sat_y; sat_z];
ground_r = [ground_x; ground_y; ground_z];

% vector arithmetic to get the normalized ECEF pointing vector 
% from ground point to the satellite (uu), and the pseudo-zenith angle
% from the ground point local ellipsoid normal (nu)
% (I think this is not technically the local normal, or the zenith
% angle, because it is computed directly from the geodetic angles)
coslat = cos(deg2rad(ground_lat));
% converts to shape [3,1]
ground_n = [coslat .* cos(deg2rad(ground_lon)); ...
            coslat .* sin(deg2rad(ground_lon)); ...
            sin(deg2rad(ground_lat))];

u = sat_r - ground_r;
uu = u ./ norm(u);

cos_nu = sum(uu .* ground_n, 1);
sin_nu = sqrt(1 - cos_nu.*cos_nu);

% The geoid surface is very smooth, so we make the approximation
% that it is constant in the domain where the terrain correction is
% done. This, we only use this one value.
h_geoid = geoid.lookup(ground_lon, ground_lat);

% start iterative loop, working down from a large altitude.
% (this should be modified to the max height within some region
% around scene poly.)
% notation here: LOS_h2geoid, LOS_h2ellip meaning height relative
% to geoid or ellipsoid, for a point along LOS vector (u)
% the underscore suffix denotes iterating value
% (_iter), previous iter (_prev) or final (_final).
% the terrain_h2geoid, follows the same notation, but this is for
% the terrain h sampled from the DEM.
%
% compared to the VIIRS ATBD, the LOS_h is their "h prime", and the
% terrain_h is their "h"
start_h2geoid = start_h;

% convert to geodetic lat/lon/height along LOS vector r (Line Of Sight)
% here the LOS_h is the height above ellipsoid for the point along
% the Line Of Sight (converted from ECEF vector to geodetic llh).
start_h2ellip = start_h2geoid + h_geoid;
Delta_u = start_h2ellip ./ cos_nu;
ground_r_iter = ground_r + Delta_u * uu;
[lat_iter, lon_iter, LOS_h_iter] = ecef2geodetic( ...
                      E, ground_r_iter(1), ground_r_iter(2), ground_r_iter(3));
terrain_h2geoid_iter = dem.get_nearest_neighbor_h(lon_iter, lat_iter);
terrain_h2ellip_iter = terrain_h2geoid_iter + h_geoid;

% increment size: the fixed increment ds is a horizontal scale:
% this is converted into an increment along the LOS (du) that will
% move the intersection point in the x-y plane by ds. This is the
% VIIRS ATBD method, which the rationale that we should move along
% u by an amount corresponding to the DEM grid spacing (so ds
% should be at the scale of the DEM grid)
% This has a side effect that near-nadir views could have enormous
% du lengths (since sin_nu -> 0), which seems risky.
% instead, we choose a fixed 1km increment. by the VIIRS rationale,
% this is risky since we could skip over peaks (since our DEM is
% much finer than 1 km spacing), but this only is a concern for
% large zenith angles. Since TIRS will be near nadir, that should
% not be a concern.

% VIIRS ATBD
%ds = 500.0;
%du = ds ./ sin_nu;
% instead use 1 km (note we don't really use ds, it is only set as
% the horizontal scale which in turn sets du.)
ds = 0.0;
du = 1000.0;

% Try to refine LOS_h_iter and terrain_h2ellip_iter, in order to minimize the
%  amount of iterations needed in the next code section:
t_start_h2ellip = terrain_h2geoid_iter+du*0.5 + h_geoid;
t_Delta_u = t_start_h2ellip ./ cos_nu;
t_ground_r_iter = ground_r + t_Delta_u * uu;
[t_lat_iter, t_lon_iter, t_LOS_h_iter] = ecef2geodetic( ...
                E, t_ground_r_iter(1), t_ground_r_iter(2), t_ground_r_iter(3));
t_terrain_h2geoid_iter = dem.get_nearest_neighbor_h(t_lon_iter, t_lat_iter);
t_terrain_h2ellip_iter = t_terrain_h2geoid_iter + h_geoid;
if t_LOS_h_iter > t_terrain_h2ellip_iter
   start_h2ellip = t_start_h2ellip;
   Delta_u = t_Delta_u;
   ground_r_iter = t_ground_r_iter;
   lat_iter = t_lat_iter;
   lon_iter = t_lon_iter;
   LOS_h_iter = t_LOS_h_iter;
   terrain_h2geoid_iter = t_terrain_h2geoid_iter;
   terrain_h2ellip_iter = t_terrain_h2ellip_iter;
end

if verbose
    fprintf('geoid_h: %6.1f\n', h_geoid);
    fprintf('start_h: %6.1f\n', start_h);
    fprintf('Du, ds, du: %8.1f %8.1f %8.1f\n', Delta_u, ds, du);
    fprintf('ground lon/lat:     %10.5f %10.5f\n', ...
            ground_lon, ground_lat);
    fprintf('iter  0 lon/lat/hLOS/hTer: %10.5f %10.5f %8.1f %8.1f\n', ...
            lon_iter, lat_iter, LOS_h_iter, terrain_h2ellip_iter);
end

% Repeat increments of size du until the LOS is below/at the terrain.
% The 20-iteration limit is sized to deal with up to 20*du ~= 20 km in
%  altitude change.
iter_ct = 1;
while (LOS_h_iter > terrain_h2ellip_iter) & (iter_ct <= 20)
   % save previous iteration.
   LOS_h_prev = LOS_h_iter;
   terrain_h2ellip_prev = terrain_h2ellip_iter;
   ground_r_prev = ground_r_iter;

   % update positions by moving downward by a du increment.
   ground_r_iter = ground_r_iter - du * uu;
   [lat_iter, lon_iter, LOS_h_iter] = ecef2geodetic( ...
                      E, ground_r_iter(1), ground_r_iter(2), ground_r_iter(3));
   terrain_h2geoid_iter = dem.get_nearest_neighbor_h(lon_iter, lat_iter);
   terrain_h2ellip_iter = terrain_h2geoid_iter + h_geoid;

   if verbose
      fprintf('iter %2d lon/lat/hLOS/hTer: %10.5f %10.5f %8.1f %8.1f\n', ...
              iter_ct, lon_iter, lat_iter, LOS_h_iter, terrain_h2ellip_iter);
   end

   iter_ct = iter_ct + 1;
end
if iter_ct == 21
   this_geo_qf = bitset(this_geo_qf, 1, 'uint16')
   if verbose
      fprintf('WARNING: Topographic correction is unfinished (too many iterations)\n');
      fprintf('         for lon/lat/hLOS/hTer: %10.5f %10.5f %8.1f %8.1f\n', ...
              lon_iter, lat_iter, LOS_h_iter, terrain_h2ellip_iter);
   end
end

% final result is the linear interp between last two steps
% (assuming we did not hit iter limit, in which case this is a
% failure)
% weighting, a, given by 3.3-104 in VIIRS ATBD.
a1 = LOS_h_prev - terrain_h2ellip_prev;
a2 = terrain_h2ellip_iter - terrain_h2ellip_prev - LOS_h_iter + LOS_h_prev;
a = a1/a2;

ground_r_final = a*ground_r_iter + (1-a)*ground_r_prev;
[lat_final, lon_final, LOS_h_final] = ecef2geodetic( ...
    E, ground_r_final(1), ground_r_final(2), ground_r_final(3));

terrain_h2geoid_final = dem.get_nearest_neighbor_h(lon_final, lat_final);
terrain_h2ellip_final = terrain_h2geoid_final + h_geoid;

if verbose
    fprintf('prev  hLOS, h2ellip, difference: %8.1f, %8.1f, %8.3f\n', ...
            LOS_h_prev, terrain_h2ellip_prev, ...
            LOS_h_prev - terrain_h2ellip_prev);
    fprintf('end   hLOS, h2ellip, difference: %8.1f, %8.1f, %8.3f\n', ...
            LOS_h_iter, terrain_h2ellip_iter, ...
            LOS_h_iter - terrain_h2ellip_iter);
    fprintf('linear weighting a: %9.6f\n', a);
    fprintf('final  lon/lat/hLOS:  %10.5f %10.5f %8.1f\n', ...
            lon_final, lat_final, LOS_h_final);
    fprintf('intp. hLOS, h2ellip, difference: %8.1f, %8.1f, %8.3f\n', ...
            LOS_h_final, terrain_h2ellip_final, ...
            LOS_h_final - terrain_h2ellip_final);
end

lon_corr = lon_final;
lat_corr = lat_final;
h_corr = terrain_h2geoid_final;

end
