function new_P = expand_polyshape(P, fraction)
%expand_polyshape expand a polyshape representing a lon/lat boundary by a multiplticative factor
%
% new_P = expand_polyshape(P, fraction)
%
% the input P is a MATLAB polyshape object containing longitude as
% the first (x) axis and latitude as the second (y) axis.
% The fraction is some expansion factor, e.g. 0.1 will increase the
% size by 10%. A negative number will shrink the polyshape. The
% mean position of the polyshape is the same after expansion.
%
% the expansion is done assuming cartesian coordinates, so this
% will only work correctly for small polygons. Note also that
% regions near the poles will not work (this will simply create
% values outside the [-90, 90] range value for latitude.)
% Longitude wraps at 180 are correctly handled.
%
% The output new_P is another polyshape object, containing the
% expanded coordinates.
%
    lon = P.Vertices(:,1);
    lat = P.Vertices(:,2);
    % if contains an antimeridan wrap, unfold the longitude onto
    % positive values. longitude will be in range:
    % 0 to beyond 360 degrees, after this 'unfolding'
    if max(lon) - min(lon) > 180.0
        lon(lon < 0) = lon(lon < 0) + 360.0;
    end
    lon_center = mean(lon);
    lat_center = mean(lat);
    dlon = lon - lon_center;
    dlat = lat - lat_center;
    new_lon = (1+fraction)*dlon + lon_center;
    new_lat = (1+fraction)*dlat + lat_center;
    % refold any longitude back to [-180,180].
    % we need to check both endpoints, regardless if the longitude
    % was unfolded above: this is because the expanded poly could
    % now cross the antimeridian, when the original might not have.
    new_lon(new_lon > 180.0) = new_lon(new_lon > 180.0) - 360;
    new_lon(new_lon < -180.0) = new_lon(new_lon < -180.0) + 360;
    new_P = polyshape(new_lon, new_lat);
end