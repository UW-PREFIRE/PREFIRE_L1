function [satint] = oblate_range(satint, varargin)
if(length(varargin)>0)
    a=varargin{1};
else
    a=6.378137000000000e+03; %default value
end

if(length(varargin)>1)
    esq=varargin{2};
else
    esq=0.0067; %default value
end
if(length(varargin)>2)
    iterations=varargin{3};
else
    iterations=3; %default value
end

newlat = deg2rad(satint.P); %initial guess for ground latitude on oblate spheroid is assumed sphere 

for i = 1:iterations %literature says that two-three iterations is enough to correct latitude for oblateness
    C = 1./(sqrt(1-esq*sin(deg2rad(satint.P)).*sin(newlat))); %removed an extra conversion factor 6/28/2020 in second sin
    newlat = atan((satint.Z+a*C*esq.*sin(newlat))./satint.R);
end
satint.nadir_range = (satint.R/cos(newlat)) - a*C;  %test case shows range of 22 km, thats good, but why 485 - 507 km?

end
