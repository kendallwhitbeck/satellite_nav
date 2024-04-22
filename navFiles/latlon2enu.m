function [R] = ecef2enu(lat,lon)
% ecef2enu : Generate the rotation matrix used to express a vector written in
% ECEF coordinates as a vector written in local east, north, up
% (ENU) coordinates at the position defined by geodetic latitude
% and longitude.
%
% INPUTS
% lat ----- geodetic latitude in radians
% lon ----- longitude in radians
%
% OUTPUTS
% R ------- 3-by-3 rotation matrix that maps a vector v_ecef expressed in the
% ECEF reference frame to a vector v_enu expressed in the local
% east, north, up (vertical) reference frame as follows: v_enu = R*v_ecef.

R = [-sin(lon), cos(lon), 0 ;
    -sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat);
    cos(lon)*cos(lat), sin(lon)*cos(lat), sin(lat)];

end
%+------------------------------------------------------------------------------+
% References: ASE 372N Module 2 - Moriba K. Jah, Ph.D.
%
%
% Author: Kendall Whitbeck
%+==============================================================================+