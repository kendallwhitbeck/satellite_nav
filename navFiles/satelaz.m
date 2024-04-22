%%
function [elRad, azRad] = satelaz(rSvEcef,rRxEcef)
% satelaz : Compute satellite elevation and azimuth angles in radians with
% respect to receiver location.
%
%
% INPUTS
%
% rSvEcef ---- 3-by-1 satellite location in ECEF coordinates, in meters.
%
% rRxEcef ---- 3-by-1 receiver location in ECEF coordinates, in meters.
%
%
% OUTPUTS
%
% elRad ------ Satellite elevation angle (the angle between the WGS-84 local
% ENU tangent plane and the receiver-to-satellite vector),
% in radians.
%
% azRad ------ Satellite azimuth angle (the angle in the ENU tangent plane
% between North and the receiver-to-satellite vector projection,
% measured positive clockwise), in radians.
%
%+------------------------------------------------------------------------------+
% References:
% https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
%
% Author: Kendall Whitbeck
%+==============================================================================+

% Latitude (rad) & Longitude (rad) Values from ECEF receiver position
lla = ecef2lla(rRxEcef');
latRx = lla(1)*pi/180;
lonRx = lla(2)*pi/180;

% Transformation matrix from ECEF to East-North-Up (ENU) frame, calculated
% using latitude & longitude of receiver
[T_ECEF_ENU] = latlon2enu(latRx,lonRx);

% Position of Sat & Receiver in ENU frame (meters)
rSv_ENU = T_ECEF_ENU * rSvEcef;
rRx_ENU = T_ECEF_ENU * rRxEcef;

% Line of Sight (relative position) unit vector btwn Receiver & Satellite
% in ECEF frame
rho_ENU = (rSv_ENU-rRx_ENU) / norm(rSv_ENU-rRx_ENU); % unit vector

% Elevation of satellite wrt receiver
elRad = asin(rho_ENU(3)); % rad

% Azimuth of satellite wrt receiver
azRad = atan2(rho_ENU(1),rho_ENU(2)); % rad

end

