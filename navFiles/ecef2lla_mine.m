function [lat,lon,alt] = ecef2lla_mine(rRxEcef)
% ecef2lla : Convert from a position vector in the Earth-centered, Earth-fixed
% (ECEF) reference frame to latitude, longitude, and altitude
% (geodetic with respect to the WGS-84 ellipsoid).
%
% INPUTS
% rRxEcef ---- 3-by-1 receiver position coordinate vector in the ECEF
% reference frame, in meters.
%
% OUTPUTS
% lat ----- latitude in radians
% lon ----- longitude in radians
% alt ----- altitude (height) in meters above the ellipsoid

%% Givens
% pVec components
x = rRxEcef(1);
y = rRxEcef(2);
z = rRxEcef(3);

% Flattening Factor of the Earth
f_inv = 298.257223563;
f = 1 / f_inv;

% Square of the Earth's eccentricity
e_sq = 2*f - f^2;

% Semi-major axis
a = 6378137.0; % meters

%% Calculations
lon = atan2(y,x);

rho = sqrt(x^2 + y^2);

r = sqrt(x^2 + y^2 + z^2);

% Initial Guess
lat = asin(z / r);

% Predefining tolerance & flag for latitude convergence
tol = 1e-8;
flag = 0;
% Iterative Loop for finding accurate latitude value
while flag ~= 1
    
    R_N = a / sqrt(1 - e_sq * sin(lat)^2);
    
    lat_new = atan2(z + R_N*e_sq*sin(lat), rho);
    
    if abs(lat_new - lat) <= tol
        flag = 1;
    end
    
    % Assigning more accurate latitude value to "lat"
    lat = lat_new;
end

% Updates value of R_N using more accurate latitude value
R_N = a / sqrt(1 - e_sq * sin(lat)^2);

% Altitude
alt = rho / lat - R_N;
end
%+------------------------------------------------------------------------+
% References: ASE 372N Module 2 - Moriba K. Jah, Ph.D.
%
%
% Author: Kendall Whitbeck
%+========================================================================+