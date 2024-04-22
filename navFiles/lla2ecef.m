function [pVec] = lla2ecef(lat,lon,alt)
% lla2ecef : Convert from latitude, longitude, and altitude (geodetic with
% respect to the WGS-84 ellipsoid) to a position vector in the
% Earth-centered, Earth-fixed (ECEF) reference frame.
%
% INPUTS
% lat ----- latitude in radians
% lon ----- longitude in radians
% alt ----- altitude (height) in meters above the ellipsoid
%
% OUTPUTS
% pVec ---- 3-by-1 position coordinate vector in the ECEF reference frame,
% in meters.

%% Givens
% Flattening Factor of the Earth
f_inv = 298.257223563;
f = 1 / f_inv;

% Semi-major axis
a = 6378137.0; % meters

%% Calculations

% Square of the Earth's eccentricity
e_sq = 2*f - f^2;

R_N = a / sqrt(1 - e_sq*sin(lat)^2);

% Output
pVec = [(R_N + alt) * cos(lat) * cos(lon);
        (R_N + alt) * cos(lat) * sin(lon);
        (R_N*(1-e_sq) + alt) * sin(lat)];

end
%+------------------------------------------------------------------------------+
% References: ASE 372N Module 2 - Moriba K. Jah, Ph.D.
%
%
% Author: Kendall Whitbeck
%+==============================================================================+
