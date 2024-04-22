function [gpsWeek, gpsSec] = utc2gps(n)
% utc2gps : Convert UTC time to GPS time expressed in GPS week and GPS
% second of week.

% INPUTS
%
% n ------------ The UTC time and date expressed as a Matlab datenum. Use
% the Matlab function datenum() to generate n; use datestr()
% to render n in a standard format.

% OUTPUTS
%
% gpsWeek -------- The unambiguous GPS week number where zero corresponds
% to midnight on the evening of 5 January/morning of 6 January, 1980. By
% unambiguous is meant the full week count since the 1980 reference time
% (no rollover at 1024 weeks).
%
% gpsSec --------- The GPS time of week expressed as GPS seconds from
% midnight on Saturday.

%+------------------------------------------------------------------------------+
% References:
% GPS Signals, Measurements & Performance 2nd ed. - Misra & Enge
% getLeapSecondsUTC.m by Todd Humphreys
% Author: Kendall Whitbeck
%+==============================================================================+

% datetime array in UTC frame from datenum type
date_t = datetime(n, 'ConvertFrom', 'datenum');

% Converting MATLAB datetime to Julian date in UTC
JD_UTC = juliandate(date_t); % days

% Getting JD in UTC since GPS Epoch (0hr 6 January 1980, i.e. midnight
% between January 5 & 6, 1980)
JD_UTC_GPS = JD_UTC - 2444244.500; % days

% Getting the GPS lead time
leapSec = getLeapSecondsUTC(n); % sec

% Converting UTC time since GPS epoch in Julidan Days to total GPS time
% since GPS epoch in seconds
gpsSecTotal = leapSec + JD_UTC_GPS*24*3600; % sec

% Converting GPS seconds to GPS Weeks since GPS epoch
spw = 604800; % Seconds Per Week
gpsWeekTotal = gpsSecTotal / spw; % total GPS weeks including fraction of week
gpsWeek = floor(gpsWeekTotal); % not including fraction of week

% GPS Seconds since start of week
gpsSec = (gpsWeekTotal - gpsWeek)*spw;

end