function [n] = gps2utc(gpsWeek,gpsSec)
% gps2utc : Convert GPS time expressed in GPS week and GPS second of week to
% UTC.

% INPUTS
%
% gpsWeek ------- The unambiguous GPS week number where zero corresponds to
% midnight on the evening of 5 January/morning of 6 January,
% 1980. By unambiguous is meant the full week count since
% the 1980 reference time (no rollover at 1024 weeks).
%
% gpsSec -------- The GPS time of week expressed as GPS seconds from midnight
% on Saturday.

% OUTPUTS
%
% n -------------- The UTC time and date expressed as a Matlab datenum. Use
% the Matlab function datestr() to read n in a standard
% format.

%+------------------------------------------------------------------------------+
% References:
% getLeapSecondsGPS.m by Todd Humphreys
%
% Author: Kendall Whitbeck
%+==============================================================================+

% function determing leapSeconds given input time
leapSec = getLeapSecondsGPS(gpsWeek,gpsSec);

% Converting gpsWeek & gpsSec to total seconds since GPS epoch
spw = 604800; % Seconds Per Week
gpsSecTotal = gpsSec + gpsWeek * spw;

% Converting total GPS time in seconds to UTC time since GPS epoch in seconds
UTC_GPS_sec = gpsSecTotal - leapSec;

% UTC time since GPS epoch in Julian Days
JD_UTC_GPS = UTC_GPS_sec / 3600 / 24;

% UTC time in JD since the original JD epoch (days)
JD_UTC = JD_UTC_GPS + 2444244.500;

% converting UTC time in JD to type 'datetime'
time = datetime(JD_UTC, 'ConvertFrom', 'juliandate');

% converting UTC time in JD since original JD epoch to type 'datenum'
n = datenum(time);

end





