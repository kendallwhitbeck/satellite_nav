%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASE 372N: Satellite-Based Navigation
% HOMEWORK 2 MAIN FILE
close all; format long g; clc; %clear all;
set(0, 'DefaultAxesFontSize',18, 'DefaultLineLineWidth',3,...
    'DefaultLineMarkerSize',16)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Physical constants used for navigation
navConstants;

%% Part (a): GPS State Vectors for September 1, 2013 at 12:00:00 UTC

% Assigning Calendar Date/Time in UTC
Y_a = 2013;
M_a = 9;
D_a = 1;
H_a = 12;
MN_a = 0;
S_a = 0;

% Converting Calendar UTC time to MATLAB data-type "datenum"
n_a = datenum(Y_a,M_a,D_a,H_a,MN_a,S_a);

% Function "utc2gps()" converting UTC datenum to GPS Week & Seconds
[gpsWeek_a, gpsSec_a] = utc2gps(n_a);

% Checking if gps2utc returns same time & date as input in utc2gps
[ncheck_a] = gps2utc(gpsWeek_a,gpsSec_a);

% Retrieving the satellite data for the desired GPS time epoch
[satdata_a, ionodata_a] = retrieveNavigationData(gpsWeek_a,gpsSec_a,0);

% Using satloc function to compute Location (m) & Velocity (m/s) of SV#2
% at desired time expressed in the ECEF reference frame
[rEcef_SV2_a, vEcef_SV2_a] = satloc(gpsWeek_a, gpsSec_a, satdata_a(2));

% Using satloc function to compute Location (m) & Velocity (m/s) of SV#5
% at desired time expressed in the ECEF reference frame
[rEcef_SV5_a, vEcef_SV5_a] = satloc(gpsWeek_a, gpsSec_a, satdata_a(5));

% Displaying SV#2 State at 12:00:00
disp('SV#2 Location (m) & Velocity (m/s) at 12:00:00')
disp([rEcef_SV2_a, vEcef_SV2_a])

% Displaying SV#5 State at 12:00:00
disp('SV#5 Location (m) & Velocity (m/s) at 12:00:00')
disp([rEcef_SV5_a, vEcef_SV5_a])

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
%%  Part (b.): GPS State Vectors for September 1, 2013 at 12:00:01 UTC

% Assigning Calendar Date/Time in UTC
Y_b = 2013;
M_b = 9;
D_b = 1;
H_b = 12;
MN_b = 0;
S_b = 1;

% Converting Calendar UTC time to MATLAB data-type "datenum"
n_b = datenum(Y_b,M_b,D_b,H_b,MN_b,S_b);

% Function "utc2gps()" converting UTC datenum to GPS Week & Seconds
[gpsWeek_b, gpsSec_b] = utc2gps(n_b);

% Checking if gps2utc returns same time & date as input in utc2gps
[ncheck_b] = gps2utc(gpsWeek_b,gpsSec_b);

% Retrieving the satellite data for the desired GPS time epoch
[satdata_b, ionodata_b] = retrieveNavigationData(gpsWeek_b,gpsSec_b,0);

% Using satloc function to compute Location (m) & Velocity (m/s) of SV#2
% at desired time expressed in the ECEF reference frame
[rEcef_SV2_b, vEcef_SV2_b] = satloc(gpsWeek_b, gpsSec_b, satdata_b(2));

% Using satloc function to compute Location (m) & Velocity (m/s) of SV#5
% at desired time expressed in the ECEF reference frame
[rEcef_SV5_b, vEcef_SV5_b] = satloc(gpsWeek_b, gpsSec_b, satdata_b(5));

% Displaying SV#2 State at 12:00:01
disp('SV#2 Location (m) & Velocity (m/s) at 12:00:01')
disp([rEcef_SV2_b, vEcef_SV2_b])

% Displaying SV#5 State at 12:00:01
disp('SV#5 Location (m) & Velocity (m/s) at 12:00:01')
disp([rEcef_SV5_b, vEcef_SV5_b])

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
%% Part (c): Showing Position Change Consistent with Velocity

% Change in Time
dt = 1; % sec

% Position of SV#2 & SV#5 after time, dt
rEcef_SV2_a_new = rEcef_SV2_a + vEcef_SV2_a * dt;
rEcef_SV5_a_new = rEcef_SV5_a + vEcef_SV5_a * dt;

% Checking convergence for SV#2
diff_SV2 = abs(rEcef_SV2_a_new - rEcef_SV2_b);
conv_SV2 = abs(floor(log10(diff_SV2)));
fprintf('\nSV#2 x-Component\nChange in position converges w/ velocity to:\n     %d order(s) of magnitude\n',conv_SV2(1))
fprintf('\nSV#2 y-Component\nChange in position converges w/ velocity to:\n     %d order(s) of magnitude\n',conv_SV2(2))
fprintf('\nSV#2 z-Component\nChange in position converges w/ velocity to:\n     %d order(s) of magnitude\n\n',conv_SV2(3))

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
% Checking convergence for SV#5
diff_SV5 = abs(rEcef_SV5_a_new - rEcef_SV5_b);
conv_SV5 = abs(floor(log10(diff_SV5)));
fprintf('\nSV#5 x-Component\nChange in position converges w/ velocity to:\n     %d order(s) of magnitude\n',conv_SV5(1))
fprintf('\nSV#5 y-Component\nChange in position converges w/ velocity to:\n     %d order(s) of magnitude\n',conv_SV5(2))
fprintf('\nSV#5 z-Component\nChange in position converges w/ velocity to:\n     %d order(s) of magnitude\n',conv_SV5(3))


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Testing SatLoc Function
% ~~ DISREGARD THIS SECTION ~~ %
% clc
navConstants;

% GPS Time Test Values
gpsWeekTest = 1653;
gpsSecTest = 570957.338101566;

% Satellite Data Test Values
[satdataTest, ionodataTest] = retrieveNavigationData(gpsWeekTest,gpsSecTest,0);

% Position & Velocity Test Values
[rSvEcefTest1, vSvEcefTest1] = satloc(gpsWeekTest, gpsSecTest, satdataTest(1));
[rSvEcefTest5, vSvEcefTest5] = satloc(gpsWeekTest, gpsSecTest, satdataTest(5));
[rSvEcefTest31, vSvEcefTest31] = satloc(gpsWeekTest, gpsSecTest, satdataTest(31));

% Checking the calendar time
[UTC_since_GPS_Test] = gps2utc(gpsWeekTest,gpsSecTest); % datenum
timeTest = datestr(UTC_since_GPS_Test);