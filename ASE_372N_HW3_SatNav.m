%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASE 372N: Satellite-Based Navigation
% LAB 3 MAIN FILE
close all; format long g; clc; clear all;
set(0, 'DefaultAxesFontSize',18, 'DefaultLineLineWidth',.75,...
    'DefaultLineMarkerSize',4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem 1

% Inputs
elMaskDeg = 0; % degrees
lat = (30+18/60+42.08823/3600)*pi/180; % rad
lon = -(097+45/60+22.71226/3600)*pi/180; % rad
alt = 164.9433; % meters

% Assigning Calendar Date/Time in UTC
Year = 2018;
Month = 10;
Day = 12;
Hour = 0;
Min = 0;
Sec = 0;

% Converting Calendar UTC time to MATLAB data-type "datenum"
n_a = datenum(Year,Month,Day,Hour,Min,Sec);

% Function "utc2gps()" converting UTC datenum to GPS Week & Seconds
[gpsWeek, gpsSec] = utc2gps(n_a);

% Receiver's position vector in ECEF frame
rRxEcef = lla2ecef(lat,lon,alt); % meters

% Retrieving satellite navigation data w/ zero ephemeris hour offset
[satdata, ~] = ...
      retrieveNavigationData(gpsWeek,gpsSec,0);
  
% Creating vector of GPS seconds
gpsSecEnd = gpsSec + 2 * 3600;
gpsSecVec = (gpsSec : 30 : gpsSecEnd)';
  
% Plotting satelitte trajectories at specific GPS epoch w/ elevation mask
[svIdVec, plotMat] = ...
    satmap(satdata,rRxEcef,elMaskDeg,gpsWeek,gpsSecVec,1);

%% Problem 2
% Part (b): using retrieveObservationData to download channel.mat data
if exist('channel.mat','file') == 2
    load('channel.mat');
    channel = channel;
else
    [channel] = retrieveObservationData(gpsWeek,gpsSec,'txau',pwd);
    channel = channel';
end

% Specifying GPS L1 CA as desired signal type
signalType = 0; % specfifying signal type desired, i.e. GPS L1 CA

% filtering out signal types other than GPS L1 CA
C = channel(channel(:, 13) == signalType, :);

% Determine increment of channel rows (drow) to get desired time inteval (dt)
dt = 30; % sec
drow = mod(C(1,4),dt); % number of rows by which to skip

% Redefines Channel matrix w/ values for every dt seconds (desired time
% interval)
jj = 1;
% C_dt = zeros(length(C(:,1)) , length(C(1,:)));
for ii = 1:length(C(:,1))
    if mod(C(ii,4),dt) == drow
        Cdt(jj,:) = C(ii,:); % Sampling data every 30 seconds
        jj = jj + 1;
    end%if
end%for

% Specifying which GPS Week the Channel epoch occurs in
gpsWeek2 = Cdt(1,3); % col 3 same for each row, i.e. all sats in same GPS Week

% Creating vector of GPS seconds
gpsSecVec2 = Cdt(:,4) + Cdt(:,5); % sec

% Creating vector of SVIDs
SvIdVec2 = Cdt(:,14);

% Creating vector of Carrier-to-Noise Ratios (C/N0)
CN0vec = Cdt(:,9); % dB-Hz

elRadVec = zeros(length(Cdt(:,1)),1); % preallocating for speed
% for-loop calculating elevation of each satellite at each time
for ii = 1:length(Cdt(:,1))
    
    % if-clause updating satdata under 2 conditions only:
    % % 1) the current time is not the same as the previous time (avoids
    % %    recalculating satdata for multiple sats at same time/epoch)
    % % 2) difference between current & initial time is a multiple of 2
    % %    hours (7200 seconds)
    if ii == 1 % avoids unreal indices on 1st iteration 
        [satdata2, ~] = retrieveNavigationData(gpsWeek2,gpsSecVec2(ii),0);
    elseif gpsSecVec2(ii) ~= gpsSecVec2(ii-1) && mod(gpsSecVec2(ii)-gpsSecVec2(1),7200) == 0
        [satdata2, ~] = retrieveNavigationData(gpsWeek2,gpsSecVec2(ii),0);
    end%if
    
    % Position of satellite in ECEF frame (meters)
    [rSvEcef, ~] = satloc(gpsWeek2, gpsSecVec2(ii), satdata2(SvIdVec2(ii)));
    
    % Azimuth and Elevation of sat wrt receiver (radians)
    [elRadVec(ii,1), ~] = satelaz(rSvEcef,rRxEcef);
    
end%for

% Plot Carrier-to-Noise Ratio vs. Elevation
figure
hold on
scatter(elRadVec*180/pi,CN0vec)
grid on
title('Carrier-to-Noise Ratio vs. Elevation for Satellites Broadcasting GPS L1 CA Signal Type Over 24 Hours')
xlabel('Satellite Elevation Angle (deg)')
ylabel('Carrier-to-Noise Ratio [C/N_0] (dB-Hz)')
hold off




