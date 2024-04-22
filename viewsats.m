%% viewsats script
% Given inputs, plots AzEl maps of satellite positions
clc; close all; format long g; clear all
navConstants;

% Inputs
gpsWeek = 1712; % weeks
gpsSec = 143800; % seconds
elMaskDeg = 15; % degrees
cdtRx = 0; % sec
% lat = 30.2875*pi/180; % rad
% lon = -97.7357*pi/180; % rad
% alt = 164.9433; % meters

% Receiver's position vector in ECEF frame
% % % % rRxEcef = lla2ecef(lat,lon,alt); % meters
rRx = getAntLoc('WRW0'); % meters

% Retrieving satellite navigation data w/ zero ephemeris hour offset
[satdata, ionodata] = ...
      retrieveNavigationData(gpsWeek,gpsSec,0);
  
% Creating vector of GPS seconds
gpsSecEnd = gpsSec + 0; % only using one point in time for now
gpsSecVec = (gpsSec : 1 : gpsSecEnd);

% Plotting satelitte positions at specific GPS epoch w/ elevation mask
[svIdVec, plotMat] = ...
    satmap(satdata,rRx,elMaskDeg,gpsWeek,gpsSecVec,0);

for ii = 1:length(svIdVec)
% Generate Pseudorange models for sats in view
[rhoBarVec(ii,1),Hrho(ii,:),tBar_OWLT(ii,1),T(ii,1),Ip(ii,1),...
    rSvEcefRx(ii,:),dtrel(ii,1),dtSvTx(ii,1)] =...
    satpr(gpsWeek, gpsSec, cdtRx, rRx, satdata(svIdVec(ii)), ionodata,[1;1]);
end%for

rhoBarVec
format short g
Hrho

format long g
% Test values
lla2 = ecef2lla( rRx' );
lat2 = lla2(1)*pi/180
lon2 = lla2(2)*pi/180
alt2 = lla2(3)
Ip_SV2 = Ip(1)
T_SV2 = T(1)
c_tOWLT_SV2 = tBar_OWLT(1)*cLight
r_SV2_Ecef_atRx = rSvEcefRx(1,:)
c_dtrel = dtrel(1)*cLight
c_dts = dtSvTx(1)*cLight

%% satpr
% clear tBar_OWLT; clear T;
% clear Ip
% clear rSvEcefRx
% clear dtrel
% clear dtSvTx
% 
% for ii = 1:8
% % Generate Pseudorange models for sats in view
% [rhoBarVec(ii,1),Hrho(ii,:),tBar_OWLT(ii,1),T(ii,1),Ip(ii,1),...
%     rSvEcefRx(ii,:),dtrel(ii,1),dtSvTx(ii,1)] =...
%     satpr_sarah(gpsWeek, gpsSec, cdtRx, rRx, satdata(svIdVec(ii)), ionodata,[1;1]);
% end%for
% 
% rhoBarVec
% format short g
% Hrho
% 
% format long g
% % Test values
% lla2 = ecef2lla( rRx' );
% lat2 = lla2(1)*pi/180
% lon2 = lla2(2)*pi/180
% alt2 = lla2(3)
% Ip_SV2 = Ip(1)
% T_SV2 = T(1)
% c_tOWLT_SV2 = tBar_OWLT(1)*cLight
% r_SV2_Ecef_atRx = rSvEcefRx(1,:)
% c_dtrel = dtrel(1)*cLight
% c_dts = dtSvTx(1)*cLight


%% Replotting sats w/ current GPS time
% 
% % Generating datenum-type, current UTC time
% [n] = nowUtc();
% 
% % Converting UTC datenum to GPS Week and Seconds
% [gpsWeekNow, gpsSecNow] = utc2gps(n);
% 
% % Retrieving satellite navigation data w/ zero ephemeris hour offset
% [satdata, ionodata] = ...
%       retrieveNavigationData(gpsWeekNow,gpsSecNow,-3);
% 
%   % Plotting satelitte positions at specific GPS epoch w/ elevation mask
% [svIdVec, plotMat] = ...
%     satmap(satdata,rRxEcef,elMaskDeg,gpsWeekNow,gpsSecNow,1);
% 
% 
% 
% CORS 'tbon', 3m acc, 1.25m precision, update initial guess%
