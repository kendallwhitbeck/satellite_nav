function [rhoBar,Hrho,t_OWLT,T,Ip,rSatEcefTx,dtrel,dtSvTx,dXspin] = ...
    satpr(gpsWeekRx, gpsSecRx, cdtRx, rRxEcef, sd, ionodata, tiFlags)
% satpr : Calculate the modeled pseudorange between the satellite specified in
% sd and the receiver located at rRxEcef at the given time of signal
% reception.
%
% INPUTS
%
% gpsWeekRx -- GPS week of signal reception event in receiver time.
%
% gpsSecRx --- GPS seconds of week of signal reception event in receiver
% time.
%
% cdtRx ------ Receiver clock error scaled by the speed of light: cdtRx =
% c*dtRx. True time t is related to receiver time tRx by t = tRx - dtRx.
%
% rRxEcef ---- 3-by-1 ECEF coordinates of the receiver antenna's phase center,
% in meters in receiver timme 
%
% sd --------- Ephemeris structure array for a single SV. Let ii be the
% numerical identifier (PRN identifier) for the SV whose location
% is sought. Then sd = satdata(ii) has the following fields:
% % SVID - satellite number
% % health - satellite health flag (0 = healthy; otherwise unhealthy)
% % we - week of ephemeris epoch (GPS week, unambiguous)
% % te - time of ephemeris epoch (GPS seconds of week)
% % wc - week of clock epoch (GPS week)
% % tc - time of clock epoch (GPS seconds of week)
% % e - eccentricity (unitless)
% % sqrta - sqrt of orbit semi-major axis (m^1/2)
% % omega0 - argument of perigee (rad.)
% % M0 - mean anomaly at epoch (rad.)
% % L0 - longitude of ascending node at beginning of week (rad.)
% % i0 - inclination angle at epoch (rad.)
% % dOdt - longitude rate (rad / sec.)
% % dn - mean motion difference (rad / sec.)
% % didt - inclination rate (rad / sec.)
% % Cuc - cosine correction to argument of perigee (rad.)
% % Cus - sine correction to argument of perigee (rad.)
% % Crc - cosine correction to orbital radius (m)
% % Crs - sine correction to orbital radius (m)
% % Cic - cosine correction to inclination (rad.)
% % Cis - sine correction to inclination (rad.)
% % af0 - 0th order satellite clock correction (s)
% % af1 - 1st order satellite clock correction (s / s)
% % af2 - 2nd order satellite clock correction (s / s^2)
% % TGD - group delay time for the satellite (s)
%
% ionodata -- Ionospheric data structure array with the following fields:
% % alpha0, alpha1, alpha2, alpha3 - power series expansion coefficients
% % for amplitude of ionospheric TEC
% % beta0, beta1, beta2, beta3 - power series expansion coefficients for
% % period of ionospheric plasma density cycle
%
% tiFlags --- 2-by-1 vector of flags that indicate whether to enable
% tropospheric delay correction (tiFlag(1)) or ionospheric
% delay correction (tiFlags(2)).
%
%
% OUTPUTS
%
% rhoBar ---- Modeled pseudorange between the receiver at time of the signal
% reception event (the time specified by gpsWeekRx and gpsSecRx)
% and the satellite at time of the signal transmission event, in
% meters.
%
% Hrho ------ 1 x 4 Jacobian matrix containing partials of pseudorange with
% respect to rRxEcef and c*dtRx
%
%+------------------------------------------------------------------------------+
% References: ASE 372N Lab 4 Assignment
% ASE 372N Module 6 and 8
%
% Author: Kendall Whitbeck
%+==============================================================================+

%% Find Guess of One Way Light Time "tBar_OWLT" given tR, r_uBar_tR, & del_t_uBar_tR
navConstants; % sets values for multivarious Navigation Constants
tR = gpsSecRx; % time of Rx since GPS week epoch in seconds

% Receiver (Rx) clock error in Rx time (sec)
dtBarRx = cdtRx / cLight;

% Initial Guess for: 
t_OWLT = 0.075; % One Way Light Time (seconds)
rSatEcefTx = zeros(3,1); % "Current" position of sat at Tx in ECEF frame (m)
rSvEcef_old = ones(3,1); % old position of sat at Tx in ECEF frame (m)
count = 1;
cm_test = 0.01e-8*ones(3,1); % centimeter tests (checks if diff in estimates is less than 1 cm)

while abs(rSatEcefTx - rSvEcef_old) > cm_test == ones(3,1) % meters
    
    % Update old value of sat's position vec at Tx in ECEF frame
    rSvEcef_old = rSatEcefTx; % m
    
    % Guess of sat clock error (true time of Transmission event)
    tBar = tR - dtBarRx - t_OWLT;
    
    % Get sat position in ECEF frame using GUESS of true time of Tx
    [rSatEcefTx, ~,~] = satloc(gpsWeekRx, tBar, sd);
    
    % Updating One Way Light Time
    t_OWLT = 1/cLight*norm(rRxEcef - rSatEcefTx);
    
    % Checking to make sure convergence doesn't take too long
    if count == 10
        error('More than 10 Iterations Performed')
    end%if
    
    count = count + 1;
    
end%while loop

% Update true time of Tx
tTx = tR - dtBarRx - t_OWLT; % sec

% Get sat position in ECEF frame using GUESS of true time of Tx
[rSatEcefTx, ~,~] = satloc(gpsWeekRx, tTx, sd);
rNoSpin = rSatEcefTx;

% Accounting for Earth's rotation rate
rSatEcefTx = R3(OmegaE*t_OWLT) * rSatEcefTx;
rSpin = rSatEcefTx;

% Pseudorange Distance Earth's spin rate accounts for
dXspin = norm(rSpin - rNoSpin);

%% Relativistic Effects & Group Delay Differential

% Satdata values
tc  = sd.tc; % time of clock epoch (GPS seconds of week)
af0 = sd.af0; % af0 - 0th order satellite clock correction (s)
af1 = sd.af1; % af1 - 1st order satellite clock correction (s / s)
af2 = sd.af2; % af2 - 2nd order satellite clock correction (s / s^2)
tGD = sd.TGD; % Group Delay Differential Time (sec)

% Eccentric anomaly of GPS sat calculated from eccentricity & mean anomaly
[~, ~,E] = satloc(gpsWeekRx, tR, sd);

% Residual Relativistic Correction
dtrel = -2 * sqrt(GM) * sd.sqrta / (cLight^2) * sd.e * sin(E);

% Satellite Clock Delay at true transmission time
dtSvTx = af0 + af1*(tTx-tc) + af2*(tTx-tc)^2 + dtrel - tGD;

% Magnitude of Distance from Rx to Sat
R = abs(norm(rRxEcef - rSatEcefTx)); % meters

%% Atmospheric Delays

tGPS.week = gpsWeekRx;
tGPS.seconds = tTx;

% GNSS signal delay due to Troposphere (is 0 if tiFlags(1) = 0)
[delTauNa] = ...
    getTropoDelay(rRxEcef,rSatEcefTx,tGPS,'Saastamoinen_MSP_Neill')*tiFlags(1); % sec
T = delTauNa * cLight;

% GNSS signal delay due to Ionosphere (is 0 if tiFlags(2) = 0)
fc = 1575.42e6; % L1 carrier frequency (Hz)
[delTauG] = ...
    getIonoDelay(ionodata,fc,rRxEcef,rSatEcefTx,tGPS,'broadcast')*tiFlags(2); % sec
Ip = delTauG * cLight;

%% rhoBar

% Modeled Pseudorange btwn Rx at time of signa Rx event (the time specified
% by gpsWeekRx and gpsSecRx) the sat at time of the signal Tx event (m)
rhoBar = R + cLight * (dtBarRx - dtSvTx + delTauG + delTauNa);

%% Getting Hrho: partials of pseudorange wrt rRxEcef & c*dtRx 
% Components of the rxer & sat positions (m)
xr = rRxEcef(1); xs = rSatEcefTx(1);
yr = rRxEcef(2); ys = rSatEcefTx(2);
zr = rRxEcef(3); zs = rSatEcefTx(3);

% Patials of pseudorange wrt to rRxEcef and cdtRx
dh_dxr = -(xs-xr)/R;
dh_dyr = -(ys-yr)/R;
dh_dzr = -(zs-zr)/R;
dh_dcdtRx = 1;

% Hrho
Hrho = [dh_dxr dh_dyr dh_dzr dh_dcdtRx];

end%function