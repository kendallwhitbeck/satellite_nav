%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASE 372N: Satellite-Based Navigation
% HOMEWORK 1 MAIN FILE
close all; format long g; clc; clear all;
set(0, 'DefaultAxesFontSize',18, 'DefaultLineLineWidth',3,...
    'DefaultLineMarkerSize',16)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 4: First Session, Mustang Statue

% Reading the NMEA file for the first site, first session
[sod_M1_v,lat_M1_v,lon_M1_v,alt_M1_v,sats_M1_v] = readNMEAposV3('mustang1.nmea');

% Calculating the mean of each value from the NMEA file
sod_M1 = mean(sod_M1_v);
lat_M1 = mean(lat_M1_v);
lon_M1 = mean(lon_M1_v);
alt_M1 = mean(alt_M1_v);
sats_M1 = mean(sats_M1_v);

% Mean ECEF position vector for Mustang1 in meters
[M1_ECEF] = lla2ecef(lat_M1*pi/180,lon_M1*pi/180,alt_M1)

%% Problem 4: Second Session, Mustang Statue
% Reading the NMEA file for the first site, first session
[sod_M2_v,lat_M2_v,lon_M2_v,alt_M2_v,sats_M2_v] = readNMEAposV3('mustang2.nmea');
% Calculating the mean of each value from the NMEA file
sod_M2 = mean(sod_M2_v);
lat_M2 = mean(lat_M2_v);
lon_M2 = mean(lon_M2_v);
alt_M2 = mean(alt_M2_v);
sats_M2 = mean(sats_M2_v);

% Mean ECEF position vector for Mustang2 in meters
[M2_ECEF] = lla2ecef(lat_M2*pi/180,lon_M2*pi/180,alt_M2)

%% Problem 4: First Session, Dinosaur Footprint

% Reading the NMEA file for the first site, first session into vectors
[sod_D1_v,lat_D1_v,lon_D1_v,alt_D1_v,sats_D1_v] = readNMEAposV3('dino1.nmea');

% Calculating the mean of each value from the NMEA file
sod_D1 = mean(sod_D1_v);
lat_D1 = mean(lat_D1_v);
lon_D1 = mean(lon_D1_v);
alt_D1 = mean(alt_D1_v);
sats_D1 = mean(sats_D1_v);

% Mean ECEF position vector for Dino1 in meters
[D1_ECEF] = lla2ecef(lat_D1*pi/180,lon_D1*pi/180,alt_D1)

%% Problem 4: Second Session, Dinosaur Footprint

% Reading the NMEA file for the first site, first session
[sod_D2_v,lat_D2_v,lon_D2_v,alt_D2_v,sats_D2_v] = readNMEAposV3('dino2.nmea');

% Calculating the mean of each value from the NMEA file
sod_D2 = mean(sod_D2_v);
lat_D2 = mean(lat_D2_v);
lon_D2 = mean(lon_D2_v);
alt_D2 = mean(alt_D2_v);
sats_D2 = mean(sats_D2_v);

% Mean ECEF position vector for Dino2 in meters
[D2_ECEF] = lla2ecef(lat_D2*pi/180,lon_D2*pi/180,alt_D2)

%% Problem 4: First Session, ETC Bridge

% Reading the NMEA file for the first site, first session
[sod_E1_v,lat_E1_v,lon_E1_v,alt_E1_v,sats_E1_v] = readNMEAposV3('etc1.nmea');

% Calculating the mean of each value from the NMEA file
sod_E1 = mean(sod_E1_v);
lat_E1 = mean(lat_E1_v);
lon_E1 = mean(lon_E1_v);
alt_E1 = mean(alt_E1_v);
sats_E1 = mean(sats_E1_v);

% Mean ECEF position vector for ETC1 in meters
[E1_ECEF] = lla2ecef(lat_E1*pi/180,lon_E1*pi/180,alt_E1)

%% Problem 4: Second Session, ETC Bridge

% Reading the NMEA file for the first site, first session
[sod_E2_v,lat_E2_v,lon_E2_v,alt_E2_v,sats_E2_v] = readNMEAposV3('etc2.nmea');

% Calculating the mean of each value from the NMEA file
sod_E2 = mean(sod_E2_v); % sec
lat_E2 = mean(lat_E2_v); % deg
lon_E2 = mean(lon_E2_v); % deg
alt_E2 = mean(alt_E2_v); % meters
sats_E2 = mean(sats_E2_v); % number of sats

% Mean ECEF position vector for ETC2 in meters
[E2_ECEF] = lla2ecef(lat_E2*pi/180,lon_E2*pi/180,alt_E2)

%% Problem 5:
% Averaged distance btwn the mustang & dino markers in meters in ECEF frame
dist_M1_D1 = norm(M1_ECEF - D1_ECEF); % from Mustang1 to DinoFoot1
dist_M1_D2 = norm(M1_ECEF - D2_ECEF); % from Mustang1 to DinoFoot2
dist_M2_D2 = norm(M2_ECEF - D2_ECEF); % from Mustang2 to DinoFoot2
dist_M_D = (dist_M1_D1 + dist_M1_D2 + dist_M2_D2)/3 % averaged estimate

% Averaged distance btwn the mustang & ETC markers in meters in ECEF frame
dist_M1_E1 = norm(M1_ECEF - E1_ECEF); % from Mustang1 to ETC1
dist_M1_E2 = norm(M1_ECEF - E2_ECEF); % from Mustang1 to ETC2
dist_M2_E2 = norm(M2_ECEF - E2_ECEF); % from Mustang2 to ETC2
dist_M_E = (dist_M1_E1 + dist_M1_E2 + dist_M2_E2)/3 % averaged estimate

% Averaged distance btwn the dino & ETC markers in meters in ECEF frame
dist_D1_E1 = norm(D1_ECEF - E1_ECEF); % from Mustang1 to ETC1
dist_D1_E2 = norm(D1_ECEF - E2_ECEF); % from Mustang1 to ETC2
dist_D2_E2 = norm(D2_ECEF - E2_ECEF); % from Mustang2 to ETC2
dist_D_E = (dist_D1_E1 + dist_D1_E2 + dist_D2_E2)/3 % averaged estimate

%% Problem 6, Part (a): Mustang Statue delta-lat vs. delta-long
% Calculating delta-longitude & delta-latitude for Mustang 1
dlon_M1 = lon_M1_v - lon_M1; % deg
dlat_M1 = lat_M1_v - lat_M1; % deg

% Calculating delta-longitude & delta-latitude for Mustang 2
dlon_M2 = lon_M2_v - lon_M2; % deg
dlat_M2 = lat_M2_v - lat_M2; % deg

% Scatter Plots of delta-lat vs. delta-lon for both Mustang data samples
figure
hold on
grid on
axis square
scatter(dlon_M1, dlat_M1, 'filled','MarkerEdgeColor', 'k')
scatter(dlon_M2, dlat_M2, 'filled', 'Marker', 'd','MarkerEdgeColor', 'k')
title('6a. Change in Longitude & Latitude for Mustang Statue at Two Times', 'FontSize', 17)
xlabel('delta-Longitude (degrees)')
ylabel('delta-Latitude (degrees)')
legend('Mustang 1', 'Mustang 2')
hold off

%% Problem 6, Part (b): Mustang Statue delta-east vs. delta-north

% Rotation matrix from ECEF to ENU using mean lat & mean long for:
% Mustang1
[R_M1] = ecef2enu(lat_M1,lat_M1);
% Mustang2
[R_M2] = ecef2enu(lat_M2,lat_M2);

% Pre-allocating For-loop variables for speed
n = length(lat_M1_v);
M1_ECEF_mat = zeros(3,n);
M1_ECEF_diff_mat = zeros(3,n);
M1_ENU_diff_mat = zeros(3,n);
dEast_M1 = zeros(1,n);
dNorth_M1 = zeros(1,n);
M1_ENU_pos_mat = zeros(3,n);
M2_ECEF_mat = zeros(3,n);
M2_ECEF_diff_mat = zeros(3,n);
M2_ENU_diff_mat = zeros(3,n);
dEast_M2 = zeros(1,n);
dNorth_M2 = zeros(1,n);
M2_ENU_pos_mat = zeros(3,n);
% For-loop computing an ECEF pos vec at each point in time for M1 & M2
for i = 1:n
    
    % Mustang 1
    
    % Matrix of ECEF pos vecs (meters) for M1 where each row is a position 
    % component (x, y, or z) & each col is a new point in time
    M1_ECEF_mat(:,i) = lla2ecef(lat_M1_v(i)*pi/180,lon_M1_v(i)*pi/180,alt_M1_v(i));

    % Matrix of ECEF position difference vecs (meters) for M1 where each
    % row is a pos component (x, y, or z) & each col is a new point in time
    M1_ECEF_diff_mat(:,i) = M1_ECEF_mat(:,i) - M1_ECEF;
    
    % Matrix of ENU position difference vecs (meters) for M1 where each
    % row is a pos component (E, N, or U) & each col is a new point in time
    M1_ENU_diff_mat(:,i) = R_M1 * M1_ECEF_diff_mat(:,i);
    
    % Vector of East-component differences for Mustang 1
    dEast_M1(i) = M1_ENU_diff_mat(1,i); % meters
    
    % Vector of North-component differences for Mustang 1
    dNorth_M1(i) = M1_ENU_diff_mat(2,i); % meters

    % Matrix of ENU position vecs (meters) for M1 where each row is a pos
    % component (E, N, or U) & each col is a new point in time
    M1_ENU_pos_mat(:,i) = R_M1 * M1_ECEF_mat(:,i);
    
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
    
    % Mustang 2
    
    % Matrix of ECEF pos vecs (meters) for M2 where each row is a position 
    % component (x, y, or z) & each col is a new point in time
    M2_ECEF_mat(:,i) = lla2ecef(lat_M2_v(i)*pi/180,lon_M2_v(i)*pi/180,alt_M2_v(i));
    
    % Matrix of ECEF position difference vecs (meters) for M2 where each
    % row is a pos component (x, y, or z) & each col is a new point in time
    M2_ECEF_diff_mat(:,i) = M2_ECEF_mat(:,i) - M2_ECEF;
        
    % Matrix of ENU position difference vecs (meters) for M2 where each
    % row is a pos component (E, N, or U) & each col is a new point in time
    M2_ENU_diff_mat(:,i) = R_M2 * M2_ECEF_diff_mat(:,i);

    % Vector of East-component differences for Mustang 2
    dEast_M2(i) = M2_ENU_diff_mat(1,i); % meters
    
    % Vector of North-component differences for Mustang 2
    dNorth_M2(i) = M2_ENU_diff_mat(2,i); % meters

    % Matrix of ENU position vecs (meters) for M1 where each row is a pos
    % component (E, N, or U) & each col is a new point in time
    M2_ENU_pos_mat(:,i) = R_M2 * M2_ECEF_mat(:,i);

end

% Scatter Plots of delta-north vs. delta-east for both Mustang data samples
figure
hold on
grid on
axis square
scatter(dEast_M1, dNorth_M1, 'filled','MarkerEdgeColor', 'k')
scatter(dEast_M2, dNorth_M2, 'filled', 'Marker', 'd','MarkerEdgeColor', 'k')
title('6b. Delta-North vs. Delta-East for Mustang Statue at Two Different Times', 'FontSize', 18)
xlabel('delta-East (meters)')
ylabel('delta-North (meters)')
legend('Mustang 1', 'Mustang 2')
hold off

%% Problem 6, Part (c): Error for Mustang 1 & 2

% Defining the collection & averages of position coordinates for East,
% North, & Up for:
% Mustang 1 
M1_E_pos = M1_ENU_pos_mat(1,:); % meters
M1_N_pos = M1_ENU_pos_mat(2,:); % meters
M1_U_pos = M1_ENU_pos_mat(3,:); % meters
% Mustang 2 
M2_E_pos = M2_ENU_pos_mat(1,:); % meters
M2_N_pos = M2_ENU_pos_mat(2,:); % meters
M2_U_pos = M2_ENU_pos_mat(3,:); % meters

% Computing the averaged position coordinates for East, North, & Up for:
% Mustang 1 
M1_E_avg = mean(M1_E_pos); % meters
M1_N_avg = mean(M1_N_pos); % meters
M1_U_avg = mean(M1_U_pos); % meters
% Mustang 2 
M2_E_avg = mean(M2_E_pos); % meters
M2_N_avg = mean(M2_N_pos); % meters
M2_U_avg = mean(M2_U_pos); % meters

% Calculate square of standard deviation of ENU components for:
% Mustang 1
std_E_sq_M1 = 1/(length(M1_E_pos)-1) * sum((M1_E_pos - M1_E_avg).^2);
std_N_sq_M1 = 1/(length(M1_N_pos)-1) * sum((M1_N_pos - M1_N_avg).^2);
std_U_sq_M1 = 1/(length(M1_U_pos)-1) * sum((M1_U_pos - M1_U_avg).^2);
% Mustang 2
std_E_sq_M2 = 1/(length(M2_E_pos)-1) * sum((M2_E_pos - M2_E_avg).^2);
std_N_sq_M2 = 1/(length(M2_N_pos)-1) * sum((M2_N_pos - M2_N_avg).^2);
std_U_sq_M2 = 1/(length(M2_U_pos)-1) * sum((M2_U_pos - M2_U_avg).^2);

% Computing 2-D rms error to get Circular Error Probable (CEP) for
% Mustang 1 & 2
rms_2D_error_M1 = sqrt(std_E_sq_M1 + std_N_sq_M1); % meters
rms_2D_error_M2 = sqrt(std_E_sq_M2 + std_N_sq_M2); % meters
CEP_M1 = rms_2D_error_M1 / 1.2 % meters
CEP_M2 = rms_2D_error_M2 / 1.2 % meters

% Computing 3-D rms error to get Spherical Error Probable (SEP) for
% Mustang 1 & 2
rms_3D_error_M1 = sqrt(std_E_sq_M1 + std_N_sq_M1 + std_U_sq_M1); % meters
rms_3D_error_M2 = sqrt(std_E_sq_M2 + std_N_sq_M2 + std_U_sq_M2); % meters
SEP_M1 = rms_3D_error_M1 / 1.3 % meters
SEP_M2 = rms_3D_error_M2 / 1.3 % meters

%% Problem 7, Part (b)
% Array of Altitude Residuals compared to mean altitude for Mustang 1 & 2
alt_res_M1_v = alt_M1_v - alt_M1;
alt_res_M2_v = alt_M2_v - alt_M2;

% Computing the averaged Altitude Residuals for Mustang 1 & 2
alt_res_M1_avg = mean(abs(alt_res_M1_v)) % meters
alt_res_M2_avg = mean(abs(alt_res_M2_v)) % meters

% Computing the standard deviation of altitude residuals for Mustang 1 & 2
std_alt_res_M1 = sqrt(1/(length(alt_res_M1_v)-1) * sum((alt_res_M1_v - alt_res_M1_avg).^2))
std_alt_res_M2 = sqrt(1/(length(alt_res_M2_v)-1) * sum((alt_res_M2_v - alt_res_M2_avg).^2))

%% Problem 7, Part (c): Precision
% Geometric Dilution of Precision (GDOP) for Mustang 1 & 2
PDOP_M1 = 1.7; TDOP_M1 = 1;
PDOP_M2 = 1.7; TDOP_M2 = 1;
GDOP = sqrt(PDOP_M1^2 + TDOP_M1^2) % for both Mustang 1 & 2

%% +----------------------------------------------------------------------+
% References:
% Moriba K. Jah, Ph.D. - ASE 372N: Module 2
% Misra & Enge - Global Positioning System: Signals, Measurements, and Performance, 2nd Ed.
% National Academies Press - The Global Positioning System: A Shared National Asset
%
% Author: Kendall Whitbeck
%+========================================================================+