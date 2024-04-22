% topsol.m
%
% Top-level script for computing a navigation solution from
% channel.mat-formatted observables
tic
clear;clc;fclose('all'); close all
%% Setup
% Directory where observables file is stored
obsPath = '../pprx/src';

% Elevation mask angle, in deg
% 35 -> 4 SVs overhead
elMaskDeg = -90; % deg

% Exclude the following svIds from participating in the solution. SVs marked
% unhealthy will be automatically added to this list.
svIdExclude = [];

% Set start and stop indices
iiStart = 1;
num_minutes = 10; % minutes
iiStop = num_minutes*60/5;
iiStop = 100 / 5 + 1; % runs for 100 seconds
iiStop = 720; % (4 hours: iiStop = 2880)

% Set approximate (BAD) receiver antenna position
[ra,va,ya] = getAntLoc('Fiction_Island');
rRAppx = ra + va*(2016.8 - ya);
[LLAFictionIsland] = ecef2lla(rRAppx');

% % Set really GOOD guess for receiver antenna position (determined from
% % previous runs starting w/ Fiction Island as the te guess)
% latAvg = 1.067787396014510; % rad
% lonAvg = -2.614243633619133; % rad
% altAvg = 89.712965105660260; % m
% [rRAppx] = lla2ecef(latAvg,lonAvg,altAvg);

% % Does NOT converge at:
% [rRAppx] = lla2ecef(0,0,15000e3);

%% Run the Solution
initial_run = true;
num_runs = 1;
for ii = 1:num_runs
    
    if initial_run == false
        [rRAppx] = lla2ecef(latAvg,lonAvg,altAvg);
    end%if
    
    % Enable/disable tropospheric and ionospheric corrections.
    tiFlags(1) = 1; tiFlags(2) = 1;
    % Name of output file for Google Earth rendering
    kmlFileOutName = 'kmlOut.kml';
    % Change to -3 if data are so recent that no matching ephemeris is found
    ephemHourOffset = 0;
    % Ephemeris/iono data and elevation masking refresh interval
    subSolutionBlockIntervalSec = 50;
    % Observation epoch stride length
    epochStride = 1;
    
    %% Load and prepare data
    navConstants; load('channel.mat'); channelMat = channel;
   
    % ensures channelMat is in nDataPts-by-14 format
    [~, nChannel] = size(channelMat);
    if nChannel ~= 14
        channelMat = channelMat';
    end%if
    % Prepare Time History
    [tRVec,obsValidMat,svIdVec,prMat,fDMat,thetaMat] = ...
        prepareTimeHistory(channelMat,epochStride);
    delii = ceil(subSolutionBlockIntervalSec/min(diff(tRVec.s)));
    
    %% Iterate on sub-solution blocks
    % Initial Conditions
    solutionMat = []; tRVecSolution.w = []; tTrueVecSolution.w = [];
    tRVecSolution.s = []; tTrueVecSolution.s = [];
    residualVec = []; badSvIdVec = []; thetaNominalMat = [];
    P_covarMat{1} = []; P_covarD = []; covarCount = 0;
    Ip = []; T = []; dtrel = []; dXspin = [];
    
    iiStop = min(iiStop,length(tRVec.s));
    nii = length(iiStart:delii:iiStop);
    fprintf('Processing ... \n');
    for(mm = 1:nii) fprintf('='); end; fprintf('\n'); %#ok<NO4LP,SEPEX>
    for iiA = iiStart:delii:iiStop
        % Refresh ephemeris/iono data and masking
        iiB = min(iiStop,iiA + delii - 1); iiM = round(mean([iiA,iiB])); %iiMstar = 2880
        [satdata, ionodata] = retrieveNavigationData(tRVec.w(iiM),tRVec.s(iiM),...
            ephemHourOffset,'navFiles');
        [svIdAllow,] = satmap(satdata,rRAppx,elMaskDeg,tRVec.w(iiM),tRVec.s(iiM),1);
        % Exclude unhealthy PRNs
        svIdExcludeLocal = svIdExclude;
        for ii=1:length(satdata) %#ok<FXSET>
            if(satdata(ii).health)
                svIdExcludeLocal = [svIdExcludeLocal,ii];
            end%if
        end%for
        svIdAllow = setdiff(svIdAllow,svIdExcludeLocal(:)); fprintf('=');
        % Perform nav solution
        [solutionMatD,tRVecSolutionD,tTrueVecSolutionD,...
            residualVecD,badSvIdVecD,P_covar,thetaNominalMatD,IpD,TD,dtrelD,dXspinD] = ...
            performNavigationSolution(tRVec,svIdVec,obsValidMat,prMat,fDMat,...
            thetaMat,iiA,iiB,rRAppx,satdata,ionodata,...
            svIdAllow,lambdaGPSL1,tiFlags);
        % Store data
        solutionMat = [solutionMat;solutionMatD];              
        tRVecSolution.w = [tRVecSolution.w;tRVecSolutionD.w];
        tRVecSolution.s = [tRVecSolution.s;tRVecSolutionD.s];
        tTrueVecSolution.w = [tTrueVecSolution.w;tTrueVecSolutionD.w];
        tTrueVecSolution.s = [tTrueVecSolution.s;tTrueVecSolutionD.s];
        residualVec = [residualVec;residualVecD];            
        badSvIdVec = [badSvIdVec;badSvIdVecD];                
        thetaNominalMat = [thetaNominalMat;thetaNominalMatD];  
        % Task 7: Perturbations
        Ip = [Ip; IpD];                                   
        T = [T; TD];                                         
        dtrel = [dtrel; dtrelD];                             
        dXspin = [dXspin; dXspinD]; %#ok<*AGROW>
        % Covariance outputs
        covarCount = covarCount + 1;
        P_covarMat{covarCount} = P_covar;
    end%for
    
    %% Solve for average position
    rRSolAvg = mean(solutionMat(:,1:3),1)';
    cdtRSolAvg = mean(solutionMat(:,4),1)';
    [LLAavg] = ecef2lla(rRSolAvg');
    latAvg = deg2rad(LLAavg(1));
    lonAvg = deg2rad(LLAavg(2));
    altAvg = LLAavg(3);
    
    initial_run = false; % forces subsequent runs to use new inital guess
end%for

%% Solve for ENU Distribution
[Recef_to_enu]=ecef2enu(latAvg,lonAvg);
dSolution = [solutionMat(:,1)-rRSolAvg(1),solutionMat(:,2)-rRSolAvg(2),...
    solutionMat(:,3)-rRSolAvg(3)];
dSolutionENU = (Recef_to_enu*dSolution')';
dSolutionAvgFromPriorENU = Recef_to_enu*(rRAppx - rRSolAvg);
[m,n] = size(dSolutionENU);

%% Finding Dilution of Precision (DOP) Values
% Making Recef_to_enu4 4 dimensional w/ no rotation wrt time
Recef_to_enu4 = zeros(4);
Recef_to_enu4(1:3,1:3) = Recef_to_enu;
Recef_to_enu4(4,4) = 1; % don't rotate wrt time

% 3D variance in ENU frame
variance_east = sum((dSolutionENU(:,1)).^2) / (m-1); % East
variance_north = sum((dSolutionENU(:,2)).^2) / (m-1); % North
variance_up = sum((dSolutionENU(:,3)).^2) / (m-1); % Up

% Scalar variance in ENU frame
variance_enu_mag_list = zeros(m,1);
for ii = 1:m
    variance_enu_mag_list(ii,1) = norm(dSolutionENU(ii,:));
end%for
variance_enu_mag = mean(variance_enu_mag_list);

% creating cell array of H_enu (DOPS matrix) at each solution time
for ii = 1:covarCount
    P_enu = Recef_to_enu4 * (P_covarMat{ii} * Recef_to_enu4'); %Covariance Matrix ENU frame
    H_enu_list{ii} = P_enu / variance_enu_mag; % DOPS matrix ENU frame
end%for

% Concatenate same-indexed H_enu matrix components into vectors
for ii = 1:covarCount
    % First Row
    H11(ii) = H_enu_list{ii}(1,1);    H12(ii) = H_enu_list{ii}(1,2);
    H13(ii) = H_enu_list{ii}(1,3);    H14(ii) = H_enu_list{ii}(1,4);
    % Second Row
    H21(ii) = H_enu_list{ii}(2,1);    H22(ii) = H_enu_list{ii}(2,2);
    H23(ii) = H_enu_list{ii}(2,3);    H24(ii) = H_enu_list{ii}(2,4);
    % Third Row
    H31(ii) = H_enu_list{ii}(3,1);    H32(ii) = H_enu_list{ii}(3,2);
    H33(ii) = H_enu_list{ii}(3,3);    H34(ii) = H_enu_list{ii}(3,4);
    % Fourth Row
    H41(ii) = H_enu_list{ii}(4,1);    H42(ii) = H_enu_list{ii}(4,2);
    H43(ii) = H_enu_list{ii}(4,3);    H44(ii) = H_enu_list{ii}(4,4);
end%for

% Averaging H_enu components to creat normalized H_enu matrix
H_enu = [mean(H11) mean(H12) mean(H13) mean(H14);
    mean(H21) mean(H22) mean(H23) mean(H24);
    mean(H31) mean(H32) mean(H33) mean(H34);
    mean(H41) mean(H42) mean(H43) mean(H44)];

% DOP values from H_enu
VDOPavg = sqrt(H_enu(3,3));
VDOP = sqrt(H33(end));
HDOPavg = sqrt(H_enu(1,1) + H_enu(2,2)); % averaged
HDOP = sqrt(H11(end) + H22(end)); % for a single epoch
PDOP = sqrt(H_enu(1,1) + H_enu(2,2)+ H_enu(3,3));
GDOP = sqrt(trace(H_enu));

% Getting best sigma_rho value to fit data (makes HDOP & VDOP fit
% respective variance calculated by difference-from-mean-soultion
sig_rho_vec(1) = variance_east / HDOP; % compared to up direction
sig_rho_vec(2) = variance_north / HDOP; % compared to up direction
sig_rho_vec(3) = variance_up / VDOP; % compared to up direction
sig_rho_best = mean(sig_rho_vec);

%% Perturbations
Ip_avg = mean(abs(Ip)); % Ionospheric (meters)
T_avg = mean(abs(T)); % Tropospheric (meters)
cxdtrel_avg = mean(abs(dtrel))*cLight; % Relativistic (meters)
dXspin_avg = mean(abs(dXspin)); % Earth spin rate(meters)

% Percent of perturbation due to Earth spin rate
percent_dXspin = dXspin_avg/sum([Ip_avg,T_avg,cxdtrel_avg,dXspin_avg])*100;

%% L1 vs L2

% SVID 19
% creates col vescs of carrier phase values at L1 & L2 for SVID #19
phi1VecSV19 = channelMat((channelMat(:, 13) == 0) + (channelMat(:, 14) == 19)==2, 7);
phi2VecSV19 = channelMat((channelMat(:, 13) == 2) + (channelMat(:, 14) == 19)==2, 7);
% creates col vescs of pseudorange values at L1 & L2 for SVID #6
rho1VecSV19 = channelMat((channelMat(:, 13) == 0) + (channelMat(:, 14) == 19)==2, 8);
rho2VecSV19 = channelMat((channelMat(:, 13) == 2) + (channelMat(:, 14) == 19)==2, 8);
% Ionospheric Delay on L1 for carrier phase & pseudorange
I_phi1VecSV19 = (fL2GPS^2/(fL1GPS^2-fL2GPS^2)*(lambdaGPSL1*phi1VecSV19 - lambdaGPSL2*phi2VecSV19));
I_rho1VecSV19 = fL2GPS^2/(fL1GPS^2-fL2GPS^2)*(rho2VecSV19 - rho1VecSV19);
% Bias (mean of the difference)
I_diff1VecSV19 = I_rho1VecSV19 - I_phi1VecSV19; % meters
I_bias1SV19 = mean(I_diff1VecSV19);  % meters
I_phi1UnbiasedSV19 = I_phi1VecSV19 + I_bias1SV19; % meters
% Root Mean Square of pseudorange & corrected carrier phase
RMSrho1SV19 = rms(I_rho1VecSV19); % pseudorange (m)
RMSphi1SV19 = rms(I_phi1UnbiasedSV19); % carrier phase corrected (m)
PercentRMSdiff1SV19 = (RMSrho1SV19 - RMSphi1SV19)/RMSphi1SV19*100; % percent

% SVID 25
% creates col vescs of carrier phase values at L1 & L2 for SVID #6
phi1VecSV3 = channelMat((channelMat(:, 13) == 0) + (channelMat(:, 14) == 3)==2, 7);
phi2VecSV3 = channelMat((channelMat(:, 13) == 2) + (channelMat(:, 14) == 3)==2, 7);
% creates col vescs of pseudorange values at L1 & L2 for SVID #6
rho1VecSV3 = channelMat((channelMat(:, 13) == 0) + (channelMat(:, 14) == 3)==2, 8);
rho2VecSV3 = channelMat((channelMat(:, 13) == 2) + (channelMat(:, 14) == 3)==2, 8);
% Ionospheric Delay on L1 for carrier phase & pseudorange
I_phi1VecSV3 = (fL2GPS^2/(fL1GPS^2-fL2GPS^2)*(lambdaGPSL1*phi1VecSV3 - lambdaGPSL2*phi2VecSV3));
I_rho1VecSV3 = fL2GPS^2/(fL1GPS^2-fL2GPS^2)*(rho2VecSV3 - rho1VecSV3);
% Bias (mean of the difference)
I_diff1VecSV3 = I_rho1VecSV3 - I_phi1VecSV3; % meters
I_bias1SV3 = mean(I_diff1VecSV3);  % meters
I_phi1UnbiasedSV3 = I_phi1VecSV3 + I_bias1SV3; % meters
% Root Mean Square of pseudorange & corrected carrier phase
RMSrho1SV3 = rms(I_rho1VecSV3); % pseudorange (m)
RMSphi1SV3 = rms(I_phi1UnbiasedSV3); % carrier phase corrected (m)
PercentRMSdiff1SV3 = (RMSrho1SV3 - RMSphi1SV3)/RMSphi1SV3*100; % percent

% number of data points recorded for iono delay values
nData = length(I_phi1VecSV19);

%% Plot results
iiVec = iiStart:(iiStart + length(tTrueVecSolution.s) - 1);

figure(5);clf;
plot(dSolutionENU(:,1),dSolutionENU(:,2), 'x','markersize', 6);
hold on;
% plot(dSolutionAvgFromPriorENU(1),dSolutionAvgFromPriorENU(2), 'r.', ...
%     'markersize', 30);
xlabel('East (meters)'); ylabel('North (meters)');
title('Horizontal spread of navigation solution points with origin at the mean');
grid on; axis equal

figure(6);clf;
plot(iiVec,dSolutionENU(:,3));
hold on;
% plot(iiVec,ones(length(iiVec),1)*dSolutionAvgFromPriorENU(3), 'r', ...
%     'linewidth', 2);
xlabel('Solution Number'); ylabel('displacement (meters)');
title('Vertical displacement from mean');

figure(7);clf;
plot(iiVec,residualVec);
xlabel('Solution Number'); ylabel('max \Delta z (meters)');
title('Worst case residual at each solution epoch');

figure(8);clf;
plot(iiVec,badSvIdVec,'.')%,'markersize', 30);
title('SV with largest residual at each epoch');
xlabel('Solution Number');
grid on;

%% New Plots

figure(9); clf;
c = categorical({'Ionospheric Delay', 'Neutral Atmosperic Delay', 'Relativistic', 'Earth Spin Rate'});
b = [Ip_avg, T_avg, cxdtrel_avg, dXspin_avg];
bar(c,b)
title('Average Perturbation Distance on Pseudorange')
ylabel('Perturbation (m)')
grid on

figure(11);clf; hold on;
plot((1:nData)*5,I_phi1UnbiasedSV19);
plot((1:nData)*5,I_rho1VecSV19);
xlim([0 nData*5])
title('Ionospheric Delay from Pseudorange & Carrier Phase for SVID #19');
xlabel('GPS Seconds');
ylabel('Ionospheric Delay (m)');
legend('Carrier Phase (Unbiased)','Pseudorange')
grid on; hold off;

figure(12);clf; hold on;
plot((1:nData)*5,I_phi1UnbiasedSV3);
plot((1:nData)*5,I_rho1VecSV3);
xlim([0 nData*5])
title('Ionospheric Delay from Pseudorange & Carrier Phase for SVID #3');
xlabel('GPS Seconds');
ylabel('Ionospheric Delay (m)');
legend('Carrier Phase (Unbiased)','Pseudorange')
grid on; hold off;

%% Print solution
fprintf('\n\nMean navigation solution:\n\n');
% ECEF position Solutions
fprintf(' x-Position ECEF (meters): %+f\n', rRSolAvg(1));
fprintf(' y-Position ECEF (meters): %+f\n', rRSolAvg(2));
fprintf(' z-Position ECEF (meters): %+f\n', rRSolAvg(3));
fprintf(' c*dtR           (meters): %+f\n\n',cdtRSolAvg);
% Topocentric solution
fprintf(' Latitude (deg)   : %+f\n', latAvg*180/pi);
fprintf(' Longitude (deg)  : %+f\n', lonAvg*180/pi);
fprintf(' Altitude (meters): %+f\n\n', altAvg);
% DOPs
fprintf(' HDOP: %+f\n', HDOP);
fprintf(' VDOP: %+f\n', VDOP);
% Comparing DOPS & Variance
fprintf('  Let \x03c3_\x03c1 = 2\n')
fprintf(' HDOP*\x03c3_\x03c1 = %+f vs. \x03c3_E^2 = %+f\n', HDOP*2,variance_east);
fprintf('                          \x03c3_N^2 = %+f\n',variance_north);
fprintf(' VDOP*\x03c3_\x03c1 = %+f vs. \x03c3_U^2 = %+f\n', VDOP*2,variance_up);
fprintf(' Best Fitting Std. Dev.:  \x03c3_\x03c1   = %+f\n\n', sig_rho_best);
% L1 vs L2 for Pseudorange & Carier Phase SV3
fprintf(' SV19 RMS for I_\x03c1 (m): %+f\n', RMSrho1SV19);
fprintf(' SV19 RMS for I_\x03A6 (m): %+f\n', RMSphi1SV19);
fprintf(' Percent Difference  : %+f%%\n\n', PercentRMSdiff1SV19)
% L1 vs L2 for Pseudorange & Carier Phase SV3
fprintf(' SV3 RMS for I_\x03c1 (m): %+f\n', RMSrho1SV3);
fprintf(' SV3 RMS for I_\x03A6 (m): %+f\n', RMSphi1SV3);
fprintf(' Percent Difference  : %+f%%\n\n', PercentRMSdiff1SV3)

% Other Default stuff
fprintf('Mean epoch spacing  (seconds): %f\n', ...
    mean(diff((tTrueVecSolution.w - tTrueVecSolution.w(1))*sec_in_week + ...
    tTrueVecSolution.s)));
fprintf('Total time interval (seconds): %f\n', ...
    (tTrueVecSolution.w(end) - tTrueVecSolution.w(1))*sec_in_week + ...
    tTrueVecSolution.s(end) - tTrueVecSolution.s(1));
fprintf('Solution distance from initial guess (meters): %f\n', ...
    norm(dSolutionAvgFromPriorENU));
fprintf('Solution horizontal distance from initial guess (meters): %f\n', ...
    norm(dSolutionAvgFromPriorENU(1:2)));
fprintf('Solution vertical distance from initial guess   (meters): %f\n', ...
    abs(dSolutionAvgFromPriorENU(3)));

%% Comparing first 6 iterations of a solution epoch
xMat = [-4956383.91578323,-2726959.48532797,-2663918.84217939,-2663848.94836890,-2663848.94893216,-2663848.94893216;-1982306.40652170,-1570602.95031319,-1551349.86436850,-1551322.61280169,-1551322.61322245,-1551322.61322244;6904346.63771916,5608884.61253962,5565161.01478375,5565102.50411840,5565102.50527626,5565102.50527624;2986387.03002877,132698.241085657,128.532245651733,0.0126398086660287,0.0136149293944008,0.0136149148665373];

figure
plot3(xMat(1,:),xMat(2,:),xMat(3,:))
title('First Six Iterations of Position Components of State Vector')
xlabel('x-Position (m)')
ylabel('y-Position (m)')
zlabel('z-Position (m)')
grid on

fprintf('\n\nFirst Iteration\n')
fprintf(' x-Position ECEF (meters): %+f\n', xMat(1,1));
fprintf(' y-Position ECEF (meters): %+f\n', xMat(2,1));
fprintf(' z-Position ECEF (meters): %+f\n', xMat(3,1));
fprintf('Second Iteration\n')
fprintf(' x-Position ECEF (meters): %+f\n', xMat(1,2));
fprintf(' y-Position ECEF (meters): %+f\n', xMat(2,2));
fprintf(' z-Position ECEF (meters): %+f\n', xMat(3,2));
fprintf('Third Iteration\n')
fprintf(' x-Position ECEF (meters): %+f\n', xMat(1,3));
fprintf(' y-Position ECEF (meters): %+f\n', xMat(2,3));
fprintf(' z-Position ECEF (meters): %+f\n', xMat(3,3));
fprintf('Fourth Iteration\n')
fprintf(' x-Position ECEF (meters): %+f\n', xMat(1,4));
fprintf(' y-Position ECEF (meters): %+f\n', xMat(2,4));
fprintf(' z-Position ECEF (meters): %+f\n', xMat(3,4));
fprintf('Fifth Iteration\n')
fprintf(' x-Position ECEF (meters): %+f\n', xMat(1,5));
fprintf(' y-Position ECEF (meters): %+f\n', xMat(2,5));
fprintf(' z-Position ECEF (meters): %+f\n', xMat(3,5));
fprintf('Sixth Iteration\n')
fprintf(' x-Position ECEF (meters): %+f\n', xMat(1,6));
fprintf(' y-Position ECEF (meters): %+f\n', xMat(2,6));
fprintf(' z-Position ECEF (meters): %+f\n', xMat(3,6));

%% Generate Google Earth KML file
% pMat = solutionMat(:,1:3);
% tVec.week = tTrueVecSolution.w;
% tVec.sec = tTrueVecSolution.s;
% genKmlFile(pMat,tVec,,kmlFileOutName);
% st = ['! cp ' kmlFileOutName ' /home/todd '];
% eval(st);

% %----- Create navsol.mat for use in CDGNSS processing; save in obsPath
% columnDefinitions;
% wkInvalid = 9999;
% nS = length(tRVecSolution.s);
% navsol = zeros(nS,12);
% navsol(:,ns_ortWkCol) = tRVecSolution.w;
% navsol(:,ns_ortWholeSecCol) = floor(tRVecSolution.s);
% navsol(:,ns_ortFracSecCol) = tRVecSolution.s - floor(tRVecSolution.s);
% navsol(:,4:7) = solutionMat;
% navsol(:,12) = 2*ones(nS,1);
% navsol = navsol';
% save([obsPath '/navsol.mat'], 'navsol');


disp(' ')
toc