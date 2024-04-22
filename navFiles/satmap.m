function [svIdVec, plotMat] = ...
    satmap(satdata,rRxEcef,elMaskDeg,gpsWeek,gpsSecVec,plotFlag)
% satmap : Generate plotting data for the SVs above a particular receiver
% position over a span of time.
%
%
% INPUTS
%
% satdata ------ Ephemeris structure array; see getephem header for details.
%
% rRxEcef ------ 3-by-1 receiver ECEF position, in meters.
%
% elMaskDeg ---- Elevation mask angle, in degrees.
%
% gpsWeek ------ GPS week number corresponding to the first element in
% gpsSecVec.
%
% gpsSecVec ---- Nt-by-1 vector of GPS seconds of week over which an SV trace
% is desired. Entries exceeding 7*86400 indicate that the time
% interval spans a GPS week boundary.
%
% plotFlag ----- Indicates whether (1) or not (0) to generate a sky plot of
% SVs.
%
%
% OUTPUTS
%
% svIdVec ------ Nsv-by-1 vector of unique SV identification numbers for the
% SVs that were above the elevation mask angle at any time
% during the interval.
%
% plotMat ------ Nt*Nsv-by-4 matrix of data corresponding to SVs in svIdVec,
% arranged as
%
% [svID gpsSec elRad azRad;
% svID gpsSec elRad azRad;
% svID gpsSec elRad azRad]
%
% where svId is the SV identification number, gpsSec is GPS
% time in seconds of week, and elRad and azRad are the SV
% elevation and azimuth angles, in radians. plotMat is
% composed of Nt batches of Nsv rows each, where the ith batch
% corresponds to time gpsSecVec(i). Each batch has data for
% all Nsv svIDs in svIdVec. Note that some elRad values may be
% below the elevation mask angle: this simply means that the
% corresponding svID rose or set during the data interval.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author: Kendall Whitbeck
%+==============================================================================+

% Converting elMask from degrees to radians
elMaskRad = elMaskDeg*pi/180; % rads

% Shortening variable name for readability
sd = satdata;

% Preallocating varibales for speed
ii_len = length(satdata);
tt_len = length(gpsSecVec);
rSvEcef = zeros(3,ii_len);
% vSvEcef = zeros(3,ii_len);
elRad = zeros(ii_len,tt_len);
azRad = zeros(ii_len,tt_len);
plotMat = zeros(1,4);
svIdVec = zeros(1);
hh = 1;

% for each iteration in time
for tt = 1:tt_len
    
    % for each SV
    for ii = 1:ii_len
        
        % Checks to see if current SV has an ephemeris (i.e. a healthy sat
        % was available at epoch measurement) & skips it if not
        if isempty(sd(ii).SVID) == true
            continue
        end%if
        
        % Calculating ECEF position & velocity vectors for satellite
        [rSvEcef(:,ii), ~ ] = satloc(gpsWeek, gpsSecVec(tt), sd(ii));

        % Calculating Elevation & Azimuth of satellite wrt receiver (rads)s
        [elRad(ii), azRad(ii)] = satelaz(rSvEcef(:,ii),rRxEcef);
        
        % Creates vector of SVID numbers for SVs above Elevation mask.
        % Duplicate SVIDs are removed from svIdVec after double for-loop
        if elRad(ii) > elMaskRad
            svIdVec(hh) = sd(ii).SVID;            
            hh = hh + 1; % increase # rows in svIdVec for next iteration
        end%if
        
        % Nt*NSv-by-4 matrix where:
        % % new row = new SV, in batches of NSv where Nt is num of time batches
        % % new column = different property
        plotMatNew(ii,:) = [sd(ii).SVID gpsSecVec(tt) elRad(ii) azRad(ii)];

    end%for
    
    % Concatenate previous plotMat time batch w/ new time batch
    plotMat = [plotMat; plotMatNew];
    
end%for

% Removing duplicate SVIDs & empty ephemeris values from svIdVec
svIdVec = unique(nonzeros(svIdVec));

% Removing empty ephemeris values from plotmat
plotMatCol1 = nonzeros(plotMat(:,1,:));
plotMatCol2 = nonzeros(plotMat(:,2,:));
plotMatCol3 = nonzeros(plotMat(:,3,:));
plotMatCol4 = nonzeros(plotMat(:,4,:));
plotMat = [plotMatCol1, plotMatCol2, plotMatCol3, plotMatCol4];

% % Checks if a plot of satellite positions/trajectories is requested
if plotFlag == 1
    jj = 1;
    for ii = 1 : length(plotMat(:,1))
        if ismember(plotMat(ii,1),svIdVec)
            plotMatMask(jj,:) = plotMat(ii,:);
            jj = jj + 1;
        end%if
    end%for
    plotsat(plotMatMask(:,:),gpsWeek,elMaskRad);
end %if
end %function
%%