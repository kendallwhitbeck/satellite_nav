function [solutionMat,tRVecSolution,tTrueVecSolution,...
    residualVec,badSvIdVec,P_covar,thetaNominalMat,Ip_vec,T_vec,dtrel_vec,dXspin] = ...
    performNavigationSolution(tRVec,svIdVec,obsValidMat,prMat,fDMat,thetaMat,...
    iiStart,iiStop,rRxAppx,satdata,ionodata,svIdAllow,lambda,tiFlags)
% performNavigationSolution : Perform the navigation solution using
% nonlinear least squares techniques.
%
% For ASE 372N Lab 4, don't need:
%
% INPUTS
%
% tRVec --------------- Structure of two Nt-by-1 vectors that jointly
% define the unique receiver time instants corresponding to all
% received observables; tRVec.w contains the GPS week
% and tRVec.s contains the GPS seconds of week.
%
% svIdVec ------------- Nsv-by-1 vector of unique SVIDs corresponding to
% SVs that were tracked at some point during the
% measurement interval spanned by tRVec.
%
% obsValidMat --------- Nt-by-Nsv matrix whose (ii,jj)th element indicates
% whether (1) or not (0) observables were valid for SVID
% svIdVec(jj) at time tRVec.s(ii).
%
% prMat --------------- Nt-by-Nsv matrix whose (ii,jj)th element is the
% measured pseudorange value for SVID svIdVec(jj) at
% time tRVec.s(ii), in meters.
%
% fDMat --------------- Nt-by-Nsv matrix whose (ii,jj)th element is the
% measured Doppler value for SVID svIdVec(jj) at time
% tRVec.s(ii), in Hz.
%
% thetaMat ------------ Nt-by-Nsv matrix whose (ii,jj)th element is the
% measured beat carrier phase value for SVID svIdVec(jj)
% at time tRVec.s(ii), in cycles.
%
% iiStart,iiStop ------ Start and stop indices into tRVec between which a
% navigation solution will be computed. These indices
% are used to isolate a specific section of data. The
% number of solutions Ns will nominally be Ns = iiStop -
% iiStart + 1 but will be less if for some epochs a
% solution could not be computed (e.g., too few SVs).
%
% rRxAppx ------------- 3-by-1 approximate location of receiver antenna, in
% ECEF meters.
%
% satdata, ionodata --- SV and ionospheric navigation data parameters.
% See getephem.m for details.
%
% svIdAllow ----------- Vector of SVIDs for SVs that are allowed to
% participate in the navigation solution.
%
% lambda -------------- Nominal wavelength of GNSS signal, in meters.
%
% tiFlags ------------- 2-by-1 vec of flags that indicate whether to enable
% tropospheric delay correction (tiFlag(1)) or
% ionospheric delay correction (tiFlags(2)).
%
%
% OUTPUTS
%
% solutionMat --------- Ns-by-4 matrix of navigation solutions whose rows
% are formatted as [xEcef,yEcef,zEcef,c*deltR].
%
% tRVecSolution ------- Structure of two Ns-by-1 vectors that jointly define
% the unique receiver time instants corresponding to
% navigation solutions in solutionMat; tRVecSolution.w
% contains the GPS week and tRVecSolution.s contains the
% GPS seconds of week.
%
% tTrueVecSolution ---- Structure of two Ns-by-1 vectors that jointly define
% the unique receiver time instants corresponding to
% navigation solutions in solutionMat. This vector is
% an estimate of ’true’ GPS time for each
% solution. tTrueVecSolution.w contains the GPS week and
% tTrueVecSolution.s contains the GPS seconds of week.
%
% residualVec --------- Ns-by-1 vector of maximum least-squares residuals
% magnitudes from each solution epoch, in meters.
%
% badSvIdVec ---------- Ns-by-1 vector of SVIDs corresponding to the SV with
% the worst residual at each solution epoch.
%
% thetaNominalMat ----- Ns-by-Nsv vector of nominal beat carrier phase values,
% in the same arrangement as those in thetaMat, based
% solely on the SV-to-rRxAppx range evaluated at
% receiver time, in cycles.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author: Kendall Whitbeck
%+==============================================================================+
navConstants; % useful constants for navigation

% Guess for state vector xBar
cdtRx = 0;
xBar = [rRxAppx; cdtRx];
xBarNew = xBar;
error = 9001*ones(3,1);

% Pre-allocating variables for speed
solutionMat = [];
tRVecSolution.w = [];
tRVecSolution.s = [];
tTrueVecSolution.w = [];
tTrueVecSolution.s = [];
residualVec = [];
badSvIdVec = [];
thetaNominalMat = [];
P_covar = [];
Ip_vec = [];
T_vec = [];
dtrel_vec = [];
dXspin_vec = [];

[Nt, Nsv] = size(prMat); % [number of time epochs, number of sats in view]
Ns = 0; % number of solutions (nominally: Ns=iiStop-iiStart+1)
%% Double for-loop extracting valid pseudorange measurements
% % 1st for-loop goes through time epoch
% % 2nd for-loop goes through tracked sats at each time epoch
for ii = iiStart:iiStop
    
    % Reset number of valid satellites
    ValidSats = 0;
    
    % Storing signal tx time data in sensible variable names
    gpsWeek = tRVec.w(ii);
    gpsSec = tRVec.s(ii);
    
    %% Determining Valid Satellites & Number of Solutions
    % For-loop generating list of valid sats
    for jj = 1:Nsv
        
        % SVID number at element jj
        SVID = svIdVec(jj);
        
        if obsValidMat(ii,jj) == 1 && sum(ismember(svIdAllow,SVID)) >= 1
            
            % Increase count of valid SVIDs
            ValidSats = ValidSats + 1;
            
            % List of valid SVIDs used in zPrime. Each row of ValidSatList
            % corresponds to the same row in zPrime.
            ValidSVIDList(ValidSats,1) = SVID;                  %#ok<AGROW>
            
            % Stacking all pseudorange measurements valid at epoch tR
            % in a column vector
            zPrime(ValidSats,1) = prMat(ii,jj);                 %#ok<AGROW>
            
            % Counts num of solns so far (nominally, Ns=iiStop-iiStart+1,
            % but if some sats aren't valid, it'll be less)
            if ValidSats == 4 % need at least 4 valid sats for a soln
                Ns = Ns + 1;
            end%if
        end%if
    end%for
    
    %% Determining State Vector
    % Condition ensuring enough valid sats were overhead at time epoch ii
    % to get a solution
    if ValidSats < 4 && ii == iiStop && Ns == 0
        error('No solutions found for this data set, timespan, & elevation mask.')
    elseif ValidSats < 4 % change back to '<'!!!!!!!!!!!!!!
        continue
    end%if
    
    % while-loop generating state vector (xBar) solution for time epoch ii
    while error > 1e-4*ones(3,1)

        % for-loop creating new zBarPrime column vec & HPrime matrix
        % for current state vec, xBar
        for kk = 1:ValidSats
            
            % Pulling the next valid SVID
            SVID = ValidSVIDList(kk);
            
            % Getting predicted pseudorange, rhoBar, for each measured
            % pseudorange & storing in zBarPrime. Also getting Jacobian row
            % vec & storing into Jacobian matrix, HrhoPrime
            [zBarPrime(kk,1), HPrime(kk,:), ~, T(kk,1), Ip(kk,1),~,dtrel(kk,1),~,dXspin(kk,1)] = ...
                satpr(gpsWeek,gpsSec,xBar(4),xBar(1:3),satdata(SVID),ionodata,tiFlags);
        
        end%for
%         zBarPrime = zBarPrime + 1000;
        % Maximum residual at time epoch ii
        [maxResidualVal, maxResidualIndex] = max(abs(zBarPrime - zPrime));
        
        % Cholesky Factorization
        sigmaRho = 2;
        P_wPrime = sigmaRho^2 * eye(ValidSats);
        Ra = chol(P_wPrime);
        
        % Standard form for starting point in linear least-squares
        zBrevePrime = zPrime - zBarPrime + HPrime*xBar;
        
        % Normalizing the measurement
        z = inv(Ra)' * zBrevePrime;
        H = inv(Ra)' * HPrime;
                    
        % Covariance Matrix
        P_covar(:,:) = inv(H'*H);

        % Least-Squares Solution at time tR
        xBarNew = inv(H'*H) * H' * z;
        
        % Updating error tolerance
        error = abs(xBarNew(1:3) - xBar(1:3));
        
        % Updating value of xBar
        xBar = xBarNew;
        
    end%while

    %% Solution Updating
    % Returning an error if no solutions are found
    if Ns == 0
        error('No solutions found for this data set, timespan, & elevation mask.')
    end
        
    % Updating solutionMat
    solutionMat(Ns,:) = xBarNew';                               %#ok<AGROW>
    
    % Updating residualVec w/ max resiudal at each time epoch
    residualVec(Ns,1) = maxResidualVal;                       %#ok<AGROW>
    
    % Updating badSvIdVec w/ SVID w/ max resiudal at each time epoch
    badSvIdVec(Ns,1) = ValidSVIDList(maxResidualIndex);         %#ok<AGROW>
    
    % Updating tRVecSolution for GPS weeks and seconds
    tRVecSolution.w(Ns,1) = tRVec.w(ii);
    tRVecSolution.s(Ns,1) = tRVec.s(ii);
    
    % Updating tTrueVecSolution for GPS weeks and seconds
    tTrueVecSolution.w = tRVecSolution.w;
    tTrueVecSolution.s(Ns,1) = tRVec.s(ii) + xBarNew(4)/cLight;
    
    % Perturbations
    T_vec = [T_vec; T];
    Ip_vec = [Ip_vec; Ip];
    dtrel_vec = [dtrel_vec; dtrel];
    dXspin_vec = [dXspin_vec; dXspin];
end%for
end%function