function [solutionMat,tRVecSolution,tTrueVecSolution,...
residualVec,badSvIdVec,thetaNominalMat, iterStore, covar, delayMat, ionoL2] = ...
performNavigationSolution(tRVec,svIdVec,obsValidMat,prMat,fDMat,thetaMat,...
iiStart,iiStop,rRxAppx, ...
satdata,ionodata,svIdAllow,lambda,tiFlags)
% performNavigationSolution : Perform the navigation solution using nonlinear
% least squares techniques.
%
%
% INPUTS
%
% tRVec --------------- Structure of two Nt-by-1 vectors that jointly define
% the unique receiver time instants corresponding to all
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
% tiFlags ------------- 2-by-1 vector of flags that indicate whether to enable
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
% an estimate of 'true' GPS time for each
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
%   Miller, I. "getTropoDelay.m" "getIonodelay.m", 2006
%   Global Positioning Systems Directorate Systems Engineering & Integration
%       Interface Specification IS-GPS-200J, 25 Apr 2018
%   Remondi, B., "Computing satellite velocity using the broadcast ephemeris,"
%   GPS Solutions, vol. 8 no. 3, 2004
%   Misra, P., Enge, P., "Global Positioning System, 2nd Ed."
%
%
% Author: Benjamin Miller
%+==============================================================================+

%% trim input variables
    %tRVec = tRVec(iiStart:iiStop); %limit time search

    xbarMat = zeros(length(tRVec.s), 4); %instantiate curve fit coefficients
    errorVec = zeros(iiStop-iiStart+1,1);
    badSV = zeros(iiStop-iiStart+1,1);
    svIdVecOld = svIdVec;
    iterStore = zeros(6,3);
    iterCount = 1;
    covar = zeros(4,4,iiStop-iiStart+1);
    delayMat = [];
    ionoL2 = zeros(iiStop-iiStart+1,2);
%% instantiate time loop
    for i = iiStart:iiStop %loop all time epochs
        
%% Initiate guess
        xbarold = [rRxAppx; 0]; %initial state guess
  
%% stack all pseudorange measurements valid for tr into z
        zprime = prMat(i,:)'; %pull pseudorange measurements
        svIdVecLocal = svIdVec;
        for j = 1:length(zprime) 
            if obsValidMat(i,j)==0 %exclude non-valid measures
                zprime(j)=0;
                svIdVecLocal(j) = 0;
            end
            if ismember(svIdVecOld(j),svIdAllow)==0
                zprime(j)=0;
                svIdVecLocal(j) = 0;
            end
        end
        zprime = nonzeros(zprime);
 
        %account for low-SV instances
        if size(zprime,1)<4
            zprime = prMat(i,:)'; %pull pseudorange measurements
            svIdVecLocal = svIdVec;
            for j = 1:length(zprime) 
                if obsValidMat(i,j)==0 %exclude non-valid measures
                    zprime(j)=0;
                    svIdVecLocal(j) = 0;
                end
                %allow all SVs
            end
        end
        
        zprime = nonzeros(zprime);
        svIdVecLocal = nonzeros(svIdVecLocal);
%% instantiate curve-fitting loop 
    tol = 1; %instantiated to start loop
    while tol>10^(-8)
    
%% calculate predicted pseudoranges and Jacobians
            Hprime = zeros(length(svIdVecLocal),4); %instantiate jacobian
            zbarprime = zeros(length(svIdVecLocal),1); %instantiate calculations
            for j = 1 : length(svIdVecLocal)
                sd = satdata(svIdVecLocal(j));
                [r, H, iono, tropo, rel] = ...
                    satpr(tRVec.w(i), tRVec.s(i),xbarold(4), xbarold(1:3), satdata(svIdVecLocal(j)), ionodata, tiFlags, 1);
                    %calculates predicted pseudoranges from model
                    [rOff, ~, ~, ~, ~] = ...
                    satpr(tRVec.w(i), tRVec.s(i),xbarold(4), xbarold(1:3), satdata(svIdVecLocal(j)), ionodata, tiFlags, 0);
                Hprime(j,:) = H;
                zbarprime(j) = r;
                delayMat = [delayMat; iono, tropo, rel, abs(r-rOff)];
                if svIdVecLocal(j)==14
                    ionoL2(i,1) = iono;
                elseif svIdVecLocal(j) ==25
                    ionoL2(i,2) = iono;
                    
                end
                    
            end
            dz = (zprime - zbarprime) + (Hprime*xbarold); %err + model = measure
%% Instantiate weighting matrix
            sigma = 2; %predicted gausian distribution of rho
            Pa = sigma^2 .* eye(size(dz,1)); % predicted covariance 
            Ra = (chol(Pa)); %Cholesky-factorized weighting matrix of Pa
    
%% normalize measurements
            z = Ra'\dz; %normalized
            H = Ra'\Hprime; %normalized Jacobian
            covar(:,:,i) = inv(H'*H);
%% mimimize cost function with QR factorization   
            [Q, R] = qr(H);
            Rest = R(1:4, :);
            zest = Q'*z;
            zest = zest(1:4);
            xbarnew = linsolve(Rest,zest); 
            tol = norm(xbarnew-xbarold); %position error
            if i == iiStart
                if iterCount<6
                    iterStore(iterCount,:) = xbarold(1:3);
                    iterStore(iterCount+1,:) = xbarnew(1:3);
                end
                iterCount = iterCount+1;
            end
            xbarold = xbarnew;

    end %escape while loop for xbarnew
    xbarMat(i,:) = xbarold; %save outside of time loop
    [errorVec(i), badSVindex] = max(abs((zprime - zbarprime)));
    badSV(i) = svIdVecLocal(badSVindex);
    end %escape time loop

%% define output 
    thetaNominalMat = []; %undefined for current algorithm
    solutionMat = zeros(iiStop-iiStart+1,4); %trim for topsol.m
    for i=1:4
        solutionMat(:,i) = nonzeros(xbarMat(:,i));
    end
    residualVec = [];
    for i = 1:size(badSV,1)
        if badSV(i) ~= 0
            residualVec = [residualVec ; errorVec(i)];
        end
    end
    badSvIdVec = nonzeros(badSV(:));
    tRVecSolution.s = tRVec.s(iiStart:iiStop);
    tRVecSolution.w = tRVec.w(iiStart:iiStop);
    tTrueVecSolution.w = tRVecSolution.w;
    tTrueVecSolution.s = tRVecSolution.s+(solutionMat(:,4)./299792458);
                %disp(svIdVecLocal)
    ionoL2 = [nonzeros(ionoL2(:,1)) , nonzeros(ionoL2(:,2))];
    delayMat = mean(delayMat); %returns row vector containing mean of each column
        %[ionospheric, tropospheric, relitivistic, earth rotation]
end