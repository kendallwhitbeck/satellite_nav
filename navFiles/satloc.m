function [rSvEcef, vSvEcef, Ek] = satloc(gpsWeek, gpsSec, sd)
% satloc : Return satellite location and velocity expressed in and relative to
% the ECEF reference frame.

% INPUTS
%
% gpsWeek ---- Week of true time at which SV location and velocity are
% desired.
%
% gpsSec ----- Seconds of week of true time at which SV location & velocity
% are desired.
%
% satdata ------- Ephemeris structure array for a single SV. Let ii be the
% numerical identifier (PRN identifier) for the SV whose location
% is sought. Then sd = satdata(ii). sd has the following
% fields:
%
% SVID - satellite number
% health - satellite health flag (0 = healthy; otherwise unhealthy)
% we - week of ephemeris epoch (GPS week, unambiguous)
% te - time of ephemeris epoch (GPS seconds of week)
% wc - week of clock epoch (GPS week)
% tc - time of clock epoch (GPS seconds of week)
% e - eccentricity (unitless)
% sqrta - sqrt of orbit semi-major axis (m^1/2)
% omega0 - argument of perigee (rad.)
% M0 - mean anomaly at epoch (rad.)
% L0 - longitude of ascending node at beginning of week (rad.)
% i0 - inclination angle at epoch (rad.)
% dOdt - longitude rate (rad / sec.)
% dn - mean motion difference (rad / sec.)
% didt - inclination rate (rad / sec.)
% Cuc - cosine correction to argument of perigee (rad.)
% Cus - sine correction to argument of perigee (rad.)
% Crc - cosine correction to orbital radius (m)
% Crs - sine correction to orbital radius (m)
% Cic - cosine correction to inclination (rad.)
% Cis - sine correction to inclination (rad.)
% af0 - 0th order satellite clock correction (s)
% af1 - 1st order satellite clock correction (s / s)
% af2 - 2nd order satellite clock correction (s / s^2)
% TGD - group delay time for the satellite (s)

% OUTPUTS
%
% rSvEcef ---- Location of SV at desired time expressed in the ECEF reference
% frame (m).
%
% vSvEcef ---- Velocity of SV at desired time relative to the ECEF reference
% frame and expressed in the ECEF reference frame (m/s). NOTE:
% vSvEcef is NOT inertial velocity, e.g., for geostationary SVs
% vSvEcef = 0.

%+------------------------------------------------------------------------------+
% References:
% https://www.gps.gov/technical/icwg/IS-GPS-200J.pdf
% https://fenrir.naruoka.org/download/autopilot/note/080205_gps/gps_velocity.pdf
% Author: Kendall Whitbeck
%+==============================================================================+
% Physical constants used for navigation
navConstants; % GM, OmegaE

% Semi-major axis of satellite
A = (sd.sqrta)^2; % meters

% Computed mean motion
n0 = sqrt(GM / A^3); % rad/sec

% time from ephemeris reference epoch where gpsSec is the desired time & te
% is the reference time ephemeris
tk = gpsSec - sd.te; % sec

% Checking for tk ambiguity
if tk > 302400
    tk = tk - 604800; % sec
elseif tk < -302400
    tk = tk + 604800; % sec
end

% Corrected Mean Motion
n = n0 + sd.dn; % rad/sec

% Mean Anomaly & its Time-Derivative
Mk = sd.M0 + n*tk; % rad
Mk_dot = n; % rad/sec

% Eccentric Anomaly (rad)
Ek0 = Mk; % initial guess
Ek_tol = 1e-9; % tolerance for Ek
for ii = 1:10
    Ek = Mk + sd.e * sin(Ek0); % calculates Ek
    if abs(Ek-Ek0) <= Ek_tol
        break; % breaks loop if Ek is suffciently accurate
    end%if    
    Ek0 = Ek; % reassigns old Ek ("Ek0")
end%for

% Time-Derivative of Eccentric Anomaly
Ek_dot = Mk_dot/(1-sd.e*cos(Ek)); % rad/sec

% True Anomaly (rad)
[tru_k, ~] = Me2tru(sd.e,Mk);

% Time-Derivative of True Anomaly (rad/sec)
tru_k_dot = sin(Ek)*Ek_dot*(1+sd.e*cos(tru_k))/((1-cos(Ek)*sd.e)*sin(tru_k));

% Argument of Latitude & its Time-Derivative
phi_k = tru_k + sd.omega0; % rad
phi_k_dot = tru_k_dot; % rad/sec

% Second Harmonic Perturbations
del_uk = sd.Cus*sin(2*phi_k) + sd.Cuc*cos(2*phi_k); % Arg of Lat Correction (rad)
del_rk = sd.Crs*sin(2*phi_k) + sd.Crc*cos(2*phi_k); % Radius Correction (meters)
del_ik = sd.Cis*sin(2*phi_k) + sd.Cic*cos(2*phi_k); % Inclination Correction (rad)

% Time-Derivatives of Second Harmonic Perturbations
del_uk_dot = 2*(sd.Cus*cos(2*phi_k) - sd.Cuc*sin(2*phi_k))*phi_k_dot; % rad/sec
del_rk_dot = 2*(sd.Crs*cos(2*phi_k) - sd.Crc*sin(2*phi_k))*phi_k_dot; % m/s
del_ik_dot = 2*(sd.Cis*cos(2*phi_k) - sd.Cic*sin(2*phi_k))*phi_k_dot; % rad/sec

% Corrected Topographic Coordinates
uk = phi_k + del_uk; % Corrected Argument of Latitude (rad)
rk = A*(1 - sd.e*cos(Ek)) + del_rk; % Corrected Radius (meters)
ik = sd.i0 + sd.didt*tk + del_ik; % Corrected Inlination (rad)

% Time-Derivative of Corrected Topographic Coordinates
uk_dot = phi_k_dot + del_uk_dot; % Time-Derivative of Corrected Argument of Latitude (rad/sec)
rk_dot = A*sd.e*sin(Ek)*Ek_dot + del_rk_dot; % Time-Derivative of Corrected Radius (m/s)
ik_dot = sd.didt + del_ik_dot; % Time-Derivative of Corrected Inlination (rad/sec)

% Positions in orbital plane
xkP = rk*cos(uk); % meters
ykP = rk*sin(uk); % meters

% Time-Derivative of positions in orbital plane
xkP_dot = rk_dot*cos(uk) - ykP*uk_dot; % m/s
ykP_dot = rk_dot*sin(uk) + xkP*uk_dot; % m/s

% Corrected longitude of ascending node & its time-derivative
Omegak = sd.L0 + (sd.dOdt - OmegaE)*tk - OmegaE*sd.te; % rad
Omegak_dot = sd.dOdt - OmegaE; % rad/sec

% Earth Fixed Coordinates
xk = xkP*cos(Omegak) - ykP*cos(ik)*sin(Omegak); % meters
yk = xkP*sin(Omegak) + ykP*cos(ik)*cos(Omegak); % meters
zk = ykP*sin(ik); % meters

% Time-Derivative of Earth Fixed Coordinates
xk_dot = xkP_dot*cos(Omegak) - ykP_dot*cos(ik)*sin(Omegak)...
    + ykP_dot*sin(ik)*sin(Omegak)*ik_dot - yk*Omegak_dot; % m/s
yk_dot = xkP_dot*sin(Omegak) + ykP_dot*cos(ik)*cos(Omegak)...
    - ykP_dot*sin(ik)*ik_dot*cos(Omegak) + xk*Omegak_dot; % m/s
zk_dot = ykP_dot*sin(ik) + ykP*cos(ik)*ik_dot; % m/s

% Location of SV at desired time expressed in the ECEF reference frame (m)
rSvEcef = [xk; yk; zk]; % meters

% Velocity of SV at desired time expressed in the ECEF reference frame (m)
vSvEcef = [xk_dot; yk_dot; zk_dot]; %  % m/s

end