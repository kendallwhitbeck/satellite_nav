% Me2Tru : Convert Mean Anomaly to True Anomaly using Newton's Method

% INPUTS
%
% ecc ---- Eccentricity of the satellite in orbit.
%
% Me ---- Mean Anomaly of the satellite in orbit.

% OUTPUTS
%
% tru ---- True Anomaly of the satellite in orbit.
%
% Ef ---- Eccentric Anomaly of the satellite in orbit.

%+------------------------------------------------------------------------------+
% References:
% Dr. Brandon Jones, ASE 366K: Spacecraft Dynamics
%
% Author: Kendall Whitbeck
%+==============================================================================+

function [tru,Ef] = Me2tru(ecc,Me)

% Determining Initial Guess for Eccentric Anomaly, E0
if Me < pi
    E0 = Me + ecc/2;
else
    E0 = Me - ecc/2;
end

% Using Newton-Raphson's to attain more accurate Eccentric Anomaly
for ii = 1:1000
    
    % Ensure that Eccentric Anomaly remains under 2*pi radians
    while E0 > 2*pi
        E0 = E0 - 2*pi;
    end
    
    % Calculate new iteration of Eccentric Anomaly guess
    Ef = E0 - (E0 - ecc*sin(E0) - Me) / (1 - ecc*cos(E0));
    
    % Check for convergence
    if abs(Ef - E0) < 10^(-8)
        break
    end
    
    % Update Eccentric Anomaly guess
    E0 = Ef;
    
    if ii == 10
        error('Eccentric Anomaly Does NOT Converge after 10 Iterations')
    end
end

% Using Eccentric anomaly to calculate True Anomaly
tru = 2 * atan2(sqrt(1 + ecc) * tan(Ef / 2), sqrt(1 - ecc));

end