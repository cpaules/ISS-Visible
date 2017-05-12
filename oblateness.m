function [RA_dot, w_dot] = oblateness(e, a, incl)

% ---------------------------------------------
%{
This function computes the rate of change of RA (rad/s) and arguement of 
perigee (rad/s) due to Earth's oblateness using equations 4.52 and 4.54 
from Curtis

mu - GM of Earth 398,600 (km^3/s^2)
J2 - Second zonal harmonic of the Earth
R - Radius of the Earth 6378 (km)
e - eccentricity (magnitude of E)
a - semimajor axis (km)
incl - inclination of the orbit (rad)
RA - right ascension of the ascending node (rad)
w - argument of perigee (rad)
User M-functions required: None

%}
% ---------------------------------------------

% Costants
mu = 398600.4418;
J2 = .00108263;
R = 6378.14;

% Change of RA calculation (rad/s)
RA_dot = -((3*sqrt(mu)*J2*R^2)/(2*(1-e^2)^2*a^(7/2)))*cos(incl);

% Change of w calculation (rad/s)
w_dot = RA_dot*((5/2)*(sin(incl)^2)-2)/cos(incl);

end

