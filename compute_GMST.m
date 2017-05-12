function [theta_G] = compute_GMST(mo, d, y, hr, min, sec, UT1)

% ---------------------------------------------
%{
This function computes Greenwich Mean Sidereal Time (GMST) using equations
from slide 32 of the lesson_17.pdf on D2L

mu - GM of Earth 398,600 (km^3/s^2)
J2 - Second zonal harmonic of the Earth
R - Radius of the Earth 6378 (km)
e - eccentricity (magnitude of E)
a - semimajor axis (km)
incl - inclination of the orbit (rad)
RA - right ascension of the ascending node (rad)
w - argument of perigee (rad)
UT1 - in hours
User M-functions required: None

%}
% ---------------------------------------------

j0 = juliandate(y,mo,d,hr,min,sec);
T0 = (j0 - 2451545)/36525;
theta_G0 = 100.460184 + 36000.77004*T0 + .000387933*T0^2 - (2.3583e-8)*T0^3;
theta_G = theta_G0 + 360.98564724*(UT1/24);
theta_G = wrapTo360(theta_G);

end

