function [ rsiteECI ] = sitepositionECI(lat, long, alt, thetaG)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Re = 6378.14;
ee = 0.081819221456;  

%=============== Compute ECEF site vector
Cearth = Re/(sqrt(1 -ee*ee*sind(lat)*sind(lat)));
Searth = Cearth*(1 - ee*ee);
rdel = (Cearth + alt)*cosd(lat);
rk = (Searth + alt)*sind(lat);
rsiteECEF = [rdel*cosd(long); rdel*sind(long); rk];
%================ Compute ECI site vector
ct = cosd(-thetaG);
st = sind(-thetaG);
ecef2eci = [ct st 0; -st ct 0; 0 0 1];
rsiteECI = ecef2eci * rsiteECEF;
end

