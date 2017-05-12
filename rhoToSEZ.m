function [rhoSEZ] = rhoToSEZ(long, lat, thetaG, rhoECI)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

ct2 = cosd(long + thetaG);
st2 = sind(long + thetaG);
r3 = [ct2 st2 0; -st2 ct2 0; 0 0 1];
x = 90 - lat;cx = cosd(x);
sx = sind(x);
r2 = [cx 0 -sx; 0 1 0; sx 0 cx];
eci2sez = r2*r3;
rhoSEZ = eci2sez*rhoECI;

end

