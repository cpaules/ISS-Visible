% Chase Paules
% AME 457 
% Final Project

% General idea:
% compute r_site and rho, propagate orbit and in each time step check if
% the ISS is visible and write if it is

clear
clc
format long

% open a file named output.txt
fileID = fopen('output.txt', 'w');
% write header line
fprintf(fileID, '%6s %6s %3s %3s %12s %12s %12s %6s\n', 'DAY', 'HR', 'MIN', 'SEC', 'RHO (KM)', 'AZ (DEG)', 'EL (DEG)', 'VIS');
    
% Constants
Re = 6378.14;
ee = 0.081819221456;  
mu = 398600.4418;

% Converting yrday.fracday to time
qx = 15337.49034855;
qy = fix(qx/1000);
[yr, mo, day, hr, min, sec] = datevec(datenum(qy,0,qx - qy*1000));
yr = yr + 2000;
% Initial UT1
UT1 =  hr + min/60 + sec/3600;
% Compute GMST
[thetaG] = compute_GMST(mo, day, yr, hr, min, sec, UT1);

% Sun vector
[rasc, decl, rsun] =  sun(juliandate(1990,12,1));

% Site Input Data
lat = 32.2317;
long = -110.9519;
alt = .71;

% ISS Data
e = .0008187;
incl = 51.6405*pi()/180;
w = 236.0206*pi()/180;
RA = 319.1241*pi()/180;
M = 124.0169*pi()/180;
m = 15.54399752; %mean motion in rev/day, convert to rad/s
m = m*2*pi()/86400;
% Converting mean motion to semi major axis
a = (mu/m^2)^(1/3);

% Converting ISS mean anomaly(M) to eccentric anomaly(E) to true
% anomaly(TA)
E = kepler_E(e, M);
sinTA = sin(E)*sqrt(1 - e^2)/(1 - e * cos(E));
cosTA = (cos(E) - e)/(1 - e * cos(E));
TA = atan2(sinTA, cosTA);

% Initial orbital elements, coe = [a e incl w RA TA]
coe = [a e incl w RA TA];

% Convert orbital elements to ECI position 
[rsat, v] = sv_from_coe(coe);
% Convert ECI position to ECEF position
[rsatECEF] = rsatToECEF(thetaG, rsat);
rsatECEF = rsatECEF';

% LST at site
theta = thetaG + long;
%=============== Compute ECEF site vector
Cearth = Re/(sqrt(1 -ee*ee*sind(lat)*sind(lat)));
Searth = Cearth*(1 - ee*ee);
rdel = (Cearth + alt)*cosd(lat);
rk = (Searth + alt)*sind(lat);
rsiteECEF = [rdel*cosd(long); rdel*sind(long); rk];
%================ Line of sight vector
rhoECEF = rsatECEF - rsiteECEF;
%================ Line of sight vector in SEZ frame
[rhoSEZ] = rhoECEFToSEZ(long, lat, rhoECEF);    
%============= angles of ISS in degrees
azimuth = atan2d(rhoSEZ(2), -rhoSEZ(1));
elevation = atan2d(rhoSEZ(3), sqrt(rhoSEZ(1)^2 + rhoSEZ(2)^2));

%============= Visibility Conditions
darkness = dot(rsun,rsiteECEF);
above_horizon = rhoSEZ(3);
psi = acosd(dot(rsun,rsat));

if (darkness < 0) && (above_horizon > 0) && (abs(norm(rsat)*sind(psi)) > Re)
    vis = 'YES';
else
    vis = 'NO';
end

% write a line of output
fprintf(fileID, '%6.0f %6.0f %3.0f %3.0f %12.3f %12.3f %12.3f %6s\n', day, hr, min, sec, norm(rhoECEF), azimuth, elevation, vis);


n = sqrt(mu/(a^3));
dt = 60;
counter = 2;

[RA_dot, w_dot] = oblateness(e, a, incl);

for j=1:1439;

    % Update M
    M = M + n*dt; 
    % Update E
    E = kepler_E(e,M);
    % Update true anomaly
    sinTA = sin(E)*sqrt(1 - e^2)/(1 - e * cos(E));
    cosTA = (cos(E) - e)/(1 - e * cos(E));
    TA = atan2(sinTA, cosTA);
    % Update RA and w
    RA = RA + RA_dot;
    w = w + w_dot;
    % Updated orbital elements to r and v
    coe(6) = TA;
    coe(5) = RA;
    coe(4) = w;
    [rsat, v] = sv_from_coe(coe);
    % Convert sat ECI position to ECEF position
    [rsatECEF] = rsatToECEF(thetaG, rsat);
    rsatECEF = rsatECEF';
    
    % Compute new UT1 time
    UT1 = UT1 + (dt/3600);
    % Compute new GMST
    [thetaG] = compute_GMST(mo, day, yr, hr, min, sec, UT1);
       
    % LST at site
    theta = thetaG + long;
    %=============== Compute ECEF site vector
    Cearth = Re/(sqrt(1 -ee*ee*sind(lat)*sind(lat)));
    Searth = Cearth*(1 - ee*ee);
    rdel = (Cearth + alt)*cosd(lat);
    rk = (Searth + alt)*sind(lat);
    rsiteECEF = [rdel*cosd(long); rdel*sind(long); rk];
    %================ Line of sight vector
    rhoECEF = rsatECEF - rsiteECEF;
    %================ Line of sight vector in SEZ frame
    [rhoSEZ] = rhoECEFToSEZ(long, lat, rhoECEF);          
    %============= angles of ISS in degrees
    azimuth = atan2d(rhoSEZ(2), -rhoSEZ(1));
    elevation = atan2d(rhoSEZ(3), sqrt(rhoSEZ(1)^2 + rhoSEZ(2)^2));
    
    %============= Visibility Conditions
    darkness = dot(rsun,rsiteECEF);
    above_horizon = rhoSEZ(3);
    psi = acosd(dot(rsun,rsat));

    if (darkness < 0) && (above_horizon > 0) && (abs(norm(rsat)*sind(psi)) > Re)
        vis = 'YES';
    else
        vis = 'NO';
    end
    
    % write a line of output
    fprintf(fileID, '%6.0f %6.0f %3.0f %3.0f %12.3f %12.3f %12.3f %6s\n', day, hr, min, sec, norm(rhoECEF), azimuth, elevation, vis);
    
    % Converting yrday.fracday to time
    qx = qx + dt/86400;
    qy = fix(qx/1000);
    [yr, mo, day, hr, min, sec] = datevec(datenum(qy,0,qx - qy*1000));
    yr = yr + 2000;
    
    counter = counter + 1;
end

% close all open files
status = fclose('all');