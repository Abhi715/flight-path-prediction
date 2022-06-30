%{
*****************************************************************************
% NOTE: backup of this script is saved in file GnssCons.m
*****************************************************************************
%}

clc;
clear all;
close all;

% variables 
rad = pi/180;
M  = 5.972 * 10^24;                                                         % mass of earth in kilogram
R = 6371;                                                                   % radius of earth in Kilometers (6378)
G = 6.673 * 10^-20;						                                    % gravitational constant in km^3 kg^-1 sec^-2
mu = 398600.4405;                                                           % km^3/sec^2

C = [1, 0; 0, 0];
S = [0, 0; 0, 0];

C_bar = [0, 0, 0, 0, 0, 0; 
        -484.165368*10^-6, -0.000187*10^-6, 2.439261*10^-6, 0, 0, 0; 
        0.957171*10^-6, 2.030137*10^-6, 0.904706*10^-6, 0.721145*10^-6, 0, 0; 
        0.539777*10^-6, -0.536244*10^-6, 0.350670*10^-6, 0.990869*10^-6, -0.188481*10^-6, 0; 
        0.068659*10^-6, -0.062727*10^-6, 0.652459*10^-6, -0.451837*10^-6, -0.295123*10^-6, 0.174832*10^-6];

S_bar = [0, 0, 0, 0, 0, 0; 
        0.000000, 0.001195*10^-6, -1.400266*10^-6, 0, 0, 0; 
        0.000000, 0.248131*10^-6, -0.618923*10^-6, 1.414204*10^-6, 0, 0; 
        0.000000, -0.473772*10^-6, 0.662571*10^-6, -0.200987*10^-6, 0.308848*10^-6, 0; 
        0.000000, -0.094195*10^-6, -0.323334*10^-6, -0.214954*10^-6, 0.049741*10^-6, -0.669393*10^-6];

syms u latitude range longitude theta
sum_a1 = 0;
sum_a2 = 0;
acceleration_altitude = [];
acceleration_latitude = [];
acceleration_longitude = [];
motion = mu/range^2;
degree = 5;                         % degree = n

% TLE data (Satellite ephhemeris data)
catalog = [];
fd = fopen("SpiraleB.txt",'r');
a = 0;
L0 = fgetl(fd);
L1 = fgetl(fd);
L2 = fgetl(fd);
while ischar(L2)
    a = a + 1;
    satnum = str2num(L1(3:7));
    if isempty(catalog) || ismember(satnum, catalog)
        fprintf('%s\n', repmat('-',1,50));
        fprintf('Satellite: %s\n', L0);
        assert(chksum(L1), 'Checksum failure on line 1')
        assert(chksum(L2), 'Checksum failure on line 2')
        fprintf('Catalog Number: %d\n', satnum)
        fprintf('Epoch time: %s\n', L1(19:32))                            % Epoch Time YYDDD.DDDDDDDD
        Incl = str2num(L2(9:16));                                         % Inclination
        fprintf('Inclination: %f deg\n', Incl)
        Omega = str2num(L2(18:25));                                       % Ra of accencding node (RAAN)
        fprintf('RA of ascending node: %f deg\n', Omega)
        ecc = str2num(['.' L2(27:33)]);                                   % Eccentricity
        fprintf('Eccentricity: %f\n', ecc)
        w = str2num(L2(35:42));                                           % argument of perigee or periapsis (perigee argument)
        fprintf('Arg of perigee: %f deg\n', w)        
        M = str2num(L2(44:51));                                           % mean anomaly
        fprintf('Mean anomaly: %f deg\n', M)
        b = str2num(L2(53:63));                                           % mean motion
        fprintf('Mean motion: %f rev/day\n', b)
        T = 86400/b;
        fprintf('Period of rev: %.0f s/rev\n', T)
		
        a = ((T/(2*pi))^2*398.6e12)^(1/3);                                % semi-major axis            
        fprintf('Semi-major axis: %.0f meters\n', a)
		
        b = a*sqrt(1-ecc^2);                                              % semi-minor axis
        fprintf('Semi-minor axis: %.0f meters\n', b);
        
		mean_anomaly = M*(pi/180);                                        % in radians
        v = (mean_anomaly + (2*ecc*(sin(mean_anomaly))));                 % in radians
        true_anomaly = v*(180/pi);
        fprintf('true anomaly: %f deg\n', true_anomaly)                   % true anomaly in 
    end
    L0 = fgetl(fd);
    L1 = fgetl(fd);
    L2 = fgetl(fd);
end

% Update the code for GNSS 
% Velocity, position, ground station, azimuth, elevation and range
startT_1 = datetime(2000, 01, 01, 12, 00, 00);
stopT_1 = datetime(2022, 01, 01, 12, 00, 00);
sampleT_1 = 60;
sc_1 = satelliteScenario(startT_1, stopT_1, sampleT_1);                % satellite scenario 
% TLE files, azimuth, elevation and range
TleFile_1 = "SpiraleB.txt";                                           % TLE file 1                                           
sat_1 = satellite(sc_1, TleFile_1);
%show(sat_1);                                                           % show the simulation of satellite orbit
%play(sc_1);
gs_1 = groundStation(sc_1);                                             % details about ground station
%sat_1_1 = satellite(sc_1, a, ecc, Incl, Omega, w, true_anomaly);
%show(sat_1_1);

%{
*****************************************************************************
% accelaration vector for altitude 
*****************************************************************************
%}

for n = 0:degree
    % Associated Legrende Polynomial
    legendre_polynomial = ((1/((2^n) *(factorial(n)))) * (diff((((u^2)-1)^n), u, n)));                 % legendre polynomial of degree n in terms of variable
    sum_a1 = (((n+1)*(R/range)^n) + sum_a1);
    for m = 0:n
        associated_legendre_poly = (((1-u^2)^(m/2)) * diff(legendre_polynomial, u, m));         % disp(associated_legendre_poly)
        associated_legendre_poly = subs(associated_legendre_poly, u, sin(latitude));            % replacing the u with sin(phi) and solve the equation

        %{
        ****************** Kronecker Symbol (δ0m) if n==m then ***************
            It will be 1 when n and m both will be 0, otherwise it will 0
                            if m is 0 then it will be 1
        **********************************************************************
        %}  
        %if n==0 && n==m
        if m == 0
            delta = 1;
        else
            delta = 0;
        end

        % geopotential coefficients
        if n < 2
            geopotential_coefficient_C = C(n+1, m+1);
            geopotential_coefficient_S = S(n+1, m+1);
        else
            geopotential_coefficient_C = (C_bar(n, m+1))/(sqrt(factorial(n+m)/((2-delta) * ((2*n)+1) * factorial(n-m))));
            geopotential_coefficient_S = (S_bar(n, m+1))/(sqrt(factorial(n+m)/((2-delta) * ((2*n)+1) * factorial(n-m))));
        end
        
        sum_a2 = (associated_legendre_poly*((geopotential_coefficient_C *cos(m*longitude))+((geopotential_coefficient_S)*sin(m*longitude))))+sum_a2;
        acceleration_altitude = (sum_a1*sum_a2); 
    end
end
acceleration_altitude = -motion*acceleration_altitude;

%{
*****************************************************************************
% accelaration vector for latitude (phi) 
*****************************************************************************
%}

for n = 1:degree
    % Associated Legrende Polynomial
    legendre_polynomial = ((1/((2^n) *(factorial(n)))) * (diff((((u^2)-1)^n), u, n)));                 % legendre polynomial of degree n in terms of variable
    %acceleration_latitude
    sum_a1 = ((R/range)^n)+sum_a1;
    for m = 0:n
        associated_legendre_poly = (((1-u^2)^(m/2)) * diff(legendre_polynomial, u, m));        % disp(associated_legendre_poly);
        associated_legendre_poly = subs(associated_legendre_poly, u, sin(latitude));
        associated_legendre_poly = diff(associated_legendre_poly, latitude, 1);
        
        %{
        ****************** Kronecker Symbol (δ0m) if n==m then ***************
            It will be 1 when n and m both will be 0, otherwise it will 0
                            if m is 0 then it will be 1
        **********************************************************************
        %}  
        %if n==0 && n==m
        if m == 0
            delta = 1;
        else
            delta = 0;
        end

        % geopotential coefficients
        if n < 2
            geopotential_coefficient_C = C(n+1, m+1);
            geopotential_coefficient_S = S(n+1, m+1);
        else
            geopotential_coefficient_C = (C_bar(n, m+1))/(sqrt(factorial(n+m)/((2-delta) * ((2*n)+1) * factorial(n-m))));
            geopotential_coefficient_S = (S_bar(n, m+1))/(sqrt(factorial(n+m)/((2-delta) * ((2*n)+1) * factorial(n-m))));
        end
        sum_a2 = ((associated_legendre_poly * (geopotential_coefficient_C*cos(m*longitude) + (geopotential_coefficient_S*sin(m*longitude)))) + sum_a2);
        acceleration_latitude = sum_a1*sum_a2;
    end
end
acceleration_latitude = motion * acceleration_latitude;

%{
*****************************************************************************
% accelaration vector for longitude (lambda) 
*****************************************************************************
%}
for n = 1:degree
    % Associated Legrende Polynomial
    legendre_polynomial = ((1/((2^n) *(factorial(n)))) * (diff((((u^2)-1)^n), u, n)));                 % legendre polynomial of degree n in terms of variable
    sum_a1 = ((R/range)^n)+sum_a1;
    for m = 1:n
        associated_legendre_poly = (((1-u^2)^(m/2)) * diff(legendre_polynomial, u, m));        % disp(associated_legendre_poly);
        associated_legendre_poly = subs(associated_legendre_poly, u, sin(latitude));
        
        %{
        ****************** Kronecker Symbol (δ0m) if n==m then ***************
            It will be 1 when n and m both will be 0, otherwise it will 0
                            if m is 0 then it will be 1
        **********************************************************************
        %}  
        %if n==0 && n==m
        if m == 0
            delta = 1;
        else
            delta = 0;
        end

        % geopotential coefficients
        if n < 2
            geopotential_coefficient_C = C(n+1, m+1);
            geopotential_coefficient_S = S(n+1, m+1);
        else
            geopotential_coefficient_C = (C_bar(n, m+1))/(sqrt(factorial(n+m)/((2-delta) * ((2*n)+1) * factorial(n-m))));
            geopotential_coefficient_S = (S_bar(n, m+1))/(sqrt(factorial(n+m)/((2-delta) * ((2*n)+1) * factorial(n-m))));
        end
        sum_a2 = ((m*(associated_legendre_poly/cos(latitude))) * (-geopotential_coefficient_C*sin(m*longitude)+(geopotential_coefficient_S*cos(m*longitude))));
        acceleration_longitude = sum_a1*sum_a2;
    end
end
acceleration_longitude = motion*acceleration_longitude;

acceleration_ecef = [acceleration_altitude acceleration_latitude acceleration_longitude];                   % acceleration in ECEF reference frame

%{
*****************************************************************************
% convert from ECEF to ECI
*****************************************************************************
%}
acceleration_cartesian = [(cos(latitude)*cos(longitude)) -(sin(latitude)*cos(longitude)) -sin(longitude); (cos(latitude)*sin(longitude)) -(sin(latitude)*sin(longitude)) cos(longitude); sin(latitude) cos(latitude) 0] * transpose(acceleration_ecef);

a = 1;
while a < 5
    time_index = input("enter the index (more than 1): ");
    predict_index = input("Enter time index to predict (shouldn't more than previous index + 10 for accuracy): ");

    % extract the position, velocity, altitude, latitude, longitude and time from file
    x = readmatrix("position.xlsx", "Range", [time_index 1 time_index 1]);
    y = readmatrix("position.xlsx", "Range", [time_index 2 time_index 2]);
    z = readmatrix("position.xlsx", "Range", [time_index 3 time_index 3]);
    
    % actual future position
    real_r_x = readmatrix("position.xlsx", "Range", [predict_index 1 predict_index 1]);
    real_r_y = readmatrix("position.xlsx", "Range", [predict_index 2 predict_index 2]);
    real_r_z = readmatrix("position.xlsx", "Range", [predict_index 3 predict_index 3]);

    % current position in ECI reference frame
    r_eci = [x y z];
    
    vx = readmatrix("position.xlsx", "Range", [time_index 4 time_index 4]);
    vy = readmatrix("position.xlsx", "Range", [time_index 5 time_index 5]);
    vz = readmatrix("position.xlsx", "Range", [time_index 6 time_index 6]);
    
    alt = readmatrix("position.xlsx", "Range", [time_index 7 time_index 7]);            % current altitude 
    % add radius of eaerth to altitude (altitude is from the surface of earth and we need from the center of earth)
    alt = alt+R;

    phi = readmatrix("position.xlsx", "Range", [time_index 8 time_index 8]);            % current latitude (phi)
    % convert latitude from degree to radian
    phi = deg2rad(phi);
    
    lambda = readmatrix("position.xlsx", "Range", [time_index 9 time_index 9]);         % current longitude (lambda)
    % convert longitude from degree to radian
    lambda = deg2rad(lambda);

    current_date = readtable("position.xlsx", "Range", [time_index 10 time_index 10]);
    predicted_date = readtable("position.xlsx", "Range", [predict_index 10 predict_index 10]);
    
    % change current time to seconds
    current_date = [current_date.(1)];
    date = [current_date.Day current_date.Month current_date.Year];
    hour = int64(current_date.Hour);
    minute = int64(current_date.Minute);
    second = int64(current_date.Second);
    time = hour*3600 + minute*60 + second;                             % current time into seconds
    
    % predict date, time and utc are used to calculate rotation matrix from ECEF to ECI
    predicted_date = [predicted_date.(1)];
    d1 = int64(predicted_date.Day);
    m1 = int64(predicted_date.Month);
    y1 = int64(predicted_date.Year);

    predict_hour = int64(predicted_date.Hour);
    predict_minute = int64(predicted_date.Minute);
    predict_second = int64(predicted_date.Second);
    predict_time = predict_hour*3600 + predict_minute*60 + predict_second;                             % predict time into seconds
    utc = predict_time/3600;

    % Calculate interval (delta_t) to calculate predicted position 
    time_interval = predict_time - time;
    
    % substitute the values 
    acceleration_cartesian = subs(acceleration_cartesian, range, alt);
    acceleration_cartesian = subs(acceleration_cartesian, latitude, phi);
    acceleration_cartesian = subs(acceleration_cartesian, longitude, lambda);
    %disp(acceleration_cartesian);

    %a_x = acceleration_cartesian(1);
    %a_y = acceleration_cartesian(2);
    %a_z = acceleration_cartesian(3);
    
    % transform acceleration in ECEF to ECI frame
    dcm = angle(y1, m1, d1, utc); 
    acceleration_eci = dcm .* acceleration_cartesian;
    acceleration_eci = double(subs(acceleration_eci));
    
    time_interval = double(time_interval);
    % position 
    r_x = x + (vx + acceleration_eci(1) * time_interval) * time_interval;
    r_y = y + (vy + acceleration_eci(2) * time_interval) * time_interval;
    r_z = z + (vz + acceleration_eci(3) * time_interval) * time_interval;
    error_x = (abs(real_r_x - r_x)/real_r_x) * 100;
    error_y = (abs(real_r_y - r_y)/real_r_y) * 100;
    error_z = abs((real_r_z - r_z)/real_r_z) * 100;
    fprintf("predicted position in ECI reference frame r_x = %f, r_y = %f, r_z = %f\n", r_x, r_y, r_z);
    fprintf("percentage error in x, y and z direction respectively is %f %f %f\n", error_x, error_y, error_z);
    a=a+1;
end

%{
%fprintf("current position in ECI reference frame x = %f, y = %f, z = %f at time %f\n", x, y, z, time);
fprintf("predicted position in ECI reference frame r_x = %f, r_y = %f, r_z = %f at time %f\n", r_x, r_y, r_z, predict_time);
plot3(r_x, r_y, r_z, '-*r');
xlabel('x');
ylabel('y');
zlabel('z');
plot3(x, y, z, '--xb');
xlabel('x');
ylabel('y');
zlabel('z');
%}

%{
*****************************************************************************
% user defined functions
*****************************************************************************
%}

% Checksum (Modulo 10)
% Letters, blanks, periods, plus signs = 0; minus signs = 1
function [result] = chksum(str)
  result = false; c = 0;
  
  for k = 1:68
    if str(k) > '0' && str(k) <= '9'
      c = c + str(k) - 48;
    elseif str(k) == '-'
      c = c + 1;
    end
  end
  if mod(c,10) == str(69) - 48
    result = true;
  end
end

% Julian date
function [jd] = julian_date(y1, m1, d1, utc)
    
    if m1>2
        y = y1;
        m = m1;
    else
        y = y1-1;
        m = m1+12;
    end
    d = d1;
    h = utc/24;
    if y1<=1582 && m1<=10 && d1 <= 4
        % Julian calendar
        b = 0;
    elseif y1 == 1582 && m1 == 10 && d1 > 4 && d1 < 15
        % Gregorian calendar reform: 10 days (5 to 14 October 1582) were skipped.
        % In 1582 after 4 October follows the 15 October.
        d = 15;
        b = -10;
    else
        a = fix(y/100);
        b = 2 - a + fix(a/4);
    end
    jd = fix(365.25*(y+4716)) + fix(30.6001*(m+1)) + d + h + b - 1525.5;
end

% Calculate sidereal time using julian date
function [st] = siderial_time(y1, m1, d1, utc, longi)
    longi = 0;
    jd = julian_date(y1, m1, d1, utc);
    t = (jd - 2451545.0)/36525;
    % Greenwich siderial time at 0h UTC (hours)
    st = (24110.54841 + 8640184.812866 * t + 0.093104 * t^2 - 0.0000062 * t^3) / 3600;
    st = st + 1.00273790935*utc;
    % Local siderial time at given UTC (longitude in degrees)
    st = st + longi/15;
    st = mod(st, 24);
end

% canculate rotation matrix to transform reference frame ECEF to ECI
function dcm = angle(y1, m1, d1, utc)
    omega_earth = 0.261799387799149;                    % rad/hour
    
    jd = julian_date(y1, m1, d1, 0);
    jd_utc = julian_date(y1, m1, d1, utc);
    gmst = siderial_time(y1, m1, d1, utc, 0);
    offset_time = 0;                                    % second
    ecef_angle = omega_earth*(gmst-offset_time);      % rad
    ecef_angle = double(ecef_angle);

    C = cos(ecef_angle);
    S = sin(ecef_angle);

    % ECI in term of ECEF
    dcm = [C -S 0;
            S C 0;
            0 0 1];
end