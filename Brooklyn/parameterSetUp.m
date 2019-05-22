function [Type, Area_size, Zone_size, Nev, N, MaxCL, T, tdeta, v, cr, range, tolerance, max_D, ws, theta, dist_h] = parameterSetUp(StationCoordinates, ArrivalList)

days = 30;
Type = 1; % 1 is station based, 2 is free floating
Area_size = 20; %km
Zone_size = 2; %km
Nev = 100; %number of EVs
Nstation = size(StationCoordinates, 1); % number of stations
MaxCL = 5; % 5 %number of charge levels of an EV
tdeta = 1;  % 1 min
T = 60 * 24 * days / tdeta; %simulation horizon in minutes
v = 1; %free flow speed of EV is 1km/min = 60km/h
ChargeTime = 100; % time needed to charge from 0 to MaxCL
cr = MaxCL / ChargeTime; %time to get fully charged is 100 min
range = 200; % maximum km that can be covered by full battery
tolerance = 30 / tdeta; % the maximum time a customer is willing to wait
max_D = 3; %km, the maximum distance that a customer is willing to walk
ws = 5; %walking speed = 5km/h
theta = 0.02;
dist_h = ChargeTime / MaxCL;


% ave_using_time = mean(ArrivalList(:,5));
% ave_driving_dis = mean(ArrivalList(:,4));
% ave_arrival_rate = size(ArrivalList, 1) / T;
% Nev = round((ave_using_time + ave_driving_dis / range * 100) * ave_arrival_rate / 2);
% e = 0; % for non-EV


if Type == 1
    N = Nstation; %number of Stations
else
    N = (Area_size / Zone_size) ^ 2; %number of Zones
end