function [Total_Number_of_Customers, Customers_Unserved, Customers_Served, Customers_Finished, Customers_Lost, RebalanceList, CT_Node_t, betaMatrix, bMatrix, CT_Number_inS, TimeSeries, EV_statistic, EVs_t, EV_Node_t, initial_loc_t, station_capacity_t, station_waiting_number_t, RebalanceResult_t] = initialization(N, Nev, MaxCL, Zones)

Total_Number_of_Customers = 0; %initial total number of customers

Customers_Unserved = {[]};
Customers_Served = {[]};
Customers_Finished = {[]};
Customers_Lost = [];
RebalanceList = [];

CT_Node_t = {[]};
betaMatrix = {[]};
bMatrix = {[]};
CT_Number_inS = [];
TimeSeries = [];
EV_statistic = [];
initial_loc_t = {[]};
station_capacity_t = {[]};
station_waiting_number_t = {[]};
RebalanceResult_t = {[]};

% initialize the EVs_t
% we have Nev EVs in total in the system
No = (1 : Nev)';
Xcoor = [];
Ycoor = [];
for i = 1 : Nev
    r = rem(i, N) + 1;
    Xcoor = [Xcoor; Zones(r, 2)];
    Ycoor = [Ycoor; Zones(r, 3)];
end
CL = MaxCL * ones(Nev, 1);
C_Node = NaN(Nev, 1);
D_x = NaN(Nev, 1);
D_y = NaN(Nev, 1);
in_EV_Node = ones(Nev, 1);
InCharge = zeros(Nev, 1);
Occu = zeros(Nev, 1);
Customer_No = NaN(Nev, 1);
idleTime = zeros(Nev, 1);
occupiedKm = zeros(Nev, 1);
culmulatedCustomer = zeros(Nev, 1);
StartStation = NaN(Nev, 1);
ChargingStation = NaN(Nev, 1);
EndStation = NaN(Nev, 1);
DesiredCL = NaN(Nev, 1);
ChargingStationArrivalTime = NaN(Nev, 1);
Speed = NaN(Nev, 1);
RebalanceStartTime = NaN(Nev, 1);
RebalanceEndTime = NaN(Nev, 1);
RebalanceDis = NaN(Nev, 1);
ChargeWaitingTime = NaN(Nev, 1);
TimeToNextDestination = NaN(Nev, 1);
Pre_Xcoor = NaN(Nev, 1);
Pre_Ycoor = NaN(Nev, 1);
EVs_t = {[No, Xcoor, Ycoor, CL, C_Node, D_x, D_y, in_EV_Node, InCharge, Occu, Customer_No, idleTime, occupiedKm, culmulatedCustomer, StartStation, ChargingStation, EndStation, DesiredCL, ChargingStationArrivalTime, Speed, RebalanceStartTime, RebalanceEndTime, RebalanceDis, ChargeWaitingTime, TimeToNextDestination, Pre_Xcoor, Pre_Ycoor]};

% initialize the EV_Node_t
Number_of_available_EV = length(find(EVs_t{1, 1}(:, 8) == 1));
if Number_of_available_EV ~= 0
    temp = find(EVs_t{1, 1}(:, 8) == 1);
    Xcoor = [];
    Ycoor = [];
    CL = [];
    for i = 1 : Number_of_available_EV
        Xcoor = [Xcoor; EVs_t{1, 1}(temp(i), 2)];
        Ycoor = [Ycoor; EVs_t{1, 1}(temp(i), 3)];
        CL = [CL; ceil(EVs_t{1, 1}(temp(i), 4))]; %e.g. an ev with CL 1.1 belongs to EV_Node of CL 2
    end
    EV_Node_temp = [Xcoor, Ycoor, CL];
    unique_test = unique(EV_Node_temp, 'rows'); %test whether we have more than 1 EV in the same EV Node
    if size(EV_Node_temp,1) == size(unique_test,1)
        No = (1 : Number_of_available_EV)';
        Number_of_EVs = ones(Number_of_available_EV, 1);
        EV_Node_t = {[No, Xcoor, Ycoor, CL, Number_of_EVs]}; 
    else
        No = (1: size(unique_test, 1))';
        Xcoor = unique_test(:, 1);
        Ycoor = unique_test(:, 2);
        CL = unique_test(:, 3);
        Number_of_EVs = [];
        for i = 1 : size(unique_test, 1)
            Number_of_EVs = [Number_of_EVs; length(find(EV_Node_temp(:, 1) == unique_test(i, 1) & EV_Node_temp(:, 2) == unique_test(i, 2) & EV_Node_temp(:, 3) == unique_test(i, 3)))];
        end
        EV_Node_t = {[No, Xcoor, Ycoor, CL, Number_of_EVs]}; 
    end
else
    EV_Node_t = {[]};
end







