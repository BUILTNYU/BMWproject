for t = 1 : T
    EVs_t_temp = EVs_t{t, 1};
    if ~isempty(find(EVs_t_temp(:, 4) <= 0))
        break
    end
end

for t = 1 : T
    EVs_t_temp = EVs_t{t, 1};
    if EVs_t_temp(48, 8) == 0
        break
    end
end

%debug the customers with MaxCL nerver be served
Customers_Unserved_temp = Customers_Unserved{T + 1, 1};
CT_index = find(Customers_Unserved_temp(:, 8) == 5 & isnan(Customers_Unserved_temp(:, 15)));
CT = Customers_Unserved_temp(CT_index, :);
ZonesWithChargers_index = find(Capacity > 0);
ct_No = [];
zone = [];
charger = [];
dis2NearestChargingStation = [];
for i = 1 : size(CT, 1)
    ct_No = [ct_No; CT(i, 1)];
    zone = [zone; Nodes(CT(i, 3), 2)];
    zone_index = find(Zones(:, 1) == Nodes(CT(i, 3), 2));
    charger = [charger; Capacity(zone_index)];
    dis = DrivingDistanceMatrix(zone_index, ZonesWithChargers_index);
    dis_min = min(dis);
    dis2NearestChargingStation = [dis2NearestChargingStation; dis_min];
end
CT_check = [ct_No, zone, charger, dis2NearestChargingStation];

%check whether the customers that have not reservedany ev really have no available ev to reserve?
Customers_Unserved_temp = Customers_Unserved{T + 1, 1};
EVs_t_temp = EVs_t{T + 1, 1};
CT_index = find(Customers_Unserved_temp(:, 13) == 0);
CT = Customers_Unserved_temp(CT_index, :);
ZonesWithChargers_index = find(Capacity > 0);
ct_No = [];
zone = [];
availableEV = [];
for i = 1 : size(CT, 1)
    ct_No = [ct_No; CT(i, 1)];
    zone = [zone; Nodes(CT(i, 3), 2)];
    zone_index = find(Zones(:, 1) == Nodes(CT(i, 3), 2));
    availableEV = [availableEV; length(find(EVs_t_temp(:, 8) == 1 & EVs_t_temp(:, 2) == Zones(zone_index, 2) & EVs_t_temp(:, 3) == Zones(zone_index, 3)))];
end
CT_check = [ct_No, zone, charger, dis2NearestChargingStation, availableEV];




















