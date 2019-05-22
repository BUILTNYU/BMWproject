function [Zones, Nodes, minReturnCL, TTMatrix, CSMatrix, adj, dist_tij, charg_stations] = geometrySetUp(StationCoordinates, DrivingTimeMatrix, DrivingDistanceMatrix, Area_size, Zone_size, N, Type, Capacity, range, MaxCL, cr)

if Type == 1
    Zones = [StationCoordinates, Capacity];
else
    No = (1 : (round(Area_size / Zone_size) * round(Area_size / Zone_size)))';
    Xcoor = (Zone_size / 2): Zone_size: (Area_size - Zone_size / 2); %Xcoor of each Zone center (station)
    Ycoor = (Zone_size / 2): Zone_size: (Area_size - Zone_size / 2); %Ycoor of each Zone center
    Xcoor = repmat(Xcoor', round(Area_size / Zone_size), 1);
    Ycoor = kron(Ycoor', ones(round(Area_size / Zone_size), 1));
    Total_chargers = randi([0,5], N, 1);
    Zones = [No, Xcoor, Ycoor, Total_chargers]; %Station / Zone list
end

No = (1 : N * MaxCL)';
CL = kron((1 : MaxCL)',ones(N, 1));
Zone_No = Zones(:, 1);
Zone_No = repmat(Zone_No, MaxCL, 1);
Nodes = [No, Zone_No, CL];

minReturnCL = []; % the min CL that an EV should have when it is returned to a station, this CL should make sure that the EV could be rebalanced to the farest charging station
ChargingStations = find(Capacity(:) > 0);
for i = 1 : N
    farest_ChargingStation = ChargingStations(find(DrivingDistanceMatrix(i, ChargingStations) == max(DrivingDistanceMatrix(i, ChargingStations)), 1));
    min_CL = MaxCL * DrivingDistanceMatrix(i, farest_ChargingStation) / range;
    minReturnCL = [minReturnCL; min_CL];
end

TTMatrix = zeros(N * MaxCL, N * MaxCL); % Travel Time matrix when rebalancing
CSMatrix = NaN(N * MaxCL, N * MaxCL); % Charging Station matrix
for i = 1 : size(Nodes, 1)
    for j = 1 : size(Nodes, 1)
        Origin_ZoneNo_index = find(Zones(:, 1) == Nodes(i, 2));
        Dest_ZoneNo_index = find(Zones(:, 1) == Nodes(j, 2));
        if Nodes(i, 3) > Nodes(j, 3) % if the origin CL >= dest CL, no need to charge during rebalance
            TTMatrix(i, j) = DrivingTimeMatrix(Origin_ZoneNo_index, Dest_ZoneNo_index);
        else % if charging is needed during rebalance
            if Zones(Origin_ZoneNo_index, 4) > 0 % if origin is a charging station
                TTMatrix(i, j) = DrivingTimeMatrix(Origin_ZoneNo_index, Dest_ZoneNo_index) + (Nodes(j, 3) - Nodes(i, 3)) / cr;
                CSMatrix(i, j) = Nodes(i, 2);
            elseif Zones(Dest_ZoneNo_index, 4) > 0 % if dest is a charging station
                TTMatrix(i, j) = DrivingTimeMatrix(Origin_ZoneNo_index, Dest_ZoneNo_index) + (Nodes(j, 3) - Nodes(i, 3)) / cr;
                CSMatrix(i, j) = Nodes(j, 2);
            else % if neither origin nor dest is a charging station
                cs = find(Zones(:, 4) > 0);
                dis = [];
                for k = 1 : length(cs)
                    dis = [dis; cs(k), (DrivingTimeMatrix(Origin_ZoneNo_index, cs(k)) + DrivingTimeMatrix(cs(k), Dest_ZoneNo_index))];
                end
                cs_min = find(dis(:, 2) == min(dis(:, 2)), 1);
                TTMatrix(i, j) = dis(cs_min, 2) + (Nodes(j, 3) - Nodes(i, 3)) / cr;
                CSMatrix(i, j) = Zones(dis(cs_min, 1), 1);
            end
        end
    end
end

adj = zeros(size(Nodes, 1), size(Nodes, 1));
dist_tij = zeros(size(Nodes, 1), size(Nodes, 1));
for i = 1 : size(Nodes, 1)
    for j = 1 : size(Nodes, 1)
        if Nodes(i, 3) == Nodes(j, 3) % adj = 1 if in the same layer
            adj(i, j) = 1;
        elseif abs(Nodes(j, 3) - Nodes(i, 3)) == 1 && abs(Nodes(i, 1) - Nodes(j, 1)) == N && Zones(find(Zones(:, 1) == Nodes(i, 2)), 4) > 0 % adj= 1 if  i is one layer lower than j, or one layer higher than j, i and j are in the same station, charging station
            adj(i, j) = 1;
        end
        dist_tij(i, j) = DrivingDistanceMatrix(find(Zones(:, 1) == Nodes(i, 2)), find(Zones(:, 1) == Nodes(j, 2)));
    end
end

charg_stations = find(Capacity > 0);

