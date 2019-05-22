function [StationCoordinates, Capacity, DrivingDistanceMatrix, DrivingTimeMatrix, ArrivalList] = LoadRealData()
load StationCoordinates;
load Capacity;
load DrivingDistanceMatrix;
load DrivingTimeMatrix;
StationCoordinates = StationCoordinates;
Capacity = Capacity;
DrivingDistanceMatrix = DrivingDistanceMatrix;
DrivingTimeMatrix = DrivingTimeMatrix;

% load Demand_ST14400_EM10;
% ArrivalList = Dem;

load ArrivalList
ArrivalList = ArrivalList;


% delete the data of the stations with 0 distance and travel time to all other stations (those stations have centroids in the sea)
row_all_zeros = find(all(DrivingDistanceMatrix == 0,2));
col_all_zeros = find(all(DrivingDistanceMatrix == 0,1));
station_no_index = StationCoordinates(row_all_zeros, 1);
DrivingDistanceMatrix(row_all_zeros, :) = [];
DrivingDistanceMatrix(:, col_all_zeros) = [];
DrivingTimeMatrix(row_all_zeros, :) = [];
DrivingTimeMatrix(:, col_all_zeros) = [];

StationCoordinates(row_all_zeros, :) = [];
Capacity(row_all_zeros, :) = [];

delete_list = [];
for i = 1 : size(ArrivalList, 1)
    for j = 1 : length(row_all_zeros)
        if ArrivalList(i, 2) == station_no_index(j)
            delete_list = [delete_list, i];
        elseif ArrivalList(i, 3) == station_no_index(j)
            delete_list = [delete_list, i];
        end
    end
end
ArrivalList(delete_list, :) = [];