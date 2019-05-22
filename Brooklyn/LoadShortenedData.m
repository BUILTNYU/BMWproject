function [StationCoordinates, Capacity, DrivingDistanceMatrix, DrivingTimeMatrix, ArrivalList] = LoadShortenedData()
load StationCoordinates_s;
load Capacity_s;
load DrivingDistanceMatrix_s;
load DrivingTimeMatrix_s;
load ArrivalList_s;
StationCoordinates = StationCoordinates_s;
Capacity = Capacity_s;
DrivingDistanceMatrix = DrivingDistanceMatrix_s;
DrivingTimeMatrix = DrivingTimeMatrix_s;
ArrivalList = ArrivalList_s;