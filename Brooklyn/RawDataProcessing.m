load ZoneList
DrivingDistanceMatrixRaw = csvread('DrivingDistanceMatrix.csv');
DrivingDistanceMatrix = DrivingDistanceMatrixRaw(ZoneList, ZoneList);
DrivingDistanceMatrix = DrivingDistanceMatrix * 0.0003048;
save ('DrivingDistanceMatrix.mat', 'DrivingDistanceMatrix')

DrivingTimeMatrixRaw = csvread('DrivingTimeMatrix.csv');
DrivingTimeMatrix = DrivingTimeMatrixRaw(ZoneList, ZoneList);
DrivingTimeMatrix = round(DrivingTimeMatrix / 60);
save ('DrivingTimeMatrix.mat', 'DrivingTimeMatrix')

StationCoordinatesRaw = csvread('StationCoordinates.csv');
StationCoordinates = StationCoordinatesRaw(ZoneList, :);
save ('StationCoordinates.mat', 'StationCoordinates')