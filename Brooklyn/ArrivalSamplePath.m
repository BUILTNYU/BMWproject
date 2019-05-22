function [New_Arrivals_with_OD_t, Lamda_t, mu_t] = ArrivalSamplePath(ArrivalList, T, MaxCL, range, Nodes, Zones, minReturnCL)

arrival_temp = [];
for t = 1 : T
    New_Arrivals_with_OD_temp = [];
    temp = find(ArrivalList(:, 1) == t);
    if ~isempty(temp)
        for i = 1 : length(temp)
            Origin_zone = ArrivalList(temp(i), 2);
            Dest_zone = ArrivalList(temp(i), 3);
            Origin_zone_index = find(Zones(:, 1) == Origin_zone);
            Dest_zone_index = find(Zones(:, 1) == Dest_zone);
            DrivingDis = ArrivalList(temp(i), 4);
            CL_temp = min(ceil(DrivingDis / range * MaxCL + minReturnCL(Dest_zone_index)), MaxCL); %if the drivingdis > range, then EV will get charged while driving, we make such customers have CL of MaxCL
            % + minReturnCL(Dest_zone_index) because we want to make sure any EV to be returned to its destination with enough CL
            Origin_node = Nodes(find(Nodes(:, 2) == Origin_zone & Nodes(:, 3) == CL_temp), 1);
            O_x = Zones(Origin_zone_index, 2);
            O_y = Zones(Origin_zone_index, 3);
            D_x = Zones(Dest_zone_index, 2);
            D_y = Zones(Dest_zone_index, 3);
            BookingTime = ArrivalList(temp(i), 5);
            New_Arrivals_with_OD_temp = [New_Arrivals_with_OD_temp; Origin_node, O_x, O_y, D_x, D_y, Dest_zone, DrivingDis, BookingTime];
            arrival_temp = [arrival_temp; t, Origin_node];
        end
    end
    New_Arrivals_with_OD_t{t, 1} = New_Arrivals_with_OD_temp;
end

for t = 1 : T
    for i = 1 : size(Nodes, 1)
%         %fixed lamda
%         Lamda_temp(i) = length(find(arrival_temp(:, 1) > t & arrival_temp(:, 1) <= (t + 60) & arrival_temp(:, 2) == i));
%         mu = sum(Lamda_temp) / 1.5;
        Lamda_temp(i) = length(find(arrival_temp(:, 2) == i)) / (ArrivalList(end, 1) / 60);
        mu = 100000000000;
    end
    Lamda_t{t, 1} = Lamda_temp';
    mu_t{t, 1} = mu * ones(size(Nodes, 1), 1);
end




