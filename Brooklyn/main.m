
[StationCoordinates, Capacity, DrivingDistanceMatrix, DrivingTimeMatrix, ArrivalList] = LoadRealData();
% [StationCoordinates, Capacity, DrivingDistanceMatrix, DrivingTimeMatrix, ArrivalList] = LoadShortenedData();
[Type, Area_size, Zone_size, Nev, N, MaxCL, T, tdeta, v, cr, range, tolerance, max_D, ws, theta, dist_h] = parameterSetUp(StationCoordinates, ArrivalList);
[Zones, Nodes, minReturnCL, TTMatrix, CSMatrix, adj, dist_tij, charg_stations] = geometrySetUp(StationCoordinates, DrivingTimeMatrix, DrivingDistanceMatrix, Area_size, Zone_size, N, Type, Capacity, range, MaxCL, cr);
[New_Arrivals_with_OD_t, Lamda_t, mu_t] = ArrivalSamplePath(ArrivalList, T, MaxCL, range, Nodes, Zones, minReturnCL);
[Total_Number_of_Customers, Customers_Unserved, Customers_Served, Customers_Finished, Customers_Lost, RebalanceList, CT_Node_t, betaMatrix, bMatrix, CT_Number_inN, TimeSeries, EV_statistic, EVs_t, EV_Node_t, initial_loc_t, station_capacity_t, station_waiting_number_t, RebalanceResult_t] = initialization(N, Nev, MaxCL, Zones);
progress = 0;
for t = 1 : T
    [Total_Number_of_Customers, Customers_Unserved] = update_customers_unserved_based_on_newArrivals(t, Nodes, ws, Total_Number_of_Customers, New_Arrivals_with_OD_t, Customers_Unserved);
     
    %%% --- MIP heuristic --- %%%
    [EVs_t, Customers_Unserved] = update_EVs(Nev, MaxCL, tdeta, max_D, cr, range, t, Zones, Nodes, minReturnCL, DrivingDistanceMatrix, DrivingTimeMatrix, Customers_Served, Customers_Unserved, EVs_t);
    if rem(t, 60) == 0  % the rebalance is optimized every 60 min
        [initial_loc_t, station_capacity_t, station_waiting_number_t] = generate_input(Nev, N, t, Zones, Nodes, EVs_t, initial_loc_t, station_capacity_t, station_waiting_number_t);
        if size(initial_loc_t{t + 1, 1}, 1) ~= 0 && max(station_capacity_t{t + 1, 1}(:, 2)) >= 1
            [RebalanceDecision,RebalancePath] = heuristic_mrk(Lamda_t{t, 1}, N, MaxCL, adj, dist_tij, theta, dist_h, station_capacity_t{t + 1, 1}(:, 2)', charg_stations', sort(initial_loc_t{t + 1, 1}(:, 2)'), mu_t{t, 1});
            [EVs_t, RebalanceResult_t] = generate_output(Nev, v, t, Zones, Nodes, EVs_t, RebalancePath, initial_loc_t, DrivingTimeMatrix, RebalanceResult_t);
        else
            RebalanceResult_t{t + 1, 1} = [];
        end
    end 
         
    %%% --- NoRebalance --- %%%
%     [EVs_t, Customers_Unserved] = update_EVs(Nev, MaxCL, tdeta, max_D, cr, range, t, Zones, Nodes, minReturnCL, DrivingDistanceMatrix, DrivingTimeMatrix, Customers_Served, Customers_Unserved, EVs_t);
    
    %%% --- MaxWeight --- %%%
%     [EVs_t, Customers_Unserved] = update_EVs(Nev, MaxCL, tdeta, max_D, cr, range, t, Zones, Nodes, minReturnCL, DrivingDistanceMatrix, DrivingTimeMatrix, Customers_Served, Customers_Unserved, EVs_t);
%     [EVs_t, Customers_Unserved] = MaxWeight(MaxCL, Nev, v, t, Zones, Nodes, TTMatrix, CSMatrix, DrivingTimeMatrix, Customers_Unserved, EVs_t);
    
    [EVs_t, RebalanceList] = RebalanceInfoCollection(Nev, t, EVs_t, RebalanceList);
    [Customers_Unserved, Customers_Served, Customers_Finished, Customers_Lost] = update_customers_unserved_served_finished(Nev, t, tolerance, EVs_t, Customers_Unserved, Customers_Served, Customers_Finished, Customers_Lost);
    [CT_Number_inN, TimeSeries, EV_statistic] = update_evaluations(N, MaxCL, t, Customers_Unserved, Customers_Served, Customers_Finished, Customers_Lost, EVs_t, CT_Number_inN, TimeSeries, EV_statistic, RebalanceList);
    if t / T >= progress
        [num2str(round(t / T, 3) * 100), '%']
        progress = progress + 0.1;
    end
end

[CN_ave, WT_ave, LostCN, sum_idleKm, sum_occupiedKm] = output(t-1, CT_Number_inN, Customers_Served, Customers_Finished, Customers_Lost, EVs_t);
