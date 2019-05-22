function [EVs_t, RebalanceResult_t] = generate_output(Nev, v, t, Zones, Nodes, EVs_t, RebalancePath, initial_loc_t, DrivingTimeMatrix, RebalanceResult_t)

EVs_t_temp = EVs_t{t + 1, 1};
initial_loc = initial_loc_t{t + 1, 1};

if size(RebalancePath,2) > 0
    for i = 1 : size(RebalancePath,2)
        RebalancePath{1, i} = flip(RebalancePath{1, i});
    end
end

for i = 1 : Nev
    if ismember(i, initial_loc(:, 1)) % if the ev is an input for MIP
        temp = find(initial_loc(:, 1) == i);
        if length(RebalancePath{1, temp}) ~= 1 % if the ev needs to be rebalanced
            EVs_t_temp(i, 5) = RebalancePath{1, temp}(end);
            EVs_t_temp(i, 8) = 0;
            EVs_t_temp(i, 15) = Nodes(RebalancePath{1, temp}(1), 2);
            EVs_t_temp(i, 17) = Nodes(RebalancePath{1, temp}(end), 2);
            initialCL = Nodes(RebalancePath{1, temp}(1), 3);
            desiredCL = Nodes(EVs_t_temp(i, 5), 3);
            EVs_t_temp(i, 18) = desiredCL;
            EVs_t_temp(i, 20) = 1 * v; % assume free flow speed
            EVs_t_temp(i, 21) = t + 1;
            EVs_t_temp(i, 23) = 0;
            if initialCL ~= desiredCL % if ev gets charged during rebalance
                CL = [];
                for j = 1 : length(RebalancePath{1, temp})
                    CL = [CL, Nodes(RebalancePath{1, temp}(j), 3)];
                end
                EVs_t_temp(i, 16) = Nodes(RebalancePath{1, temp}(find(CL == desiredCL, 1)), 2);
                zone_no_index = find(Zones(:, 1) == EVs_t_temp(i, 16));
                EVs_t_temp(i, 6) = Zones(zone_no_index, 2);
                EVs_t_temp(i, 7) = Zones(zone_no_index, 3);
                if EVs_t_temp(i, 15) ~= EVs_t_temp(i, 16) % if the start station is not the charging station
                    EVs_t_temp(i, 9) = 0;
                    EVs_t_temp(i, 19) = NaN;
                end
            else % if ev does not get charged during rebalance
                EVs_t_temp(i, 9) = 0;
                EVs_t_temp(i, 19) = NaN;
                zone_no_index = find(Zones(:, 1) == EVs_t_temp(i, 17));
                EVs_t_temp(i, 6) = Zones(zone_no_index, 2);
                EVs_t_temp(i, 7) = Zones(zone_no_index, 3);
            end
            % update the ev's time to next destination
            if ~isnan(EVs_t_temp(i, 16)) && EVs_t_temp(i, 15) ~= EVs_t_temp(i, 16) % if the ev need to get charged, but the start station is not the charging station
                row_temp = find(Zones(:, 1) == EVs_t_temp(i, 15));
                column_temp = find(Zones(:, 1) == EVs_t_temp(i, 16));
                EVs_t_temp(i, 25) = t + 1 + DrivingTimeMatrix(row_temp, column_temp);
                EVs_t_temp(i, 26) = EVs_t_temp(i, 2);
                EVs_t_temp(i, 27) = EVs_t_temp(i, 3);
                EVs_t_temp(i, 2) = NaN;
                EVs_t_temp(i, 3) = NaN;
            elseif isnan(EVs_t_temp(i, 16)) && EVs_t_temp(i, 15) ~= EVs_t_temp(i, 17) % if the ev does not need to get charged and the start station is not the end station
                row_temp = find(Zones(:, 1) == EVs_t_temp(i, 15));
                column_temp = find(Zones(:, 1) == EVs_t_temp(i, 17));
                EVs_t_temp(i, 25) = t + 1 + DrivingTimeMatrix(row_temp, column_temp);
                EVs_t_temp(i, 26) = EVs_t_temp(i, 2);
                EVs_t_temp(i, 27) = EVs_t_temp(i, 3);
                EVs_t_temp(i, 2) = NaN;
                EVs_t_temp(i, 3) = NaN;
            end
        end
    end
end

EVs_t{t + 1, 1} = EVs_t_temp;
RebalanceResult_t{t + 1, 1} = RebalancePath;