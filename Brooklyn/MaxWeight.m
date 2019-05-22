function [EVs_t, Customers_Unserved] = MaxWeight(MaxCL, Nev, v, t, Zones, Nodes, TTMatrix, CSMatrix, DrivingTimeMatrix, Customers_Unserved, EVs_t)
EVs_t_temp_pre = EVs_t{t, 1};
EVs_t_temp = EVs_t{t + 1, 1};
Customers_Unserved_temp = Customers_Unserved{t + 1, 1};

if size(Customers_Unserved_temp, 1) >= 1
    unreserved_customers = Customers_Unserved_temp(find(Customers_Unserved_temp(:, 13) == 0  & Customers_Unserved_temp(:, 14) == 0), :);
    number_of_unreserved_customers = size(unreserved_customers, 1);
else
    number_of_unreserved_customers = 0;
end
if number_of_unreserved_customers > 0
    % for each available EV, find the customer that is nearest to it, then
    % mark the EV with this customer
    for i = 1 : Nev
        if EVs_t_temp_pre(i, 8) == 1 && EVs_t_temp(i, 8) == 1 % if the ev is available
            ev_zone = Zones(find(Zones(:, 2) == EVs_t_temp(i, 2) & Zones(:, 3) == EVs_t_temp(i, 3)), 1);
            ev_node = Nodes(find(Nodes(:, 2) == ev_zone & Nodes(:, 3) == ceil(EVs_t_temp(i, 4))), 1);
            dis = [];
            for j = 1 : number_of_unreserved_customers
                dis = [dis; unreserved_customers(j, 1), TTMatrix(ev_node, unreserved_customers(j, 3))];
            end
            customerNo_min_dis = dis(find(dis(:, 2) == min(dis(:, 2)), 1), 1); % the customer No. that is nearest to this EV
            EVs_t_temp(i, 11) = customerNo_min_dis; % mark the EV with this customer, note that multiple EVs may be marked by the same customer
            %customer_min_dis_CTnode = unreserved_customers(unreserved_customers(:, 1) == customerNo_min_dis, 3);
        end
    end
    for i = 1 : Nev
        Cust_No = EVs_t_temp(i, 11);
        Cust_No_index = find(unreserved_customers(:, 1) == Cust_No);
        if EVs_t_temp(i, 8) == 1 && ~isempty(Cust_No_index) % if the vehicle's customer could be found in the unserved customer list, means this customer has not been researved by other vehicles
            ev_temp = EVs_t_temp(find(EVs_t_temp(:, 8) == 1 & EVs_t_temp(:, 11) == EVs_t_temp(i, 11)), :);
            if size(ev_temp, 1) > 1 % if one customer is marked by more than one 1 EV
                dis = [];
                for j = 1 : size(ev_temp, 1)
                    ev_zone = Zones(find(Zones(:, 2) == ev_temp(j, 2) & Zones(:, 3) == ev_temp(j, 3)), 1);
                    ev_node = Nodes(find(Nodes(:, 2) == ev_zone & Nodes(:, 3) == ceil(ev_temp(j, 4))), 1);
                    dis = [dis; ev_temp(j, 1), TTMatrix(ev_node, unreserved_customers(Cust_No_index, 3))];
                end
                evNo_min_dis = dis(find(dis(:, 2) == min(dis(:, 2)), 1), 1); % find the EV No. with the minimum distance to the customer
            else
                evNo_min_dis = EVs_t_temp(i, 1);
            end
            if EVs_t_temp(i, 1) == evNo_min_dis % if this ev is the nearest ev
                ct_Node = unreserved_customers(Cust_No_index, 3);
                ct_No_index = find(Customers_Unserved_temp(:, 1) == Cust_No);
                Customers_Unserved_temp(ct_No_index, 14) = 1;
                EVs_t_temp(i, 15) = Zones(find(Zones(:, 2) == EVs_t_temp(i, 2) & Zones(:, 3) == EVs_t_temp(i, 3)), 1);
                EVs_t_temp(i, 17) = Nodes(ct_Node, 2);
                initialCL = EVs_t_temp(i, 4);
                desiredCL = min(Nodes(ct_Node, 3) + 1, MaxCL);
                EVs_t_temp(i, 5) = Nodes(find(Nodes(:, 2) == EVs_t_temp(i, 17) & Nodes(:, 3) == desiredCL), 1);
                EVs_t_temp(i, 18) = desiredCL;
                EVs_t_temp(i, 20) = 1 * v; % assume free flow speed
                EVs_t_temp(i, 21) = t + 1;
                EVs_t_temp(i, 23) = 0;
                if initialCL < desiredCL % if ev gets charged during rebalance
                    EVs_t_temp(i, 16) = CSMatrix(find(Nodes(:, 2) == EVs_t_temp(i, 15) & Nodes(:, 3) == ceil(initialCL)), EVs_t_temp(i, 5));
                    zone_index = find(Zones(:, 1) == EVs_t_temp(i, 16));
                    EVs_t_temp(i, 6) = Zones(zone_index, 2);
                    EVs_t_temp(i, 7) = Zones(zone_index, 3);
                    if EVs_t_temp(i, 15) == EVs_t_temp(i, 16) % if the start station is the charging station
                        if EVs_t_temp(i, 9) == 0 && isnan(EVs_t_temp(i, 19))
                            EVs_t_temp(i, 19) = t + 1;
                        end
                    else  % if the start station is not the charging station
                        EVs_t_temp(i, 19) = NaN;
                    end
                else % if ev does not get charged during rebalance
                    EVs_t_temp(i, 9) = 0;
                    EVs_t_temp(i, 19) = NaN;
                    zone_index = find(Zones(:, 1) == EVs_t_temp(i, 17));
                    EVs_t_temp(i, 6) = Zones(zone_index, 2);
                    EVs_t_temp(i, 7) = Zones(zone_index, 3);
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
                unreserved_customers(Cust_No_index, :) = []; %the customer has already find his/her nearest ev and been reserved, delete the customer from unreserved customer list
            end
        end
    end
    for i = 1 : Nev
        if EVs_t_temp(i, 8) == 1 && ~isnan(EVs_t_temp(i, 11)) && isnan(EVs_t_temp(i, 5))
            EVs_t_temp(i, 11) = NaN;
        elseif EVs_t_temp(i, 8) == 1 && ~isnan(EVs_t_temp(i, 11)) && ~isnan(EVs_t_temp(i, 5))
            EVs_t_temp(i, 8) = 0;
        end
    end
end

EVs_t{t + 1, 1} = EVs_t_temp;
Customers_Unserved{t + 1, 1} = Customers_Unserved_temp;
