function [EVs_t, Customers_Unserved] = update_EVs(Nev, MaxCL, tdeta, max_D, cr, range, t, Zones, Nodes, minReturnCL, DrivingDistanceMatrix, DrivingTimeMatrix, Customers_Served, Customers_Unserved, EVs_t)
EVs_t_temp = EVs_t{t, 1};
EVs_t_temp_pre = EVs_t{t, 1};
Customers_Unserved_temp = Customers_Unserved{t, 1};
Customers_Served_temp = Customers_Served{t, 1};

% C_Node is the rebalance destination, 
% when an EV is arranged for rebalance, C_Node becomes non-NaN
% Once an EV is reserved by a customer, Customer_No becomes non-NaN
% Once the EV is picked up by a customer, Occu becomes 1

%update the info of EVs
for i = 1 : Nev
    % update CL
    if EVs_t_temp(i, 9) == 1 % if the ev is being charged
        EVs_t_temp(i, 4) = min([EVs_t_temp(i, 4) + cr * tdeta, MaxCL]);
        if EVs_t_temp(i, 4) == MaxCL % if the ev is fully charged
            EVs_t_temp(i, 9) = 0;
            EVs_t_temp(i, 19) = NaN;
        end
    else % if the ev is not being charged
        if ~isnan(EVs_t_temp(i, 19)) % if the ev is still waiting to get charged
            zone_temp = find(Zones(:, 2) == EVs_t_temp(i, 2) & Zones(:, 3) == EVs_t_temp(i, 3));
            incharge_temp = length(find(EVs_t_temp(:, 2) == Zones(zone_temp, 2) & EVs_t_temp(:, 3) == Zones(zone_temp, 3) & EVs_t_temp(:, 9) == 1));
            if incharge_temp < Zones(zone_temp, 4) % if there is available charger left
                available_charger = Zones(zone_temp, 4) - incharge_temp;
                waiting_ev_No = find(EVs_t_temp(:, 2) == Zones(zone_temp, 2) & EVs_t_temp(:, 3) == Zones(zone_temp, 3) & EVs_t_temp(:, 9) == 0 & ~isnan(EVs_t_temp(:, 19)));
                waiting_ev_list = EVs_t_temp(waiting_ev_No, :);
                waiting_ev_list_sort = sortrows(waiting_ev_list, 19);
                waiting_order = find(waiting_ev_list_sort(:, 1) == i);
                if waiting_order <= available_charger % if the ev is in the enough front of the queue
                    EVs_t_temp(i, 9) = 1; % charge the ev
                    EVs_t_temp(i, 24) = t + 1 - EVs_t_temp(i, 19);
                    EVs_t_temp(i, 19) = NaN;
                end
            end
        end
    end
    if ~isnan(EVs_t_temp(i, 5)) % if this is an ev being rebalanced
        if EVs_t_temp_pre(i, 25) <= t + 1 % if the ev is about to arrive its next destination, we have here <= but not == because we have data booking time or driving time equal to 0
            row_index = find(Zones(:, 2) == EVs_t_temp(i, 26) & Zones(:, 3) == EVs_t_temp(i, 27));
            column_index = find(Zones(:, 2) == EVs_t_temp(i, 6) & Zones(:, 3) == EVs_t_temp(i, 7));
            d_traveled = DrivingDistanceMatrix(row_index, column_index);
            % update the ev's location
            EVs_t_temp(i, 2) = EVs_t_temp(i, 6);
            EVs_t_temp(i, 3) = EVs_t_temp(i, 7);
            EVs_t_temp(i, 26) = NaN;
            EVs_t_temp(i, 27) = NaN;
            EVs_t_temp(i, 25) = NaN;
            % update the ev's charge level
            EVs_t_temp(i, 4) = max(EVs_t_temp(i, 4) - MaxCL * d_traveled / range, 0.001); % make sure CL will not be <= 0
            EVs_t_temp(i, 23) = EVs_t_temp(i, 23) + d_traveled;
        else
            d_traveled = 0;
        end
%         % update location
%         d = sqrt((EVs_t_temp_pre(i, 2) - EVs_t_temp(i, 6)) ^ 2 + (EVs_t_temp_pre(i, 3) - EVs_t_temp(i, 7)) ^ 2);
%         if d <= EVs_t_temp(i, 20) * tdeta
%             EVs_t_temp(i, 2) = EVs_t_temp(i, 6);
%             EVs_t_temp(i, 3) = EVs_t_temp(i, 7);
%         else
%             EVs_t_temp(i, 2) = EVs_t_temp_pre(i, 2) + (EVs_t_temp(i, 6) - EVs_t_temp_pre(i, 2)) * EVs_t_temp_pre(i, 20) * tdeta / d;
%             EVs_t_temp(i, 3) = EVs_t_temp_pre(i, 3) + (EVs_t_temp(i, 7) - EVs_t_temp_pre(i, 3)) * EVs_t_temp_pre(i, 20) * tdeta / d;
%         end
%         d_traveled = sqrt((EVs_t_temp_pre(i, 2) - EVs_t_temp(i, 2)) ^ 2 + (EVs_t_temp_pre(i, 3) - EVs_t_temp(i, 3)) ^ 2);
%         EVs_t_temp(i, 23) = EVs_t_temp(i, 23) + d_traveled;
%         EVs_t_temp(i, 4) = EVs_t_temp(i, 4) - MaxCL * d_traveled / range;

        zone_temp_index = find(Zones(:, 2) == EVs_t_temp(i, 2) & Zones(:, 3) == EVs_t_temp(i, 3));
        if ~isempty(zone_temp_index)
            zone_temp = Zones(zone_temp_index, 1);
        end
        if ~isempty(zone_temp_index) % if the ev is at a certain station, i.e. we know the EV's current exact location
            if zone_temp == EVs_t_temp(i, 16) % if the ev is at its charging station
                if d_traveled ~= 0 % if the ev just arrived at its charging station
                    EVs_t_temp(i, 19) = t + 1;
                end
                if EVs_t_temp(i, 4) >= EVs_t_temp(i, 18) % if the ev has reached its desired CL
                    EVs_t_temp(i, 9) = 0; % remove the ev from chargers
                    if EVs_t_temp(i, 16) == EVs_t_temp(i, 17) % if the charging station is also the end station
                        if EVs_t_temp(i, 4) < MaxCL && Zones(zone_temp_index, 4) > 0 % if the ev has not reached MaxCL and it is in a charging station
                            EVs_t_temp(i, 19) = t + 1;
                        end
                        EVs_t_temp(i, 5) = NaN;
                        EVs_t_temp(i, 8) = 1;
                        EVs_t_temp(i, 11) = NaN;
                        EVs_t_temp(i, 15) = NaN;
                        EVs_t_temp(i, 16) = NaN;
                        EVs_t_temp(i, 17) = NaN;
                        EVs_t_temp(i, 18) = NaN;
                        EVs_t_temp(i, 20) = NaN;
                        EVs_t_temp(i, 22) = t + 1;
                    else % if the charging station is not the end station
                        row_temp = find(Zones(:, 1) == EVs_t_temp(i, 16));
                        column_temp = find(Zones(:, 1) == EVs_t_temp(i, 17));
                        % update the ev's destination location
                        EVs_t_temp(i, 6) = Zones(column_temp, 2);
                        EVs_t_temp(i, 7) = Zones(column_temp, 3);
                        % update the ev's pre location, pre location is current location
                        EVs_t_temp(i, 26) = EVs_t_temp(i, 2);
                        EVs_t_temp(i, 27) = EVs_t_temp(i, 3);
                        % update the ev's current location, current location becomes unknown
                        EVs_t_temp(i, 2) = NaN;
                        EVs_t_temp(i, 3) = NaN;
                        % update the ev's time to next destination
                        EVs_t_temp(i, 25) = t + 1 + DrivingTimeMatrix(row_temp, column_temp);
                    end
                end
            % if ev (does not need to charge and is at end station but not charging station) or (is at end station but not charging station after charged)
            % we do not just say if ev is at its end station because an ev is possible to pass its end station when it is on its way to its charging station
%            elseif (zone_temp == EVs_t_temp(i, 17) && isnan(EVs_t_temp(i, 16))) || (zone_temp == EVs_t_temp(i, 17) && ~isnan(EVs_t_temp(i, 16)) && ~isnan(EVs_t_temp(i, 24)))
            elseif zone_temp == EVs_t_temp(i, 17) % if the ev is at its end station but not charging station
                if Zones(zone_temp_index, 4) > 0 && EVs_t_temp(i, 4) < MaxCL % zone_temp_index here is the index for EVs_t_temp(i, 17)
                    EVs_t_temp(i, 19) = t + 1;
                end
                EVs_t_temp(i, 5) = NaN;
                EVs_t_temp(i, 8) = 1;
                EVs_t_temp(i, 11) = NaN;
                EVs_t_temp(i, 15) = NaN;
                EVs_t_temp(i, 16) = NaN;
                EVs_t_temp(i, 17) = NaN;
                EVs_t_temp(i, 18) = NaN;
                EVs_t_temp(i, 20) = NaN;
                EVs_t_temp(i, 22) = t + 1;
            end
        end
        EVs_t_temp(i, 12) = EVs_t_temp(i, 12) + tdeta;     
    else % if this is an ev not being rebalanced
        % for non-rebalancing ev, we have available ev, ev reserved by
        % customers, ev already rented by customers (out of system)
        if ~isnan(EVs_t_temp(i, 11)) && EVs_t_temp(i, 10) == 0  && EVs_t_temp(i, 8) == 0 % if the ev is actively reserved by a customer
            temp = find(Customers_Unserved_temp(:, 1) == EVs_t_temp(i, 11));
            if Customers_Unserved_temp(temp, 15) == t % if the customer is about to arrive the reserved ev
                EVs_t_temp(i, 6) = Customers_Unserved_temp(temp, 6);
                EVs_t_temp(i, 7) = Customers_Unserved_temp(temp, 7);
                EVs_t_temp(i, 9) = 0;
                EVs_t_temp(i, 10) = 1;
                EVs_t_temp(i, 19) = NaN;
                % update the ev's pre location, pre location is current location
                EVs_t_temp(i, 26) = EVs_t_temp(i, 2);
                EVs_t_temp(i, 27) = EVs_t_temp(i, 3);
                % update the ev's current location, current location becomes unknown
                EVs_t_temp(i, 2) = NaN;
                EVs_t_temp(i, 3) = NaN;
                % update the ev's time to next destination
                EVs_t_temp(i, 25) = t + 1 + Customers_Unserved_temp(temp, 10);
            end
            EVs_t_temp(i, 12) = EVs_t_temp(i, 12) + tdeta;
        elseif EVs_t_temp(i, 10) == 1 % if the ev is being used by customers
            if EVs_t_temp(i, 25) <= t + 1 % if it's time for the customer to return the ev, we have here <= but not == because we have data booking time or driving time equal to 0
                temp = find(Customers_Served_temp(:, 1) == EVs_t_temp(i, 11));
                % update the ev's current location, current location is the destination location
                EVs_t_temp(i, 2) = EVs_t_temp(i, 6);
                EVs_t_temp(i, 3) = EVs_t_temp(i, 7);
                % update the ev's pre location, pre location is NaN
                EVs_t_temp(i, 26) = NaN;
                EVs_t_temp(i, 27) = NaN;
                % update the ev's time to next destination, it becomes NaN
                EVs_t_temp(i, 25) = NaN;
                dest_temp = find(Zones(:, 2) == EVs_t_temp(i, 6) & Zones(:, 3) == EVs_t_temp(i, 7));
                drivingdis_temp = Customers_Served_temp(temp, 9);
                if ceil(drivingdis_temp / range * MaxCL + minReturnCL(dest_temp)) <= MaxCL % if no charging during driving
                    EVs_t_temp(i, 4) = max(EVs_t_temp(i, 4) - MaxCL * drivingdis_temp / range, minReturnCL(dest_temp));
                else % if charging is needed during driving, the return CL is a random number >= minReturnCL, but <= MaxCL
                    EVs_t_temp(i, 4) = min(minReturnCL(dest_temp) + MaxCL * rand(), MaxCL);
                end
                EVs_t_temp(i, 8) = 1;
                EVs_t_temp(i, 10) = 0;
                EVs_t_temp(i, 11) = NaN;
                EVs_t_temp(i, 13) = EVs_t_temp(i, 13) + drivingdis_temp;
                EVs_t_temp(i, 14) = EVs_t_temp(i, 14) + 1;
                if Zones(dest_temp, 4) > 0 % if the ev is returned to a charging station
                    EVs_t_temp(i, 19) = t + 1; % record the time when the ev just arrived the charging station
                end
            end
        elseif EVs_t_temp(i, 8) == 1 % if the ev is available (for reservation or rebalance)
            EVs_t_temp(i, 12) = EVs_t_temp(i, 12) + tdeta;
        end
    end
end

% update info of available ev based on Customer's reservation
% a customer can actively reserve an ev if he/she has not actively reserved one yet
for j = 1 : size(Customers_Unserved_temp, 1)
    if Customers_Unserved_temp(j, 13) == 0
        dl = [];
        c_Zone = Nodes(Customers_Unserved_temp(j, 3), 2);
        c_Zone_index = find(Zones(:, 1) == c_Zone);
        for i = 1 : Nev
            ev_Zone_index = find(Zones(:, 2) == EVs_t_temp_pre(i, 2) & Zones(:, 3) == EVs_t_temp_pre(i, 3));
            % the original CL requirement below no longer needed since we allow vehicles to get charged during driving
            if EVs_t_temp_pre(i, 8) == 1 && Customers_Unserved_temp(j, 8) <= (MaxCL - 1) && EVs_t_temp_pre(i, 4) >= Customers_Unserved_temp(j, 8) %(MaxCL * Customers_Unserved_temp(j, 9) / range + minReturnCL(dest)) 
                dis = DrivingDistanceMatrix(c_Zone_index, ev_Zone_index);
            elseif EVs_t_temp_pre(i, 8) == 1 && Customers_Unserved_temp(j, 8) == MaxCL && ceil(EVs_t_temp_pre(i, 4)) >= Customers_Unserved_temp(j, 8) % the customer with MaxCL requirement is special because such customers arriving at non-charging stations may nerver be served 
                dis = DrivingDistanceMatrix(c_Zone_index, ev_Zone_index);
            else
                dis = 10000;
            end
            dl = [dl; dis];
        end
        if min(dl) <= max_D % if nearest ev within max_D
            min_index = find(dl(:) == min(dl), 1); % customer will choose to reserve the nearest ev
            EVs_t_temp_pre(min_index, 8) = 0;
            EVs_t_temp(min_index, 8) = 0;
            EVs_t_temp(min_index, 11) = Customers_Unserved_temp(j, 1);
            Customers_Unserved_temp(j, 13) = 1;
            Customers_Unserved_temp(j, 15) = t + 1 + round(dl(min_index) / Customers_Unserved_temp(j, 12));
        end
    end
end

EVs_t{t + 1, 1} = EVs_t_temp;
Customers_Unserved{t + 1, 1} = Customers_Unserved_temp;

if t >= 3
   EVs_t{t - 2, 1}=[];
end
    
    