function [EVs_t] = update_EVs_NP(Nev, MaxCL, v, cr, e, t, Customers_Unserved, EVs_t)
EVs_t_temp = EVs_t{t, 1};
EVs_t_temp_pre = EVs_t{t, 1};
Customers_Unserved_temp = Customers_Unserved{t, 1};

% if C_Node is not NaN, then the EV is assigned to a C_Node, but has not picked up its customer
% Once the EV is assigned to a customer, C_Node becomes NaN, Customer_No becomes non-NaN
% Once the EV picks up a customer, Occu becomes 1
% if C_Node is NaN, and Customer_No is NaN, then EV_Node must be 1, and the EV is available for new assignment
% only when EV_Node is 1, D_X and D_Y is NaN

% update the EVs_t based on itself's previous data
for i = 1 : Nev
    %update Xcoor, Ycoor of EVs with destinations
    if ~isnan(EVs_t_temp_pre(i, 6))
        d = sqrt((EVs_t_temp_pre(i, 2) - EVs_t_temp(i, 6)) ^ 2 + (EVs_t_temp_pre(i, 3) - EVs_t_temp(i, 7)) ^ 2);
        if d <= v
            EVs_t_temp(i, 2) = EVs_t_temp(i, 6);
            EVs_t_temp(i, 3) = EVs_t_temp(i, 7);
        else
            EVs_t_temp(i, 2) = EVs_t_temp_pre(i, 2) + (EVs_t_temp(i, 6) - EVs_t_temp_pre(i, 2)) * v / d;
            EVs_t_temp(i, 3) = EVs_t_temp_pre(i, 3) + (EVs_t_temp(i, 7) - EVs_t_temp_pre(i, 3)) * v / d;
        end
    end
    d_traveled = sqrt((EVs_t_temp_pre(i, 2)-EVs_t_temp(i, 2)) ^ 2 + (EVs_t_temp_pre(i, 3)-EVs_t_temp(i, 3)) ^ 2);
    if EVs_t_temp_pre(i, 10) == 0
        EVs_t_temp(i, 12) = EVs_t_temp(i, 12) + d_traveled;
    else
        EVs_t_temp(i, 13) = EVs_t_temp(i, 13) + d_traveled;
    end
    %update CL
    if EVs_t_temp_pre(i, 9) == 1
        EVs_t_temp(i, 4) = min([EVs_t_temp_pre(i, 4) + cr, MaxCL]);
    else
        EVs_t_temp(i, 4) = EVs_t_temp_pre(i, 4) - d_traveled * e;
    end
    
    %update Occu, Customer_No, Destination, and Incharge when an EV with occupied customer just arrives its destination
    if EVs_t_temp_pre(i, 10) == 1 && EVs_t_temp(i, 2) == EVs_t_temp_pre(i, 6) && EVs_t_temp(i, 3) == EVs_t_temp_pre(i, 7)
       EVs_t_temp(i, 10) = 0;
       EVs_t_temp(i, 11) = NaN;
       EVs_t_temp(i, 6) = NaN;
       EVs_t_temp(i, 7) = NaN;
       EVs_t_temp(i, 9) = 1; % here we assume EV_charging stations everywhere. By adding the GPS for EV_charging stations, we could modify this code
       EVs_t_temp(i, 14) = EVs_t_temp(i, 14) + 1;
    end
    
    %update Occu, Customer_No, Destination, and Incharge when an EV with reserved customer just arrives its customer
    if ~isnan(EVs_t_temp_pre(i, 11)) && EVs_t_temp_pre(i, 10) == 0 && EVs_t_temp(i, 2) == EVs_t_temp_pre(i, 6) && EVs_t_temp(i, 3) == EVs_t_temp_pre(i, 7)
       EVs_t_temp(i, 5) = NaN;
       EVs_t_temp(i, 6) = Customers_Unserved_temp(find(Customers_Unserved_temp(:, 1) == EVs_t_temp_pre(i, 11)), 6);
       EVs_t_temp(i, 7) = Customers_Unserved_temp(find(Customers_Unserved_temp(:, 1) == EVs_t_temp_pre(i, 11)), 7);
       EVs_t_temp(i, 10) = 1;
       % once a customer is served, delete it from the customer_unserved list to prevent one customer from being served by more than one vehicle
       Customers_Unserved_temp(find(Customers_Unserved_temp(:, 1) == EVs_t_temp_pre(i, 11)), :) = [];    
    end
        
    %update EV_Node, an EV is in an EV_Node only when it is available
    if isnan(EVs_t_temp(i, 6)) && EVs_t_temp(i, 8) == 0
       EVs_t_temp(i, 8) = 1; %EV is in an EV node
    end
end

Number_of_available_ev = length(find(EVs_t_temp(:, 8) == 1));
Number_of_waiting_customer = length(find(Customers_Unserved_temp(:, 10) == 0));
if Number_of_available_ev <= Number_of_waiting_customer % if the customer is more than ev, than ev chooses customer
    % update Destination, EV_Node, Incharge, Occu, and Customer_No based on Vehicle's choice
    % a vehicle can actively reserve a customer if he/she has not been actively reserved yet
    for i = 1 : Nev
        if EVs_t_temp(i, 8) == 1
            dl = [];
            ev_x = EVs_t_temp(i, 2);
            ev_y = EVs_t_temp(i, 3);
            for j = 1 : size(Customers_Unserved_temp, 1)
                c_x = Customers_Unserved_temp(j, 4);
                c_y = Customers_Unserved_temp(j, 5);
                emn = sqrt((c_x - ev_x)^2 + (c_y - ev_y)^2) * e;
                if (Customers_Unserved_temp(j, 8) - 1) <= (EVs_t_temp(i, 4) - emn) && Customers_Unserved_temp(j, 10) == 0% if the ev's CL is enough and the customer has not been actively reserved yet
                    dis = sqrt((ev_x - c_x)^2 + (ev_y - c_y)^2);
                else % if the ev's CL is not enough, set the dis to be a large number
                    dis = 1000;
                end
                dl = [dl; dis];
            end
            if min(dl) ~= 1000 % if any customer is waiting, the vehicle will choose the nearest customer to serve
                min_index = find(dl(:) == min(dl), 1);
                EVs_t_temp(i, 6) = Customers_Unserved_temp(min_index, 4);
                EVs_t_temp(i, 7) = Customers_Unserved_temp(min_index, 5);
                EVs_t_temp(i, 8) = 0;
                EVs_t_temp(i, 9) = 0;
                EVs_t_temp(i, 11) = Customers_Unserved_temp(min_index, 1);
                Customers_Unserved_temp(min_index, :) = [];
            end
        end
    end
else % if the ev is more than customer, then customer chooses ev
    % update Destination, EV_Node, Incharge, Occu, and Customer_No based on Customer's choice
    % a customer can actively reserve an ev if he/she has not actively reserved one yet
    for j = 1 : size(Customers_Unserved_temp, 1)
        if Customers_Unserved_temp(j, 10) == 0
            dl = [];
            c_x = Customers_Unserved_temp(j, 4);
            c_y = Customers_Unserved_temp(j, 5);
            for i = 1 : Nev
                ev_x = EVs_t_temp(i, 2);
                ev_y = EVs_t_temp(i, 3);
                emn = sqrt((c_x - ev_x)^2 + (c_y - ev_y)^2) * e;
                % -1 because we + 1 when calculate the CL for customer
                if (Customers_Unserved_temp(j, 8) - 1) <= (EVs_t_temp(i, 4) - emn) && EVs_t_temp(i, 8) == 1 % if the ev's CL is enough and it is available
                    dis = sqrt((ev_x - c_x)^2 + (ev_y - c_y)^2);
                else % if the ev's CL is not enough, or the ev is not available, set the dis to be a large number
                    dis = 1000;
                end
                dl = [dl; dis];
            end
            if min(dl) ~= 1000 % if any ev is available, the customer could consider whether to actively reserve an ev or not
                min_index = find(dl(:) == min(dl), 1); % customer will choose to reserve the nearest ev
                EVs_t_temp(min_index, 6) = Customers_Unserved_temp(j, 4);
                EVs_t_temp(min_index, 7) = Customers_Unserved_temp(j, 5);
                EVs_t_temp(min_index, 8) = 0;
                EVs_t_temp(min_index, 9) = 0;
                EVs_t_temp(min_index, 11) = Customers_Unserved_temp(j, 1);
            end
        end
    end
end

EVs_t{t + 1, 1} = EVs_t_temp;

% if t >= 3
%    EVs_t{t-1,1}=[];
% end
    
    