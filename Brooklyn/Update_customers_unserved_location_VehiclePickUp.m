function [Customers_Unserved] = Update_customers_unserved_location_VehiclePickUp(t, tdeta, EVs_t, Customers_Unserved)
EVs_t_temp = EVs_t{t, 1};
Customers_Unserved_temp = Customers_Unserved{t, 1};

for i = 1 : size(Customers_Unserved_temp, 1)
    if Customers_Unserved_temp(i, 13) == 1
        c_x = Customers_Unserved_temp(i, 4);
        c_y = Customers_Unserved_temp(i, 5);
        temp = find(EVs_t_temp(:, 11) == Customers_Unserved_temp(i, 1) & EVs_t_temp(:, 8) == 0 & isnan(EVs_t_temp(:, 5)));
        ev_x = EVs_t_temp(temp, 2);
        ev_y = EVs_t_temp(temp, 3);
        if c_x ~= ev_x | c_y ~= ev_y
            d = sqrt((c_x - ev_x)^2 + (c_y - ev_y)^2);
            if d <= Customers_Unserved_temp(i, 12) * tdeta
                Customers_Unserved_temp(i, 4) = ev_x;
                Customers_Unserved_temp(i, 5) = ev_y;
            else
                Customers_Unserved_temp(i, 4) = c_x + (ev_x - c_x) * Customers_Unserved_temp(i, 12) * tdeta / d;
                Customers_Unserved_temp(i, 5) = c_y + (ev_y - c_y) * Customers_Unserved_temp(i, 12) * tdeta / d;
            end
        end
    end
end

Customers_Unserved{t + 1, 1} = Customers_Unserved_temp;