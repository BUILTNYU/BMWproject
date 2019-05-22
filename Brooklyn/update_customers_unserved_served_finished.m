function [Customers_Unserved, Customers_Served, Customers_Finished, Customers_Lost] = update_customers_unserved_served_finished(Nev, t, tolerance, EVs_t, Customers_Unserved, Customers_Served, Customers_Finished, Customers_Lost)
EVs_t_temp = EVs_t{t, 1};
EVs_t_temp_next = EVs_t{t + 1, 1};
Customers_Unserved_temp = Customers_Unserved{t + 1, 1};
Customers_Served_temp = Customers_Served{t, 1};
Customers_Finished_temp = Customers_Finished{t, 1};

number_of_unserved_customers = size(Customers_Unserved_temp, 1);
if number_of_unserved_customers > 0
    for i = 1 : number_of_unserved_customers
        Customers_Unserved_temp(i, 11) = t + 1 - Customers_Unserved_temp(i, 2);
    end
end

for i = 1 : Nev   
    if EVs_t_temp(i, 10) == 0 && EVs_t_temp_next(i, 10) == 1
       temp = find(Customers_Unserved_temp(:, 1) == EVs_t_temp_next(i, 11));
       Customers_Served_temp = [Customers_Served_temp; Customers_Unserved_temp(temp, :), i, t + 1];
       Customers_Unserved_temp(temp, :) = [];
    end
    
    if EVs_t_temp(i, 10) == 1 && EVs_t_temp_next(i, 10) == 0
       temp = find(Customers_Served_temp(:, 1) == EVs_t_temp(i, 11));
       Using_time = t + 1 - Customers_Served_temp(temp, 17);
       Customers_Finished_temp = [Customers_Finished_temp; Customers_Served_temp(temp, :), t + 1, Using_time];
       Customers_Served_temp(temp, :) = [];
    end
end

if number_of_unserved_customers > 0
    abandon_list = find(Customers_Unserved_temp(:, 11) > tolerance & Customers_Unserved_temp(:, 13) == 0);
    if ~isempty(abandon_list)
        for i = 1 : length(abandon_list)
            Customers_Lost = [Customers_Lost; Customers_Unserved_temp(abandon_list(i), :), t + 1];
        end
    end
    Customers_Unserved_temp(abandon_list, :) = [];
end

Customers_Unserved{t + 1, 1} = Customers_Unserved_temp;
Customers_Served{t + 1, 1} = Customers_Served_temp;
Customers_Finished{t + 1, 1} = Customers_Finished_temp;

if t >= 3
    Customers_Unserved{t - 2, 1} = [];
    Customers_Served{t - 2, 1} = [];
    Customers_Finished{t - 2, 1} = [];
end

