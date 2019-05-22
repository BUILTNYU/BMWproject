function [initial_loc_t, station_capacity_t, station_waiting_number_t] = generate_input(Nev, N, t, Zones, Nodes, EVs_t, initial_loc_t, station_capacity_t, station_waiting_number_t)

EVs_t_temp_pre = EVs_t{t, 1};
EVs_t_temp = EVs_t{t + 1, 1};

ev = [];
loc = [];
for i = 1 : Nev
    if EVs_t_temp_pre(i, 8) == 1 &&  EVs_t_temp(i, 8) == 1 % if the ev is available
        zone_no_index = find(Zones(:, 2) == EVs_t_temp(i, 2) & Zones(:, 3) == EVs_t_temp(i, 3));
        zone_no = Zones(zone_no_index, 1);
        node_no = find(Nodes(:, 2) == zone_no & Nodes(:, 3) == ceil(EVs_t_temp(i, 4)));
        ev = [ev; i];
        loc = [loc; node_no];
    end
end
initial_loc = [ev, loc];

zone = [];
chargers = [];
waiting_number = [];
for i = 1 : N
    if Zones(i, 4) > 0
        zone = [zone; i];
        chargers = [chargers; Zones(i, 4) - length(find(EVs_t_temp(:, 9) == 1 & EVs_t_temp(:, 2) == Zones(i, 2) & EVs_t_temp(:, 3) == Zones(i, 3)))];
        waiting_number = [waiting_number; length(find(EVs_t_temp(:, 9) == 0 & EVs_t_temp(:, 2) == Zones(i, 2) & EVs_t_temp(:, 3) == Zones(i, 3) & ~isnan(EVs_t_temp(:, 19))))];
    end
end
station_capacity = [zone, chargers];
station_waiting_number = [zone, waiting_number];

initial_loc_t{t + 1, 1} = initial_loc;
station_capacity_t{t + 1, 1} = station_capacity;
station_waiting_number_t{t + 1, 1} = station_waiting_number;
