function [CT_Number_inN, TimeSeries, EV_statistic] = update_evaluations(N, MaxCL, t, Customers_Unserved, Customers_Served, Customers_Finished, Customers_Lost, EVs_t, CT_Number_inN, TimeSeries, EV_statistic, RebalanceList)
Customers_Unserved_temp = Customers_Unserved{t, 1};
Customers_Served_temp = Customers_Served{t, 1};
Customers_Finished_temp = Customers_Finished{t, 1};
if size(Customers_Lost, 1) > 0
    Customers_Lost_temp = Customers_Lost(find(Customers_Lost(:, 16) <= t), :);
else
    Customers_Lost_temp = [];
end
EVs_t_temp = EVs_t{t, 1};

Total_waiting_CT = size(Customers_Unserved_temp, 1);

if size(RebalanceList, 1) > 0
    RebalanceList_temp = RebalanceList(find(RebalanceList(:, 2) + RebalanceList(:, 4) < t + 1), :);
    AccuOpeCost = sum(RebalanceList_temp(:, 3)) + sum(EVs_t_temp(find(~isnan(EVs_t_temp(:, 23))), 23));
else
    AccuOpeCost = sum(EVs_t_temp(find(~isnan(EVs_t_temp(:, 23))), 23));
end

if size(Customers_Unserved_temp, 1) > 0
    a = sum(Customers_Unserved_temp(:, 11));
else
    a = 0;
end
if size(Customers_Served_temp, 1) > 0
    b = sum(Customers_Served_temp(:, 11));
else
    b = 0;
end
if size(Customers_Finished_temp, 1) > 0
    c = sum(Customers_Finished_temp(:, 11));
else
    c = 0;
end
if size(Customers_Lost_temp, 1) > 0
    d = sum(Customers_Lost_temp(:, 11));
else
    d = 0;
end
AccuDelay = a + b + c + d;

if size(Customers_Unserved_temp, 1) + size(Customers_Served_temp, 1) + size(Customers_Finished_temp, 1) + size(Customers_Lost_temp, 1)> 0
    AveWaitingTime_allHistoryCustomer = AccuDelay / (size(Customers_Unserved_temp, 1) + size(Customers_Served_temp, 1) + size(Customers_Finished_temp, 1) + size(Customers_Lost_temp, 1));
else
    AveWaitingTime_allHistoryCustomer = 0;
end
TimeSeries = [TimeSeries; Total_waiting_CT, AccuOpeCost, AccuDelay, AveWaitingTime_allHistoryCustomer];


CT_Number_inN_temp = [];
if size(Customers_Unserved_temp, 1) >= 1
    for i = 1 : N * MaxCL
        CT_Number_inN_temp = [CT_Number_inN_temp, length(find(Customers_Unserved_temp(:, 3) == i))];
    end
else
    for i = 1 : N * MaxCL
        CT_Number_inN_temp = [CT_Number_inN_temp, 0];
    end
end
CT_Number_inN = [CT_Number_inN; CT_Number_inN_temp];

Number_of_idle_EVs = length(find(EVs_t_temp(:, 8) == 1));
Number_of_occupiedAR_EVs = length(find(EVs_t_temp(:, 10) == 1)) + length(find(~isnan(EVs_t_temp(:, 11)) & EVs_t_temp(:, 10) == 0 & isnan(EVs_t_temp(:, 5))));
Number_of_rebalanced_EVs = length(find(~isnan(EVs_t_temp(:, 5))));
EV_statistic = [EV_statistic; Number_of_idle_EVs, Number_of_occupiedAR_EVs, Number_of_rebalanced_EVs];