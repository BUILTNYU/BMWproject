figure
time = 1 : T;
subplot(3,2,1);
if size(RebalanceList, 1) > 0
    histogram(RebalanceList(:, 3))
end
xlabel ('Distance (km)','FontSize',11);
ylabel ('Number of rebalances','FontSize',11);
title ('Travel distance for rebalance','FontSize',11);

subplot(3,2,2);
plot(time, TimeSeries(:, 2))
xlim([0 T]);
% ylim([0 5000]);
xlabel ('Time (min)','FontSize',11);
ylabel ('Distance (km)','FontSize',11);
title ('Accumulated Rebalance Distance','FontSize',11);


subplot(3,3,4);
if size(Customers_Lost, 1)> 0
    histogram([Customers_Finished{T, 1}(:, 11); Customers_Lost(:, 11)])
else
    histogram(Customers_Finished{T, 1}(:, 11))
end
% xlim([0 80])
xlabel ('Waiting Time (min)','FontSize',11);
ylabel ('Number of customers','FontSize',11);
title ('Customer waiting time distribution','FontSize',11);

subplot(3,3,5);
plot(time, TimeSeries(:, 3))
xlim([0 T]);
xlabel ('Time (min)','FontSize',11);
ylabel ('Waiting time (min)','FontSize',11);
title ('Accumulated waiting time','FontSize',11);

subplot(3,3,6);
plot(time, TimeSeries(:, 4))
% ylim([0 30])
xlim([0 T]);
xlabel ('Time (min)','FontSize',11);
ylabel ('Waiting time (min)','FontSize',11);
title ('Average waiting time','FontSize',11);

subplot(3,2,5);
Number_of_idle_EVs = EV_statistic(:, 1);
Number_of_occupied_EVs = EV_statistic(:, 2);
Number_of_rebalanced_EVs = EV_statistic(:, 3);
plot(time, Number_of_idle_EVs, 'g', time, Number_of_occupied_EVs, 'r', time, Number_of_rebalanced_EVs, 'b')
xlim([0 T]);
xlabel ('Time (min)','FontSize',11);
ylabel ('Number of vehicles','FontSize',11);
title ('Number of EVs','FontSize',11);
legend('Idle EVs','Occupied EVs', 'Rebalanced EVs');
lgd = legend;
lgd.FontSize = 8;
% lgd.NumColumns = 3;

subplot(3,2,6);
% hold on
% for i = 1 : N
%     plot(time, CT_Number_inS(:, i))
% end
plot(time, TimeSeries(:, 1))
xlim([0 T]);
xlabel ('Time (min)','FontSize',11);
ylabel ('Number of customers','FontSize',11);
title ('Number of waiting customers','FontSize',11);

% legend('Station 1', 'Station 2', 'Station 3', 'Station 4', 'Station 5', 'Station 6', 'All Stations');

