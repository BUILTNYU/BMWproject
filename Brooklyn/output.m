function [CN_ave, WT_ave, LostCN, sum_idleTime, sum_occupiedKm] = output(T, CT_Number_inN, Customers_Served, Customers_Finished, Customers_Lost, EVs_t)
LostCN = size(Customers_Lost, 1);
sum_idleTime = sum(EVs_t{size(EVs_t, 1), 1}(:, 12));
sum_occupiedKm = sum(EVs_t{size(EVs_t, 1), 1}(:, 13));

CN_max = [];
CN_ave = [];
for t = 1 : T
    CN_max = max([CN_max,max(CT_Number_inN(t, :))]);
    CN_ave = mean([CN_ave,mean(CT_Number_inN(t, :))]);
end
if size(Customers_Finished{T, 1}, 1) > 0
    WT_max = max(Customers_Finished{T, 1}(:, 11));
    WT_ave = mean(Customers_Finished{T, 1}(:, 11));
    WT_CV = WT_ave / sqrt(var(Customers_Finished{T, 1}(:, 11)));

    UT_max = max(Customers_Finished{T, 1}(:, 19));
    UT_ave = mean(Customers_Finished{T, 1}(:, 19));
    UT_CV = UT_ave / sqrt(var(Customers_Finished{T, 1}(:, 19)));
elseif size(Customers_Served{T, 1}, 1) > 0
    WT_max = max(Customers_Served{T, 1}(:, 11));
    WT_ave = mean(Customers_Served{T, 1}(:, 11));
    WT_CV = WT_ave / sqrt(var(Customers_Served{T, 1}(:, 11)));

    UT_max = NaN;
    UT_ave = NaN;
    UT_CV = NaN;
else
    WT_max = NaN;
    WT_ave = NaN;
    WT_CV = NaN;

    UT_max = NaN;
    UT_ave = NaN;
    UT_CV = NaN;
end
