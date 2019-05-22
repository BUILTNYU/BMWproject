function [EVs_t, RebalanceList] = RebalanceInfoCollection(Nev, t, EVs_t, RebalanceList)

EVs_t_temp_pre = EVs_t{t, 1};
EVs_t_temp = EVs_t{t + 1, 1};

for i = 1 : Nev
    if ~isnan(EVs_t_temp_pre(i, 20)) && isnan(EVs_t_temp(i, 20))
        RebalanceNo = size(RebalanceList, 1) + 1;
        RebalanceStartTime = EVs_t_temp(i, 21);
        RebalanceDis = EVs_t_temp(i, 23);
        RebalanceTotalTime = EVs_t_temp(i, 22) - EVs_t_temp(i, 21);
        ChargeWaitingTime = EVs_t_temp(i, 24);
        RebalanceList = [RebalanceList; RebalanceNo, RebalanceStartTime, RebalanceDis, RebalanceTotalTime, ChargeWaitingTime];
        EVs_t_temp(i, 21) = NaN;
        EVs_t_temp(i, 22) = NaN;
        EVs_t_temp(i, 23) = NaN;
        EVs_t_temp(i, 24) = NaN;
    end
end

EVs_t{t + 1, 1} = EVs_t_temp;